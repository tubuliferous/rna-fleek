"""RNA-FLEEK supervisor — multi-user front for fleek subprocesses.

Listens on a public port, serves a login + signup UI, manages a SQLite
user database, and reverse-proxies every authenticated request to a
per-user fleek subprocess that's spawned on first use and reaped when
idle. Each subprocess runs `fleek --remote-mode` with a `--data-dir`
unique to that user; the user's Anthropic API key (if they supplied
one at signup) is passed via env var. The server-side endpoint
lockdowns the fleek subprocess enforces (chrooted browse / load,
cache + key writes refused, etc.) are the actual security boundary —
this module is the gate that decides whose --data-dir to use and
whose key to inject.

Architecture (single Python process, threading-based):

    [browser] ─https→ [nginx :443]  ──┐
                  ↓ optional lab gate │ supervisor :8080
                                      │   ├─ /login, /signup, /logout
                                      │   ├─ cookie validation (HMAC-SHA256)
                                      │   ├─ proxies → subprocess port
                                      │   └─ idle-reap thread
                                      ↓
                                 ┌─────────────────────┐
                                 │ fleek (alice) :19001│
                                 │ fleek (bob)   :19002│
                                 │   ...               │
                                 └─────────────────────┘

Per CLAUDE.md spirit: stdlib-only (no Flask / FastAPI / bcrypt / etc.).
PBKDF2-HMAC-SHA256 for password hashing; HMAC-signed cookies; threading
+ http.client for the proxy. Designed for ~5–50 users on a single VM,
not as a JupyterHub replacement."""

import argparse
import hashlib
import hmac
import http.client
import json
import os
import secrets
import socket
import sqlite3
import subprocess
import sys
import threading
import time
import urllib.parse
from http.server import BaseHTTPRequestHandler, HTTPServer
from socketserver import ThreadingMixIn
from pathlib import Path


# ── Globals (set at startup from CLI) ────────────────────────────────
USER_BASE = None         # Path: per-user dirs live under this
SHARED_DATA = None       # Path: read-only shared dataset dir
USER_QUOTA_MB = 100      # Per-user write quota
DB_PATH = None           # Path to SQLite users.db
COOKIE_SECRET = None     # bytes: HMAC-signing key for session cookies
PORT_RANGE = (19000, 19999)
IDLE_TIMEOUT_S = 30 * 60         # Reap subprocesses idle longer than this
SESSION_LIFETIME_S = 24 * 3600   # Cookie validity
FLEEK_PYTHON = sys.executable
FLEEK_SERVER_MODULE = "rna_fleek.server"  # python -m <this>

# Subprocess registry — username → {"process": Popen, "port": int,
# "last_activity": float}. Mutated under REGISTRY_LOCK; the reaper
# thread and incoming requests both touch it.
PROCESS_REGISTRY = {}
REGISTRY_LOCK = threading.Lock()

# Reserved usernames that signup can't claim. "admin" so a stranger
# past the lab gate can't squat the obvious name; the rest are
# operational identifiers we don't want to collide with.
RESERVED_USERNAMES = {"admin", "root", "fleek", "shared", "system"}


# ── Database ─────────────────────────────────────────────────────────
def init_db():
    """Create the users table on first run. Re-running is a no-op."""
    conn = sqlite3.connect(str(DB_PATH))
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS users (
            username       TEXT PRIMARY KEY,
            password_hash  TEXT NOT NULL,
            created_at     REAL NOT NULL,
            last_login     REAL,
            api_key        TEXT,
            assigned_port  INTEGER UNIQUE
        );
    """)
    conn.commit()
    conn.close()


def hash_password(pw, salt=None):
    """PBKDF2-HMAC-SHA256, 200k iterations. Returns 'salt_hex$hash_hex'.
    bcrypt would be the textbook choice; PBKDF2 keeps us stdlib-only and
    is the OWASP-recommended fallback."""
    if salt is None:
        salt = secrets.token_bytes(16)
    elif isinstance(salt, str):
        salt = bytes.fromhex(salt)
    h = hashlib.pbkdf2_hmac("sha256", pw.encode("utf-8"), salt, 200_000)
    return f"{salt.hex()}${h.hex()}"


def verify_password(pw, stored):
    """Constant-time compare against the stored hash."""
    try:
        salt_hex, _ = stored.split("$", 1)
    except ValueError:
        return False
    return hmac.compare_digest(stored, hash_password(pw, salt_hex))


def db_get_user(username):
    conn = sqlite3.connect(str(DB_PATH))
    cur = conn.execute(
        "SELECT username, password_hash, api_key, assigned_port FROM users WHERE username = ?",
        (username,),
    )
    row = cur.fetchone()
    conn.close()
    if not row:
        return None
    return {
        "username": row[0],
        "password_hash": row[1],
        "api_key": row[2],
        "assigned_port": row[3],
    }


def db_create_user(username, pw_hash, api_key=None):
    """Insert a new user, allocate a port for their fleek subprocess.
    Raises sqlite3.IntegrityError on username collision (caller turns
    that into a friendly UI error)."""
    conn = sqlite3.connect(str(DB_PATH))
    try:
        port = _next_free_port_in_db(conn)
        conn.execute(
            "INSERT INTO users (username, password_hash, created_at, api_key, assigned_port) VALUES (?, ?, ?, ?, ?)",
            (username, pw_hash, time.time(), api_key, port),
        )
        conn.commit()
    finally:
        conn.close()
    return port


def db_record_login(username):
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("UPDATE users SET last_login = ? WHERE username = ?", (time.time(), username))
    conn.commit()
    conn.close()


def db_update_api_key(username, api_key):
    """Set or clear (api_key=None) the user's stored Anthropic key.
    Caller is responsible for terminating the user's running fleek
    subprocess afterward so it respawns with the new ANTHROPIC_API_KEY
    in env."""
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("UPDATE users SET api_key = ? WHERE username = ?", (api_key, username))
    conn.commit()
    conn.close()


def _next_free_port_in_db(conn):
    """Lowest unallocated port in PORT_RANGE. We trust the DB as the
    source of truth — once allocated to a user, that user keeps that
    port forever (so log lines / nginx routes are stable across
    restarts). With a 1000-port range we'll never run out for lab
    use; if we do, expand the range."""
    used = {r[0] for r in conn.execute("SELECT assigned_port FROM users WHERE assigned_port IS NOT NULL").fetchall()}
    for p in range(PORT_RANGE[0], PORT_RANGE[1] + 1):
        if p not in used:
            return p
    raise RuntimeError(f"No free ports in range {PORT_RANGE}")


# ── Cookie sessions ──────────────────────────────────────────────────
def make_session_cookie(username):
    """Returns a cookie value carrying {u: username, exp: ts} signed with
    HMAC-SHA256. Format: 'v1.<payload-hex>.<sig-hex>'. Embedding the
    payload in the cookie (rather than a server-side session table)
    keeps the supervisor stateless and survives restarts cleanly."""
    payload = json.dumps(
        {"u": username, "exp": time.time() + SESSION_LIFETIME_S},
        separators=(",", ":"),
    ).encode("utf-8")
    sig = hmac.new(COOKIE_SECRET, payload, hashlib.sha256).digest()
    return f"v1.{payload.hex()}.{sig.hex()}"


def verify_session_cookie(cookie):
    """Return the username if cookie is valid + unexpired, else None."""
    try:
        if not cookie.startswith("v1."):
            return None
        _, payload_hex, sig_hex = cookie.split(".", 2)
        payload = bytes.fromhex(payload_hex)
        sig = bytes.fromhex(sig_hex)
        expected = hmac.new(COOKIE_SECRET, payload, hashlib.sha256).digest()
        if not hmac.compare_digest(sig, expected):
            return None
        data = json.loads(payload.decode("utf-8"))
        if data.get("exp", 0) < time.time():
            return None
        return data.get("u")
    except Exception:
        return None


def _load_cookie_secret():
    """Read or generate the cookie-signing secret. Lives in
    USER_BASE/.fleek_supervisor_secret (mode 600). Generated once on
    first boot, persists across restarts so existing cookies stay
    valid; rotating it is a deliberate action (delete the file)."""
    sec_path = USER_BASE / ".fleek_supervisor_secret"
    if sec_path.exists():
        return sec_path.read_bytes()
    sec = secrets.token_bytes(32)
    sec_path.write_bytes(sec)
    try:
        sec_path.chmod(0o600)
    except OSError:
        pass
    return sec


# ── Subprocess management ────────────────────────────────────────────
def ensure_user_process(username):
    """Return the port the user's fleek subprocess is bound to. Spawns
    one if not running. Idempotent — subsequent calls just bump the
    last-activity timestamp. The reaper thread reclaims idle slots."""
    with REGISTRY_LOCK:
        entry = PROCESS_REGISTRY.get(username)
        if entry and entry["process"].poll() is None:
            entry["last_activity"] = time.time()
            return entry["port"]
        # Need to spawn (either no entry or the existing process died)
        user = db_get_user(username)
        if not user:
            raise RuntimeError(f"Unknown user: {username}")
        port = user["assigned_port"]
        user_dir = USER_BASE / username
        user_dir.mkdir(parents=True, exist_ok=True)
        # Lock down to mode 700 so even if the supervisor and fleek run
        # under the same OS uid, ls /data/users/ doesn't leak content
        # to other users on the box. Fails silently on filesystems
        # that don't support unix perms (e.g. some bind mounts).
        try:
            user_dir.chmod(0o700)
        except OSError:
            pass
        env = os.environ.copy()
        if user["api_key"]:
            env["ANTHROPIC_API_KEY"] = user["api_key"]
        else:
            # Don't inherit a server-level key when the user hasn't set
            # their own — they should see "no key" in their fleek UI,
            # not silently use the admin's key (which would mix bills
            # and break per-user attribution if we ever added that).
            env.pop("ANTHROPIC_API_KEY", None)
        cmd = [
            FLEEK_PYTHON, "-m", FLEEK_SERVER_MODULE,
            "--remote-mode",
            "--data-dir", str(user_dir),
            "--shared-data", str(SHARED_DATA),
            "--user-quota-mb", str(USER_QUOTA_MB),
            "--port", str(port),
            "--host", "127.0.0.1",
            "--no-browser",
        ]
        proc = subprocess.Popen(
            cmd, env=env,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        # Wait for the subprocess to bind the port. Polls every 100ms
        # for up to 15s — covers cold-imports of scanpy / anndata which
        # can take 5–10s on a server with disk-cached wheels.
        deadline = time.time() + 15
        bound = False
        while time.time() < deadline:
            if proc.poll() is not None:
                raise RuntimeError(
                    f"fleek subprocess for '{username}' exited prematurely (rc={proc.returncode})"
                )
            try:
                with socket.create_connection(("127.0.0.1", port), timeout=0.5):
                    bound = True
                    break
            except (OSError, socket.timeout):
                time.sleep(0.1)
        if not bound:
            try:
                proc.kill()
            except Exception:
                pass
            raise RuntimeError(f"fleek subprocess for '{username}' did not bind port {port} in 15s")
        PROCESS_REGISTRY[username] = {
            "process": proc, "port": port, "last_activity": time.time(),
        }
        print(f"  [supervisor] spawned fleek for {username} on port {port} (pid={proc.pid})")
        return port


def reap_idle_processes():
    """Background thread: every minute, kill subprocesses idle longer
    than IDLE_TIMEOUT_S, and clean up entries for processes that died
    on their own (e.g. fleek's own auto-quit heartbeat). Also handles
    the case where a subprocess crashed — its slot gets removed so
    the next request triggers a fresh spawn."""
    while True:
        time.sleep(60)
        now = time.time()
        with REGISTRY_LOCK:
            stale = []
            for username, entry in PROCESS_REGISTRY.items():
                proc = entry["process"]
                if proc.poll() is not None:
                    print(f"  [supervisor] {username}'s fleek exited (rc={proc.returncode})")
                    stale.append(username)
                elif now - entry["last_activity"] > IDLE_TIMEOUT_S:
                    print(f"  [supervisor] reaping idle fleek for {username}")
                    try:
                        proc.terminate()
                        proc.wait(timeout=5)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
                    stale.append(username)
            for u in stale:
                PROCESS_REGISTRY.pop(u, None)


def shutdown_all_processes():
    """Called on supervisor exit. SIGTERM each subprocess so they get a
    chance to flush state cleanly, then KILL after a short timeout."""
    with REGISTRY_LOCK:
        for username, entry in PROCESS_REGISTRY.items():
            try:
                entry["process"].terminate()
            except Exception:
                pass
        for username, entry in list(PROCESS_REGISTRY.items()):
            try:
                entry["process"].wait(timeout=3)
            except Exception:
                try:
                    entry["process"].kill()
                except Exception:
                    pass
        PROCESS_REGISTRY.clear()


# ── HTTP handler ─────────────────────────────────────────────────────
HOP_HEADERS = {
    "connection", "keep-alive", "proxy-authenticate", "proxy-authorization",
    "te", "trailers", "transfer-encoding", "upgrade", "host",
}


class SupervisorHandler(BaseHTTPRequestHandler):
    """Routes:
       GET  /login        → login form
       GET  /signup       → signup form
       POST /auth/login   → validate creds, set cookie, 303 to /
       POST /auth/signup  → create user, set cookie, 303 to /
       POST /auth/logout  → clear cookie, 303 to /login
       *  any other path  → reverse-proxy to the user's fleek port
    """

    def do_GET(self):    self._dispatch()
    def do_POST(self):   self._dispatch()
    def do_PUT(self):    self._dispatch()
    def do_DELETE(self): self._dispatch()

    def _dispatch(self):
        path = self.path.split("?", 1)[0]
        try:
            if path == "/login":
                return self._render_form(LOGIN_HTML, error=None)
            if path == "/signup":
                return self._render_form(SIGNUP_HTML, error=None)
            if path == "/auth/login":
                return self._handle_login()
            if path == "/auth/signup":
                return self._handle_signup()
            if path == "/auth/logout":
                return self._handle_logout()
            if path == "/account":
                return self._render_account(error=None, notice=None)
            if path == "/account/api-key":
                return self._handle_account_key()
            # Authenticated routes — proxy to user's fleek
            username = self._get_authenticated_user()
            if not username:
                # API requests get 401, navigations get a redirect.
                if path.startswith("/api/"):
                    return self._send_simple(401, "application/json", b'{"error":"not authenticated"}')
                self.send_response(303)
                self.send_header("Location", "/login")
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            self._proxy_to_user(username)
        except Exception as e:
            try:
                self.send_error(500, f"Supervisor error: {e}")
            except Exception:
                pass

    # ── Auth ────────────────────────────────────────────────────────
    def _get_authenticated_user(self):
        cookie_header = self.headers.get("Cookie", "")
        for c in cookie_header.split(";"):
            c = c.strip()
            if c.startswith("fleek_session="):
                return verify_session_cookie(c[len("fleek_session="):])
        return None

    def _read_form_body(self):
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length).decode("utf-8")
        return urllib.parse.parse_qs(body)

    def _handle_login(self):
        form = self._read_form_body()
        username = (form.get("username", [""])[0] or "").strip().lower()
        password = form.get("password", [""])[0] or ""
        user = db_get_user(username) if username else None
        if not user or not verify_password(password, user["password_hash"]):
            return self._render_form(LOGIN_HTML, error="Invalid username or password.")
        db_record_login(username)
        cookie = make_session_cookie(username)
        self.send_response(303)
        self.send_header("Set-Cookie", f"fleek_session={cookie}; Path=/; HttpOnly; SameSite=Lax; Max-Age={SESSION_LIFETIME_S}")
        self.send_header("Location", "/")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _handle_signup(self):
        form = self._read_form_body()
        username = (form.get("username", [""])[0] or "").strip().lower()
        password = form.get("password", [""])[0] or ""
        password2 = form.get("password2", [""])[0] or ""
        api_key = (form.get("api_key", [""])[0] or "").strip()

        # Validate inputs. Username rules deliberately strict — only
        # alnum + underscore + dash, 3–32 chars. Anything else risks
        # path-traversal in the data dir or weird shell behaviour.
        if not username:
            return self._render_form(SIGNUP_HTML, error="Username is required.")
        if not all(c.isalnum() or c in "_-" for c in username):
            return self._render_form(SIGNUP_HTML, error="Username may only contain letters, digits, underscore, and dash.")
        if not (3 <= len(username) <= 32):
            return self._render_form(SIGNUP_HTML, error="Username must be 3–32 characters.")
        if username in RESERVED_USERNAMES:
            return self._render_form(SIGNUP_HTML, error="That username is reserved.")
        if len(password) < 8:
            return self._render_form(SIGNUP_HTML, error="Password must be at least 8 characters.")
        if password != password2:
            return self._render_form(SIGNUP_HTML, error="Passwords don't match.")
        if api_key and not api_key.startswith("sk-ant-"):
            return self._render_form(SIGNUP_HTML, error="Anthropic API key should start with sk-ant-… (or leave blank).")
        if db_get_user(username):
            return self._render_form(SIGNUP_HTML, error="That username is already taken.")
        # Create the user. The user dir is created lazily on first
        # subprocess spawn — no need to mkdir here.
        try:
            db_create_user(username, hash_password(password), api_key=(api_key or None))
        except sqlite3.IntegrityError:
            return self._render_form(SIGNUP_HTML, error="That username is already taken.")
        cookie = make_session_cookie(username)
        print(f"  [supervisor] new account: {username}")
        self.send_response(303)
        self.send_header("Set-Cookie", f"fleek_session={cookie}; Path=/; HttpOnly; SameSite=Lax; Max-Age={SESSION_LIFETIME_S}")
        self.send_header("Location", "/")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _render_account(self, error=None, notice=None):
        """Account page: shows the username + lets the user update or
        clear their Anthropic API key + has a logout button. Auth
        required; unauthenticated requests redirect to /login."""
        username = self._get_authenticated_user()
        if not username:
            self.send_response(303)
            self.send_header("Location", "/login")
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        user = db_get_user(username)
        has_key = bool(user and user.get("api_key"))
        html = ACCOUNT_HTML
        html = html.replace("{{USERNAME}}", _html_escape(username))
        html = html.replace("{{KEY_STATUS}}",
            "Key on file (last 4: …" + _html_escape(user["api_key"][-4:]) + ")"
            if has_key else "No key on file.")
        if error:
            html = html.replace("{{ERROR}}", f'<div class="err">{_html_escape(error)}</div>')
            html = html.replace("{{NOTICE}}", "")
        elif notice:
            html = html.replace("{{NOTICE}}", f'<div class="ok">{_html_escape(notice)}</div>')
            html = html.replace("{{ERROR}}", "")
        else:
            html = html.replace("{{ERROR}}", "").replace("{{NOTICE}}", "")
        body = html.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _handle_account_key(self):
        """POST /account/api-key — update or clear the user's Anthropic
        API key. Action is one of 'set' or 'clear'. After updating, the
        user's running fleek subprocess (if any) is terminated so the
        next request respawns with the new key in env. The user is
        redirected back to /account with a notice."""
        username = self._get_authenticated_user()
        if not username:
            self.send_response(303)
            self.send_header("Location", "/login")
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        form = self._read_form_body()
        action = (form.get("action", [""])[0] or "").strip()
        if action == "clear":
            db_update_api_key(username, None)
            notice = "API key cleared."
        else:  # 'set'
            new_key = (form.get("api_key", [""])[0] or "").strip()
            if not new_key:
                return self._render_account(error="Key is empty (use Clear to remove).")
            if not new_key.startswith("sk-ant-"):
                return self._render_account(error="Anthropic keys start with sk-ant-…")
            db_update_api_key(username, new_key)
            notice = "API key updated."
        # Kick the user's current fleek subprocess so next request
        # respawns with the new key in ANTHROPIC_API_KEY env.
        with REGISTRY_LOCK:
            entry = PROCESS_REGISTRY.pop(username, None)
        if entry:
            try:
                entry["process"].terminate()
            except Exception:
                pass
        self._render_account(notice=notice)

    def _handle_logout(self):
        # Also tear down the user's fleek subprocess so the next login
        # gets a fresh state (matches user expectation of "log out =
        # forget where I was"). The subprocess would idle-reap anyway,
        # but explicit logout is faster.
        username = self._get_authenticated_user()
        if username:
            with REGISTRY_LOCK:
                entry = PROCESS_REGISTRY.pop(username, None)
            if entry:
                try:
                    entry["process"].terminate()
                except Exception:
                    pass
        self.send_response(303)
        self.send_header("Set-Cookie", "fleek_session=; Path=/; Max-Age=0")
        self.send_header("Location", "/login")
        self.send_header("Content-Length", "0")
        self.end_headers()

    # ── Proxy ───────────────────────────────────────────────────────
    def _proxy_to_user(self, username):
        port = ensure_user_process(username)
        # Read the entire request body before opening the upstream
        # connection. http.client buffers it anyway; this keeps the
        # control flow linear. For uploads (large bodies) we should
        # later switch to a streaming proxy, but for Phase 2 this is
        # fine — uploads in remote mode are bounded by quota.
        body = b""
        if self.headers.get("Content-Length"):
            length = int(self.headers["Content-Length"])
            body = self.rfile.read(length)
        # Filter hop-by-hop headers per RFC 7230 §6.1.
        out_headers = {}
        for k, v in self.headers.items():
            if k.lower() in HOP_HEADERS:
                continue
            out_headers[k] = v
        # Tell the upstream who the user is — useful for logging on
        # the fleek side and a hook for future per-user customizations.
        out_headers["X-Fleek-User"] = username
        try:
            conn = http.client.HTTPConnection("127.0.0.1", port, timeout=300)
            conn.request(self.command, self.path, body=body, headers=out_headers)
            resp = conn.getresponse()
        except (ConnectionRefusedError, http.client.HTTPException, OSError) as e:
            # The user's fleek may have crashed between ensure() and
            # request(). Drop the dead entry so the next call respawns.
            with REGISTRY_LOCK:
                PROCESS_REGISTRY.pop(username, None)
            self.send_error(502, f"Upstream fleek unreachable: {e}")
            return
        try:
            self.send_response(resp.status)
            for h, v in resp.getheaders():
                if h.lower() in HOP_HEADERS:
                    continue
                self.send_header(h, v)
            self.end_headers()
            while True:
                chunk = resp.read(65536)
                if not chunk:
                    break
                try:
                    self.wfile.write(chunk)
                except (BrokenPipeError, ConnectionResetError):
                    break
        finally:
            conn.close()
        with REGISTRY_LOCK:
            if username in PROCESS_REGISTRY:
                PROCESS_REGISTRY[username]["last_activity"] = time.time()

    # ── Helpers ─────────────────────────────────────────────────────
    def _render_form(self, html_template, error=None):
        if error:
            html = html_template.replace("{{ERROR}}", f'<div class="err">{_html_escape(error)}</div>')
        else:
            html = html_template.replace("{{ERROR}}", "")
        body = html.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _send_simple(self, status, ctype, body):
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, format, *args):
        # Silence the default per-request access log; re-enable for
        # debugging by removing this override.
        pass


def _html_escape(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;").replace("'", "&#39;"))


# ── HTML templates ───────────────────────────────────────────────────
# Inline + minimal. Match FLEEK's dark / accent2 palette so the login
# page doesn't feel like a different product. No JS required — just
# plain forms that POST to /auth/{login,signup}.

_BASE_CSS = """
  *{box-sizing:border-box;margin:0;padding:0;}
  body{background:#08080e;color:#dcdee5;font-family:'JetBrains Mono','SF Mono',monospace;
       min-height:100vh;display:flex;align-items:center;justify-content:center;padding:20px;}
  .card{background:#0c0c14;border:1px solid #2a2e3a;border-radius:10px;padding:32px 28px;
        width:100%;max-width:380px;box-shadow:0 24px 80px rgba(0,0,0,0.5);}
  h1{font-size:18px;font-weight:700;margin-bottom:4px;}
  h1 .accent{color:#a5b4fc;letter-spacing:-0.5px;}
  .sub{font-size:11px;color:#7a7e8c;margin-bottom:20px;}
  label{display:block;font-size:10px;color:#7a7e8c;text-transform:uppercase;
        letter-spacing:0.05em;margin-bottom:4px;margin-top:14px;}
  input{width:100%;background:#10131c;color:#dcdee5;border:1px solid #2a2e3a;
        border-radius:4px;font-family:inherit;font-size:12px;padding:8px 10px;outline:none;
        transition:border-color 0.1s;}
  input:focus{border-color:#6366f1;}
  .hint{font-size:10px;color:#7a7e8c;margin-top:3px;}
  button{width:100%;background:#6366f1;color:#fff;border:none;border-radius:4px;
         font-family:inherit;font-size:12px;font-weight:600;padding:10px;cursor:pointer;
         margin-top:20px;transition:filter 0.1s;}
  button:hover{filter:brightness(1.08);}
  .alt{font-size:11px;color:#7a7e8c;text-align:center;margin-top:18px;}
  .alt a{color:#a5b4fc;text-decoration:none;}
  .alt a:hover{text-decoration:underline;}
  .err{background:rgba(220,80,80,0.15);border:1px solid rgba(220,80,80,0.4);
       border-radius:4px;padding:8px 10px;font-size:11px;color:#f0a0a0;margin-bottom:14px;}
"""

LOGIN_HTML = """<!doctype html>
<html><head><meta charset="utf-8"><title>RNA-FLEEK · login</title>
<style>""" + _BASE_CSS + """</style></head>
<body><div class="card">
<h1><span style="color:#7a7e8c;">RNA</span><span class="accent">-FLEEK</span></h1>
<div class="sub">Lab server · sign in to continue</div>
{{ERROR}}
<form method="POST" action="/auth/login">
  <label>Username</label>
  <input name="username" autocomplete="username" autofocus required>
  <label>Password</label>
  <input name="password" type="password" autocomplete="current-password" required>
  <button type="submit">Sign in</button>
</form>
<div class="alt">No account yet? <a href="/signup">Create one</a>.</div>
</div></body></html>"""

ACCOUNT_HTML = """<!doctype html>
<html><head><meta charset="utf-8"><title>RNA-FLEEK · account</title>
<style>""" + _BASE_CSS + """
  .ok{background:rgba(80,180,120,0.15);border:1px solid rgba(80,180,120,0.4);
      border-radius:4px;padding:8px 10px;font-size:11px;color:#90d4a8;margin-bottom:14px;}
  .row{display:flex;gap:8px;align-items:center;margin-top:18px;}
  .row form{flex:1;margin:0;}
  .row button{margin-top:0;}
  .btn-secondary{background:#2a2e3a;color:#dcdee5;}
  .btn-danger{background:#8a3a3a;color:#fff;}
  .stat{font-size:11px;color:#7a7e8c;background:#10131c;border:1px solid #2a2e3a;
        border-radius:4px;padding:6px 10px;margin-top:6px;}
  a.back{color:#a5b4fc;font-size:11px;text-decoration:none;}
  a.back:hover{text-decoration:underline;}
</style></head>
<body><div class="card">
<h1><span style="color:#7a7e8c;">RNA</span><span class="accent">-FLEEK</span></h1>
<div class="sub">Logged in as <strong style="color:#dcdee5;">{{USERNAME}}</strong></div>
{{NOTICE}}{{ERROR}}
<label>Anthropic API key</label>
<div class="stat">{{KEY_STATUS}}</div>
<form method="POST" action="/account/api-key">
  <input type="hidden" name="action" value="set">
  <input name="api_key" type="password" placeholder="sk-ant-..." autocomplete="off" style="margin-top:8px;">
  <div class="hint">Updating restarts your FLEEK session so the new key takes effect.</div>
  <button type="submit">Update key</button>
</form>
<form method="POST" action="/account/api-key" style="margin-top:6px;">
  <input type="hidden" name="action" value="clear">
  <button type="submit" class="btn-secondary">Clear key</button>
</form>
<form method="POST" action="/auth/logout" style="margin-top:24px;">
  <button type="submit" class="btn-danger">Sign out</button>
</form>
<div class="alt"><a href="/" class="back">← Back to FLEEK</a></div>
</div></body></html>"""

SIGNUP_HTML = """<!doctype html>
<html><head><meta charset="utf-8"><title>RNA-FLEEK · sign up</title>
<style>""" + _BASE_CSS + """</style></head>
<body><div class="card">
<h1><span style="color:#7a7e8c;">RNA</span><span class="accent">-FLEEK</span></h1>
<div class="sub">Create your lab account</div>
{{ERROR}}
<form method="POST" action="/auth/signup">
  <label>Username</label>
  <input name="username" autocomplete="username" autofocus required>
  <div class="hint">3–32 chars · letters, digits, _ or - only.</div>
  <label>Password</label>
  <input name="password" type="password" autocomplete="new-password" required>
  <div class="hint">8+ characters.</div>
  <label>Confirm password</label>
  <input name="password2" type="password" autocomplete="new-password" required>
  <label>Anthropic API key <span style="text-transform:none;color:#5a5e6a;">— optional</span></label>
  <input name="api_key" type="password" placeholder="sk-ant-..." autocomplete="off">
  <div class="hint">Required for cell-type AI annotation features. Add later if you don't have one yet.</div>
  <button type="submit">Create account</button>
</form>
<div class="alt">Already have an account? <a href="/login">Sign in</a>.</div>
</div></body></html>"""


# ── CLI ──────────────────────────────────────────────────────────────
class ThreadingHTTPServer(ThreadingMixIn, HTTPServer):
    daemon_threads = True
    allow_reuse_address = True


def main():
    parser = argparse.ArgumentParser(description="RNA-FLEEK supervisor — multi-user front")
    parser.add_argument("--port", type=int, default=8080,
                        help="Port the supervisor binds to (the public port). Default 8080.")
    parser.add_argument("--host", type=str, default="127.0.0.1",
                        help="Host to bind. Default 127.0.0.1; nginx fronts this.")
    parser.add_argument("--user-base", type=str, required=True,
                        help="Per-user data root (one subdir per user).")
    parser.add_argument("--shared-data", type=str, required=True,
                        help="Read-only shared dataset root.")
    parser.add_argument("--user-quota-mb", type=int, default=100,
                        help="Per-user write quota in MB. Default 100.")
    parser.add_argument("--port-range", type=str, default="19000-19999",
                        help="Port range for spawned fleek subprocesses (e.g. 19000-19999).")
    parser.add_argument("--idle-timeout-s", type=int, default=30 * 60,
                        help="Reap subprocesses idle longer than this (seconds). Default 1800.")
    parser.add_argument("--fleek-python", type=str, default=None,
                        help="Python interpreter to run fleek subprocesses with (default: same as supervisor).")
    args = parser.parse_args()

    global USER_BASE, SHARED_DATA, USER_QUOTA_MB, DB_PATH, COOKIE_SECRET
    global PORT_RANGE, IDLE_TIMEOUT_S, FLEEK_PYTHON
    USER_BASE = Path(args.user_base).expanduser().resolve()
    SHARED_DATA = Path(args.shared_data).expanduser().resolve()
    USER_QUOTA_MB = max(1, int(args.user_quota_mb))
    IDLE_TIMEOUT_S = max(60, int(args.idle_timeout_s))
    if args.fleek_python:
        FLEEK_PYTHON = args.fleek_python

    if not SHARED_DATA.exists():
        parser.error(f"--shared-data path does not exist: {SHARED_DATA}")
    USER_BASE.mkdir(parents=True, exist_ok=True)
    try:
        USER_BASE.chmod(0o755)  # supervisor needs +rw on the parent
    except OSError:
        pass

    try:
        lo, hi = (int(x) for x in args.port_range.split("-"))
        if lo < 1024 or hi > 65535 or lo > hi:
            raise ValueError
        PORT_RANGE = (lo, hi)
    except ValueError:
        parser.error(f"--port-range must be 'LOW-HIGH' (e.g. 19000-19999); got {args.port_range!r}")

    DB_PATH = USER_BASE / ".fleek_supervisor.db"
    init_db()
    COOKIE_SECRET = _load_cookie_secret()

    # Background reaper
    threading.Thread(target=reap_idle_processes, daemon=True).start()

    # Graceful shutdown
    import atexit
    atexit.register(shutdown_all_processes)

    server = ThreadingHTTPServer((args.host, args.port), SupervisorHandler)
    print("=" * 50)
    print(f"  RNA-FLEEK supervisor running at: http://{args.host}:{args.port}")
    print(f"    user data: {USER_BASE}")
    print(f"    shared:    {SHARED_DATA}  (read-only)")
    print(f"    quota:     {USER_QUOTA_MB} MB / user")
    print(f"    ports:     {PORT_RANGE[0]}–{PORT_RANGE[1]}")
    print(f"    idle:      reap after {IDLE_TIMEOUT_S}s")
    print("=" * 50)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n  [supervisor] shutting down")
    finally:
        shutdown_all_processes()


if __name__ == "__main__":
    main()
