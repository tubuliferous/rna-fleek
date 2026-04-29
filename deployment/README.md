# RNA-FLEEK · multi-user deployment

A walk-through for putting FLEEK on a shared lab server (e.g. an OpenStack VM
on UAB Cloud RC) with HTTPS, a lab-level password gate, and per-user
accounts. The architecture is described in the project root `CLAUDE.md`;
this doc is the runbook.

```
[browser] ─https→ [nginx :443]      ──┐
                  ↓ basic-auth gate   │ supervisor :8080
                                      │   ├─ /login, /signup, /account, /auth/…
                                      │   └─ proxies → per-user fleek subprocess
                                      ↓
                                 fleek (alice) :19000
                                 fleek (bob)   :19001
                                 …
```

## Layout

```
/data/
  shared/          mode 755, owner=admin, contents readable by fleek user
    datasets/      pre-baked h5ad + caches the whole lab uses
    gmt/           custom GMT files for pathway ORA
  users/           mode 755, owner=fleek
    alice/         mode 700, contains alice's caches/sessions/uploads/exports
    bob/           ...
    .fleek_supervisor.db          SQLite user database
    .fleek_supervisor_secret      auto-generated cookie-signing secret (mode 600)
```

## One-time setup

### 1. Install the package

On the server:

```bash
# Either install from PyPI (when published):
sudo pip install rna-fleek

# Or install from a checkout:
sudo pip install /path/to/RNA-FLEEK
```

That installs both `fleek` and `fleek-supervisor` console scripts.

### 2. Create the system user + dirs

```bash
sudo useradd -r -m -d /home/fleek -s /bin/bash fleek
sudo mkdir -p /data/shared/datasets /data/shared/gmt /data/users
sudo chown -R fleek:fleek /data/users
sudo chmod 755 /data/users

# Pre-bake whichever shared atlases your lab uses. /data/shared
# stays admin-owned (you keep write; fleek can only read), so users
# can't clobber the canonical caches.
sudo chown -R "$USER":"$USER" /data/shared
sudo chmod -R 755 /data/shared
```

### 3. Drop in the systemd unit

```bash
sudo cp deployment/fleek-supervisor.service /etc/systemd/system/
# Edit /etc/systemd/system/fleek-supervisor.service to match your
# paths (User=, ExecStart paths, ports). Then:
sudo systemctl daemon-reload
sudo systemctl enable --now fleek-supervisor.service
sudo journalctl -u fleek-supervisor -f
```

The supervisor binds 127.0.0.1:8080. It's not yet reachable from the
outside.

### 4. Drop in the nginx config

```bash
sudo apt install nginx apache2-utils

# Lab-level password gate (one shared password for the whole lab):
sudo htpasswd -c /etc/nginx/.fleek-lab-htpasswd lab
# (you'll be prompted for the password; share that pair with the lab)

# TLS cert. Self-signed is fine for an internal tool:
sudo mkdir -p /etc/nginx/certs
sudo openssl req -x509 -nodes -newkey rsa:4096 -days 825 \
  -keyout /etc/nginx/certs/fleek.key \
  -out    /etc/nginx/certs/fleek.crt \
  -subj   "/CN=fleek.example.uab.edu"

sudo cp deployment/nginx.conf.example /etc/nginx/sites-available/fleek
# Edit server_name + ssl_certificate paths.
sudo ln -s /etc/nginx/sites-available/fleek /etc/nginx/sites-enabled/
sudo rm /etc/nginx/sites-enabled/default
sudo nginx -t && sudo systemctl reload nginx
```

Open ports 80 and 443 in the OpenStack security group for whatever
network range your lab uses.

### 5. First-time check

Browse to `https://<your-host>/`. You should see:

1. The browser's TLS-warning interstitial (self-signed cert) — accept once.
2. The nginx basic-auth dialog — enter the lab/<password> credential.
3. The supervisor's `/login` page — click "Create one".
4. Sign up with your username + password (+ optional Anthropic key).
5. The FLEEK UI loads. The "REMOTE" chip appears in the sidebar header,
   and the "account" chip next to it links to `/account` for key-management
   and sign-out.

## Day-to-day operations

### Add a user

Users self-signup via `/signup`. No admin work.

### Reset a user's password

```bash
sudo -u fleek sqlite3 /data/users/.fleek_supervisor.db
> .schema users
> -- For now there's no in-supervisor reset flow. Quickest path:
> --   1. delete the row, tell the user to re-sign-up with same name
> --      (their dir + caches are preserved)
> DELETE FROM users WHERE username = 'alice';
> .quit
```

A future supervisor `--admin-reset` flag would do this cleanly; for now,
the SQL approach works.

### Check who's using what

```bash
# Currently-running per-user fleek subprocesses:
sudo systemctl status fleek-supervisor
sudo journalctl -u fleek-supervisor --since "10 min ago"

# All registered users:
sudo -u fleek sqlite3 /data/users/.fleek_supervisor.db \
  "SELECT username, last_login, assigned_port FROM users ORDER BY last_login DESC NULLS LAST;"

# Disk usage per user:
sudo du -shc /data/users/*/
```

### Bump quotas / restart

Edit the systemd unit's `--user-quota-mb` flag, then:

```bash
sudo systemctl daemon-reload
sudo systemctl restart fleek-supervisor
```

A restart kills all running per-user fleek subprocesses (they re-spawn
on next request from each user). User data, sessions, settings persist
on disk and survive the restart.

### Pre-bake caches on a shared atlas

```bash
# As your admin user (NOT as fleek):
fleek /data/shared/datasets/pbmc_atlas.h5ad --no-browser
# Wait for "Cache: writing to disk" lines for UMAP/PCA/CSC/markers/etc.
# Ctrl-C once everything's cached. Verify:
ls /data/shared/datasets/pbmc_atlas.*
# pbmc_atlas.h5ad
# pbmc_atlas.fleek_cache.npz   ← UMAP + PCA + PaCMAP + Leiden
# pbmc_atlas.csc_cache.npz     ← sparse column index
# pbmc_atlas.annot_markers.json
# pbmc_atlas.fleek_cluster_genes.json
```

Now every user opening that atlas reuses these caches instantly. User-
recomputed variants land in their own `/data/users/<name>/caches/` so
the canonical files stay untouched.

### Backup

```bash
# /data/users/ is the entire user-state surface (db + secret + per-user
# dirs). /data/shared/ is whatever atlases you've placed there.
sudo tar -C /data -czf /backup/fleek-$(date +%F).tar.gz users/
```

## Limits / known issues

- **Memory.** Each running per-user fleek holds its own loaded dataset
  in RAM. With 32 GB and 5 concurrent users on 8 GB datasets, you'll
  pressure RAM. Mitigations: shorten `--idle-timeout-s` so unused
  processes exit faster, add cgroup `MemoryMax` per subprocess via a
  systemd template, or get more RAM.
- **No password reset UI.** Admin-only via SQL for now.
- **Upload streaming is body-buffered** in the supervisor. Fine for
  the 100 MB default quota; if you raise quotas above ~1 GB consider
  switching to a streaming proxy.
- **Per-user uids not enforced.** The supervisor and all per-user fleek
  subprocesses run as the same uid (`fleek`), and tenancy is enforced
  by per-user dir paths. Adequate for a trusted lab; not adequate if
  someone with shell access on the box could exploit a Python deserial
  bug to read another user's files. See "Hardening" below.
- **No rate-limiting on Claude API calls** per-user. If a user goes
  rogue with annotation requests they could run up that user's API
  bill (since each user has their own key, there's at least a
  per-user blast radius).

## Hardening (optional)

If your threat model needs more than "trusted lab":

- **Real OS-user isolation.** Run each user's fleek as a separate uid.
  Requires changes to the supervisor (sudo / setuid wrapper) and a
  `useradd` per signup. Significant lift; only worth it if your users
  could realistically attack each other.
- **UAB SSO.** Replace nginx basic-auth with Shibboleth at the edge.
  Open a ticket with UAB RC.
- **Disk quotas via filesystem.** Linux `quota` package gives kernel-
  enforced caps that the supervisor's app-level check can't bypass.
  Useful if a user finds a path that writes outside the supervisor's
  quota check.
- **Memory limits.** A systemd slice or per-subprocess cgroup with
  `MemoryMax=8G` prevents one user from OOM-killing the rest.
