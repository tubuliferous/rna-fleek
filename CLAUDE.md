# RNA-FLEEK: Fast Lightweight Expression Exploration Kit

## Project Overview

RNA-FLEEK is a browser-based single-cell RNA-seq visualization and analysis tool. It consists of:

- **`fleek.html`** — Single-file frontend (~10000 lines). Three.js WebGL renderer for 2D/3D point clouds of cells. All JS is inline inside `<script>` tags. No build system, no framework.
- **`fleek_server.py`** — Python HTTP server (~4000 lines). Loads h5ad files, computes embeddings (UMAP/PCA/PaCMAP), DEG (single-cell + pseudo-bulk), clustering, serves binary data to the client.
- **Utility scripts**: `download_census.py`, `preprocess_census.py`, `download_markers.py`
- **GitHub repo**: `tubuliferous/rna-fleek`

## Architecture

### Frontend (`fleek.html`)
- Single HTML file with inline CSS (`<style>`) and JS (`<script>`)
- Three.js r128 from CDN for WebGL rendering
- Binary protocol (SCRN format) for initial data transfer
- Virtual scrolling for gene variability list and cell type list
- Session state persistence via `sessionStorage`
- All UI in a left sidebar; canvas fills remaining space
- Split Views panel for managing multiple named view splits
- Split-screen mode: two views side-by-side with independent/tethered cameras
- Custom dropdown menus (div-based) replacing native `<select>` popups for theme consistency
- 2D point picking bypasses raycaster — pure pixel-space projection for accuracy

### Server (`fleek_server.py`)
- stdlib `http.server` with `ThreadingMixIn` — no Flask/FastAPI
- Loads h5ad via anndata/scanpy
- Binary SCRN format for init payload (header JSON + float32/int32 arrays)
- Background threads for CSC index building, gene variability scoring
- API endpoints return JSON (annotations, search, DEG) or binary (init, gene expression, embeddings)

### API Endpoints
| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/init` | GET | Initial binary payload (cells, clusters, embeddings) |
| `/api/gene` | POST | Gene expression values (float32 binary) |
| `/api/search` | POST | Gene name search |
| `/api/gene-var` | POST | Gene variability ranked list |
| `/api/cluster-genes` | GET | Per-cluster gene rankings (all genes, full metrics + Cohen's d) |
| `/api/load` | POST | Load a new h5ad file |
| `/api/annotate` | POST | Marker DB cell type annotation |
| `/api/annotate-llm` | POST | Claude AI cell type annotation |
| `/api/lineage` | POST | Developmental lineage tree via Claude |
| `/api/embedding` | POST | Fetch specific embedding |
| `/api/compute-embedding` | POST | Trigger UMAP/PaCMAP computation |
| `/api/export` | POST | Export h5ad subset |
| `/api/browse` | POST | File browser for server filesystem |
| `/api/clear-cache` | POST | Delete specific cache file |
| `/api/cache-settings` | POST | Get/set cache location mode |
| `/api/progress` | GET | Polling for load/compute progress |
| `/api/abort` | POST | Cancel in-progress operation |
| `/api/unload` | POST | Free dataset from memory |
| `/api/pseudobulk-deg` | POST | Pseudo-bulk DEG (selection groups + replicate column) |
| `/api/heartbeat` | POST | Client heartbeat for auto-unload |

## Critical Conventions

### SVG in `<script>` blocks
**ALL closing HTML tags inside `<script>` must be escaped.** The HTML parser sees `</` followed by a letter as a potential end tag, which prematurely terminates the script block — silently breaking all JS that follows (Three.js controls, rendering, highlights, everything).

```javascript
// WRONG — breaks the script parser:
var icon = '<svg ...></svg>';

// CORRECT:
var icon = '<svg ...><\/svg>';
```

This applies to `</svg>`, `</button>`, `</div>`, `</span>`, etc. After any edit that adds HTML strings inside JS, verify with:
```bash
python3 -c "
import re
with open('fleek.html') as f: c=f.read()
s=c.find('<script>')+8; e=c.find('</script>',s)
print(len(re.findall(r'(?<!\\\\)</(svg|button|div|span)', c[s:e])), 'unescaped tags')
"
```

### Virtual Scrolling Pattern
Both the gene variability list and cell type list use virtual scrolling:
- A **spacer div** sets total height (`height: rowCount * rowHeight`)
- A **viewport div** (position:absolute inside spacer) contains only visible rows
- **onscroll handler** calculates visible range and re-renders only those rows
- **Range caching** (`_lastStart`/`_lastEnd`) skips re-render if range unchanged
- **clientHeight=0 fallback**: When container is just un-hidden, `clientHeight` may be 0. Use a fallback value (280px for cell types, 200px for genes)

Gene list: fixed 18px row height, simple `floor(scrollTop/rowH)` math.
Cell type list: variable heights (24px base + 16px per annotation line), pre-computed cumulative Y positions, binary search for visible range.

### Cache System
Two modes controlled by `CACHE_MODE` global:
- **"dataset"**: Caches stored alongside .h5ad file (default)
- **"server"**: Caches stored in `~/.fleek_cache/` (or `--cache-dir`)

Always use the helper functions:
- `_cache_read_path(suffix)` — returns `(path, writable)` or `(None, False)`. Checks dataset dir first, then server dir.
- `_cache_write_path(suffix)` — returns best writable path. In dataset mode, tries dataset dir, falls back to server dir if not writable.
- **Never** use `Path(LOADED_PATH).with_suffix(suffix)` directly.

Cache suffixes:
- `.fleek_cache.npz` — UMAP, PCA, PaCMAP, Leiden (numpy arrays)
- `.csc_cache.npz` — Sparse column index for fast gene lookup
- `.annot_markers.json` — Marker DB annotations
- `.annot_llm.json` — Claude annotations
- `.annot_lineage.json` — Developmental lineage tree

### Memory Safety (Large Datasets)
The DEG computation path must **subsample BEFORE copying**. On a 1.6M × 61K dataset, copying the full adata for normalization doubles memory and triggers OOM kill. The pattern is:
1. Compute valid/small clusters on indices only (no copy)
2. Build stratified subsample on indices only (~50k cells)
3. `adata[use_idx].copy()` — copy only the small subset
4. Normalize in-place on the small copy
5. Run DEG on the small copy

### Lineage Tree
- Empty clusters (0 cells) are **excluded** from the lineage prompt — sending hundreds of empty clusters wastes tokens and can cause truncation
- `max_tokens` is 16384, timeout is 300s
- Check `stop_reason == "max_tokens"` for truncation detection
- Lineage validation only checks non-empty cluster_ids

### Cell Type List Interactions
Both flat list and lineage tree should have consistent modifier key behavior:
- **Click**: Toggle visibility
- **Alt+click**: Solo/unsolo (saves & restores previous visibility state)
- **Shift+click**: Select cells into active group
- **Alt+Shift+click**: Select only this node's direct clusters (lineage only)

In lineage view, **dot** and **name** have different scopes:
- **Dot click/Alt+click**: Acts on this node's direct clusters only
- **Name click/Alt+click**: Acts on this node + all descendants
- **Caret (▸/▾) click**: Collapse/expand branch (Alt+click: recursive)
- Children sort by current `typeSortMode` (alpha or count)

### Theme Rules (Dark / Light)
- **All CSS must work in both themes.** Use `var(--text)`, `var(--bg3)`, etc. — never hardcode colors that only look right in one theme.
- **`color-scheme`**: `:root` has `color-scheme:dark`, `.light` has `color-scheme:light`. Native form elements (`select`, `input`, scrollbars) inherit from this. Never hardcode `color-scheme:dark` on individual elements — use `inherit`.
- **`<meta name="color-scheme">`**: toggled by `applyTheme()` — tells the browser which native UI style to use (affects select dropdown popups on Safari/WebKit).
- **When adding any visible element**: verify it looks correct in both themes. Use browser dev tools to toggle `.light` on `<body>`.
- **Halo/glow colors**: dark mode uses additive blending (`ONE, ONE`), light mode uses alpha blending (`SRC_ALPHA, ONE_MINUS_SRC_ALPHA`) with a dark navy color. These are set in `applyTheme()`.
- **Expression pulse target**: dark mode pushes towards white, light mode towards deep navy `(0.02, 0.05, 0.18)`.
- **Slice pill colors**: inline backgrounds set by `_syncSlicePills()` use theme-aware `onColor`/`offColor` with `isLight` check.

### UI Conventions
- Empty clusters: hidden behind collapsed "▸ N empty clusters" toggle
- Selection group badges: eye icon + visible count, slashed-eye + hidden count
- Cache badges: grouped under "techniques" (dashed border, non-deletable) and "caches" (solid border, deletable with × on hover)
- Hint bar: centered bottom of canvas (`left:50%; transform:translateX(-50%)`)
- Load time: persisted in `sessionStorage` as `fleek_infoText`, survives page refresh
- Double-click on cell: center + zoom to 35% of cloud radius (fixed distance every time)

### Gene Variability
Two-stage locally standardized trend residual:
1. **Stage 1**: Fit polynomial to log(mean) vs log(Fano) trend
2. **Stage 2**: Fit scatter envelope |residual| vs expression to get local sigma
3. **Score**: z = raw_residual / max(local_sigma, 0.1)

Red-flagged genes: bitmask (bit 0: no expression, bit 1: low mean, bit 2: high mean, bit 3: low fraction)

### Design system (tokens + primitives)
Every new style rule MUST resolve size/spacing/radius/duration values to the design tokens in `:root` at the top of `fleek.html`. If a hard-coded pixel value is tempting, the right answer is almost always "use a token"; if no existing token fits, add a new one rather than inlining.

**Token families:**
- Color: `--text / --text2..4`, `--bg / --bg2..3`, `--bg-hover`, `--bg-active`, `--border / --border2`, `--accent / --accent2`, `--danger`, `--success`. Light-mode overrides under `.light`.
- Font size: `--fs-xs` (9px) · `--fs-sm` (10px) · `--fs-md` (11px) · `--fs-lg` (12px) · `--fs-xl` (13px).
- Font weight: `--fw-regular / medium / semibold / bold`.
- Line height: `--lh-tight / base / relaxed`.
- Spacing (4px-ish base): `--sp-0..8`.
- Icon size: `--ic-xs` (10px) · `--ic-sm` (11px) · `--ic-md` (12px, default) · `--ic-lg` (14px) · `--ic-xl` (16px).
- Radius: `--rad-sm` (3px) · `--rad-md` (4px) · `--rad-lg` (5px) · `--rad-xl` (8px) · `--rad-pill`.
- Duration: `--dur-fast` (0.1s) · `--dur-med` (0.2s) · `--dur-slow` (0.3s).
- Icon interaction opacity: `--op-rest` (0.7) · `--op-active` (1) · `--op-disabled` (0.3).

**Component primitives:**
- `.btn-icon` — icon-only button. Baseline: `--op-rest` opacity, brightens to `--op-active` on hover, color bumps to `--text2`. Modifiers:
  - `.btn-icon.download` → hover = `--success`
  - `.btn-icon.danger`   → hover = `--danger`
  - `.btn-icon.accent`   → hover = `--accent2`
  - Size: `.sm` / `.md` (default) / `.lg` / `.xl`.
  - `.btn-icon--row-hover` → starts invisible, fades in with the surrounding row's hover.
- `.pill` — rounded tint chip. Variants: `.accent` / `.danger` / `.success` / `.warning`. Used for load badges, gene indicator, status chips.
- Legacy alias `.icon-btn` is kept working but NEW markup should prefer `.btn-icon`.

**Scrollbar strategy:**
- Global native scrollbars: 7px, `--border` thumb, thin (see `::-webkit-scrollbar` + `scrollbar-width:thin`).
- Two exceptions (`#sidebar-upper`, `#sidebar-right-upper`) draw their OWN custom overlay thumb via `setupCustomScrollbars()` → native is explicitly hidden on these containers with `scrollbar-width:none` + `::-webkit-scrollbar{display:none}`, otherwise macOS shows both during scroll and they visually stack.
- Rule: if a container gets a custom thumb, hide its native scrollbar. Never let both render simultaneously.

**Semantic colors for actions:**
- Destructive (delete / close / unload): `--danger`
- Download / file export: `--success`
- Neutral interactive: `--text2` (hover only)
- Active/selected state: `--accent2`
Any other ad-hoc color (`#3b82f6`, `#f59e0b`, hardcoded `#f87171`) should be replaced with the matching var.

## Testing

```bash
# Start server
python fleek_server.py path/to/data.h5ad --port 8080 --api-key sk-ant-...

# With options
python fleek_server.py data.h5ad --host 0.0.0.0 --port 8080 --cache-mode server --no-browser

# Quick validation
python3 -c "import ast; ast.parse(open('fleek_server.py').read()); print('Server syntax OK')"
```

After editing `fleek.html`, always verify:
1. Script tags balanced: `grep -c '<script' fleek.html` should equal `grep -c '</script>' fleek.html`
2. No unescaped closing tags in script (see SVG section above)
3. Server syntax: `python3 -c "import ast; ast.parse(open('fleek_server.py').read())"`

### Pseudo-bulk DEG
Uses selection groups as conditions (not obs columns). User selects a single "Replicate column" dropdown.
- Server receives cell indices per group + replicate column name
- Aggregates raw counts per (group × replicate) combination
- Runs pyDESeq2 (preferred) or Welch's t-test on CPM (fallback)
- Requires ≥2 replicates per group, ≥10 total counts per gene
- Results displayed in same volcano plot + table as single-cell DEG

### Split Views
- `_views[]` array holds named snapshots (slice filter, visibility, deleted cells)
- `_splitViewIds = [leftId, rightId]` tracks which views render in each half
- `_liveViews[]` holds GPU resources only (pts, ringPts, haloPts, cameras)
- `recolor()` in split mode: saves active → loops each view (load snapshot → `_recolorView`) → restores active
- Selections are fully global — shared across all split views
- Split Views panel in sidebar (between Selection Groups and Slice)

### Auto-quit (heartbeat)
- Renamed from "auto-unload" but variables/settings keep the `auto_unload`
  name for session-storage compatibility. Default: ON, 20-second timeout.
- Client sends heartbeat POST every `timeout/2` seconds when enabled.
- Server cancels/resets a `threading.Timer` on each heartbeat.
- If no heartbeat within `timeout + 0.5s`, server: (a) unloads any loaded
  dataset via `_reset_all()` + `gc.collect()`, then (b) `os._exit(0)` the
  process. Runs even when no dataset is loaded — the point is to kill the
  server process so a naive user who closes the browser tab doesn't leave
  a background server holding a port + RAM.
- A page reload within the grace window reconnects and cancels the quit.
- Disable via Settings → "Auto-quit server on disconnect" (for shared
  servers, dev sessions, or multi-tab users who want the server to persist).

### Organism Detection
- `_detect_organism()` returns `(name, reason)` tuple
- Primary: checks for 18 species-specific marker genes (human UPPERCASE, mouse Title case)
- Fallback: gene name casing heuristic on first 200 alpha-only names
- Used for GO database loading and displayed in Data panel info line

## Pending Features / Known Issues
- Gene set analysis (GSEA/pathway scoring)
- PLY export for Blender
- Tissue hint input for Claude annotation
- README needs updating with all new features

### Slice-specific embeddings (Phase 3)
Slice panel (Phase 2) is implemented: obs metadata columns shown as checkboxes, cells outside the filter are exiled client-side, DEG naturally respects the slice. Phase 3 is not yet built:

- **Slice-specific UMAP/PaCMAP**: When a filter is active, compute a new embedding on just the matching cells. Cache keyed by canonical slice string (e.g. `organ=D|timepoint=D8`). Cache files alongside the h5ad as `.fleek_cache_slice_<hash>.npz`.
- **Comparability constraint**: Naive per-slice UMAPs are in different coordinate spaces and cannot be directly compared. The global embedding with visual masking (current behaviour) IS comparable. Slice-specific embeddings are for revealing internal structure of a single condition, not for cross-slice positional comparison.
- **Server endpoint**: `/api/slice-embed` — accepts slice definition + method, subsets adata, computes embedding, caches, streams progress.
- **Client**: "Re-embed this slice" button in Slice panel (only when filter is active); progress indicator; load-badge showing which embedding is active (global vs slice-specific).
- **Leiden re-clustering**: Optionally re-cluster the slice subset after re-embedding to find condition-specific clusters.
