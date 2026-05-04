# RNA-FLEEK: Fast Lightweight Expression Exploration Kit

## Project Overview

RNA-FLEEK is a browser-based single-cell RNA-seq visualization and analysis tool. It consists of:

- **`rna_fleek/fleek.html`** — Single-file frontend (~22000 lines). Three.js WebGL renderer for 2D/3D point clouds of cells. All JS is inline inside `<script>` tags. No build system, no framework. (Top-level `fleek.html` is a symlink to this file.)
- **`rna_fleek/server.py`** — Python HTTP server (~7000 lines). Loads h5ad files, computes embeddings (UMAP/PCA/PaCMAP), DEG (single-cell + pseudo-bulk), clustering, serves binary data to the client. (Top-level `fleek_server.py` is a thin shim.)
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

### Server (`rna_fleek/server.py`)
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
| `/api/pathway-databases` | GET | List available pathway databases (Hallmark, Reactome, GO BP/MF/CC, user GMTs) |
| `/api/pathway-ora` | POST | ORA hypergeometric test (input genes + database id → BH-corrected results) |
| `/api/heartbeat` | POST | Client heartbeat for auto-unload |

## Critical Conventions

### Single source of truth: `rna_fleek/`
**All server and frontend code lives in `rna_fleek/`.** There are no mirrors.

| Concern | Canonical file |
|---------|---------------|
| Server | `rna_fleek/server.py` |
| Frontend | `rna_fleek/fleek.html` |

The top-level `fleek_server.py` is a **3-line shim** that does `from rna_fleek.server import main; main()` — it exists only so the historical `python fleek_server.py …` dev workflow keeps working without duplication. The top-level `fleek.html` is a **symlink** to `rna_fleek/fleek.html`. Do not edit either of the top-level files directly; edit `rna_fleek/`.

The PyInstaller specs and `pip install rna-fleek` already point at the package layout; nothing else needed to change when the mirrors were collapsed.

**Version bumps** still touch four spots (no fifth `fleek_server.py FLEEK_VERSION` constant — the shim has none): `pyproject.toml`, `rna_fleek/__init__.py`, `rna_fleek/server.py` (`FLEEK_VERSION = "x.y.z"`), and `rna_fleek/fleek.html` (`var FLEEK_VERSION="vx.y.z"`). All four must match.

### Every version bump gets a Changelog entry
**Whenever you bump the version (the four files above), you MUST also add a new entry to the Changelog section in `rna_fleek/fleek.html`.** The Changelog lives at `<div class="help-sec" id="help-changelog">` near the bottom of the help body and is surfaced via the "Changelog" button in the help header (left of "Send feedback") and the search-palette `open-changelog` entry. It's the in-app release-notes feed users open to find out what's new.

**Format:** add a fresh `<div class="changelog-entry">` block at the TOP of the changelog body (most-recent first). Use this template:

```html
<div class="changelog-entry">
  <h4>vX.Y.Z</h4>
  <ul>
    <li><strong>Headline change:</strong> One or two sentences describing the user-visible behavior. Reference UI elements by name when possible.</li>
    <li>Smaller bullets for secondary changes — keyboard shortcuts, fixes, polish.</li>
  </ul>
</div>
```

**What to include:**
- New features, panels, dialogs, or click affordances.
- Keyboard-shortcut changes (especially remaps — these break muscle memory).
- Visible behavior changes (new defaults, removed/relocated controls, semantic shifts).
- Bug fixes that users would have noticed (silent feature gating, broken-looking states).
- Performance changes that change how the app *feels* (load time wins, big interaction snappiness improvements).

**What NOT to include:**
- Pure refactors, file moves, internal renames.
- CI / build / packaging plumbing.
- Per-commit minutiae — group by theme. The commit log on GitHub has the full granular history; the in-app changelog is for "what does this mean for me as a user".

**Tone:** matter-of-fact, second-person where natural ("Click <em>Hidden N</em> to open …"), worked examples for compound interactions. Match the voice of existing entries — read the most recent ones before writing a new one. Use `<code>`, `<strong>`, `<em>`, `<span class="help-key">` for inline emphasis (same idiom as the rest of the help panel).

**The version number's job is to point at the changelog, not stand alone.** Bumping the version without writing the entry is the same kind of half-job as touching `rna_fleek/server.py` without bumping the version — it leaves the user wondering what changed.

### Stay inside the requested scope — strictly
**When the user asks for a specific change, change ONLY what was asked. Nothing else.** No silent improvements to adjacent code, no restyling nearby UI, no renaming unrelated variables, no refactoring surrounding logic, no "fixing" things that look wrong but weren't part of the request. The user is the only judge of whether a related change is welcome — and they have to find each unwelcome change manually, then ask for it to be reverted. That cycle is far more expensive than any alleged improvement is worth, and it has happened repeatedly in this project. The drift is the single largest source of regressions.

**This rule overrides any conflicting impulse.** Even when the change "obviously" needs follow-up work, even when the surrounding code is "clearly" inconsistent, even when polishing it would only take a moment — STOP. The right action is to call it out in your reply and ask the user whether to address it, not to fold it silently into the same edit.

**The only exceptions are absolutely essential to the requested change:**
- A new control needs a new state variable, event handler, and CSS rule — those are part of the change.
- A bug fix in function A requires also patching function B that calls A and would now misbehave — that's part of the change.
- A removed feature requires deleting its now-orphan call sites — that's part of the change.

If you can describe the side change as "while I was in there, I also …" — it does NOT belong in the edit. Pull it out.

**Concrete examples of forbidden drift, all of which have actually happened:**
- Renaming a function or variable "for clarity" while fixing a bug in it.
- Tightening a CSS rule to "match the design system" while wiring up a new feature nearby.
- Reformatting an HTML block (line breaks, attribute order) while editing one attribute inside it.
- Replacing inline styles with classes while moving a node.
- Adding `event.stopPropagation()` to a handler that doesn't strictly need it because "it's safer."
- Removing seemingly-unused markup or click handlers without verifying they really aren't reached.
- Swapping a `<span>` for a `<select>`, or vice versa, because the new one "fits better" — when the user just asked for a label change.

**Visual changes are especially dangerous.** Users notice button sizes, padding, spacing, icon swaps, color tweaks, and font changes immediately. Never change a control's appearance unless that's the explicit ask. Same for keyboard shortcuts and event-handler behavior — users build muscle memory around these and any silent change feels like the app broke.

**When wrapping or restructuring HTML (a frequent regression source):** any wrapper around a flex/grid item inherits the parent's layout role. If the original element had `flex:1` or any width-defining rule, the wrapper now needs that rule and the inner element needs to fill the wrapper. Test the visual result before submitting — the example here was wrapping the disabled DEG button in a `<span>` for tooltip-on-disabled support; the span didn't get `flex:1`, so the toolbar collapsed the button to text-width.

**Self-check before submitting any edit:** for each line you changed, ask "did the user explicitly request this change, or is it logically required by the user's request?" If the honest answer is "no" — revert that line. If you're unsure, revert it AND mention it in your reply for the user to weigh in.

**Why:** the user's mental model is "I asked for one thing; one thing changed". Every side effect breaks that model and makes every future change harder to trust. Trust, once spent, is hard to rebuild.

### Keep Help/About in sync with user-facing changes
**When adding or changing a user-facing feature, update the corresponding Help panel section in `fleek.html` (inside `#help-body`) as part of the same change.** This is a workflow rule, not an optional polish step — the Help panel is the app's only built-in documentation, and each drift is a small paper cut for users.

**Applies to:**
- New panels, controls, toggles, dropdowns, chips, keyboard shortcuts.
- New plot types, export formats, new modal dialogs.
- Semantic changes to existing controls (e.g., threshold flips, default changes, renames).
- New visual encodings (e.g., size / width / opacity mapping to data).
- New `obs` / metadata interpretations or filter rules.

**Doesn't apply to:**
- Pure refactors, performance tweaks, internal plumbing changes.
- Bug fixes that restore previously-documented behavior.
- Theme/styling tweaks that don't change meaning.

**How to update:**
- If an existing `<div class="help-sec" id="help-…">` covers the feature, edit that section's `<div class="help-sec-body">` inline.
- If the feature is novel, add a new `<div class="help-sec" id="help-…">` block **and** add a matching `<a href="#" onclick="event.preventDefault();helpJump('help-…')">…</a>` entry to the TOC in `#help-body` so it's navigable.
- Help sections use `<p>`, `<ul>`, `<table class="help-table">`, `<code>`, `<strong>`, `<em>` — keep the tone consistent with existing sections: matter-of-fact, explain the *why* where it's non-obvious, and include worked examples for compound interactions.

### Keyboard-shortcut symbols must be OS-aware everywhere
Whenever a shortcut hint is shown to the user — Help panel, search palette `shortcut` field, button `title` tooltips, mid-action `_flashAction` toasts, tooltips on disabled controls — the displayed text MUST resolve to the OS-correct modifier name. Mac users see `⌘` / `⌥` / `⌃` / `⇧`; Linux/Windows users see `Ctrl` / `Alt` / `Shift`. **Mixing them — e.g. showing `⌥` on a Linux box because the source string was hard-coded with the Mac glyph — is a bug, not a stylistic choice.**

**The mechanisms:**
- **Static help text inside `#help-body`** is rewritten at boot by `_patchOS()` (the function does the `⌘ → Ctrl`, `⌥ → Alt`, etc. swaps when `_isMac` is false). Authors can write the source HTML with Mac glyphs and trust the runtime swap; just keep the source readable.
- **Static button / element `title` tooltips set in HTML** are also caught by `_patchOS()` for elements explicitly listed there. When adding a new button whose tooltip mentions a modifier, either use the `_MOD` / `_ALT` globals from JS, or add the element to the `_patchOS()` patch list, or write the tooltip with `Ctrl`/`Alt` (which both Mac and Linux users can read).
- **Search-palette `shortcut` field** in entries (`_searchRegistry`) is rendered through `_osShortcutLabel(s)` (helper near `_isMac`), which converts Mac glyphs to OS-correct text at render time. **Author entries with Mac glyphs as the source-of-truth** (`shortcut: "⌥ \\"`) and let the renderer translate.
- **Mid-action toasts via `_flashAction(label, key)`** — the `key` arg is the small chip shown next to the action name. Use the OS-aware globals (`_MOD`, `_ALT`) for any modifier portion, or write `Ctrl`/`Alt`. Don't hard-code `⌘`/`⌥` here unless `_flashAction` itself does the conversion (it doesn't today).
- **Dynamically-built `title` strings in JS** (e.g. the DEG button's "needs 2 groups…" tooltip) — same rule. Use `_MOD` / `_ALT` or use OS-neutral words.

**Self-check after adding a shortcut:**
1. Does the user see the shortcut anywhere? (Help text, palette, button title, flash chip, error message, etc.)
2. Does each of those surfaces use either an OS-aware mechanism (`_patchOS`, `_osShortcutLabel`) or write plain `Ctrl`/`Alt` text?
3. If you wrote `⌘`/`⌥` literally, does the surrounding code path run through one of the converters?

If any answer is "no", you'll have Linux/Windows users seeing Mac glyphs they don't recognize. The fix is one of: route through a converter, swap to `Ctrl`/`Alt` literals, or use `_MOD`/`_ALT`.

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
- **Both themes are first-class. Every new feature MUST work and look correct in light mode AND dark mode** — not as an afterthought. When adding a UI element, plot, modal, indicator, etc., explicitly read `document.body.classList.contains("light")` for any decisions that need different colors / opacities / contrasts in the two themes (e.g., density-curve fill alpha that reads as faint on white may need to be bumped up vs the dark-mode value; SVG renderers should pick axis/grid/label/background colors per theme rather than assume one).
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

### Consistency with the existing design
**Any new UI element should look like it was always part of the app.** Before introducing a new control, scan the existing HTML/CSS for an equivalent and reuse it — don't invent parallel styling.

- **Toggles / booleans: prefer switches over checkboxes.** When a feature can be implemented with either, use the `.se-wrap` / `.se-switch` / `.se-track` switch primitive already used for Slice "Show Empty" and Settings toggles. Markup template:
  ```html
  <label class="se-wrap" title="What the toggle does">
    <span class="se-label">Label text</span>
    <span class="se-switch"><input type="checkbox" id="..." onchange="..."><span class="se-track"></span></span>
  </label>
  ```
  Reserve native `<input type="checkbox">` for rare cases where a switch is genuinely wrong (e.g. multi-select lists with many options).
- **Icons:** reuse existing SVG icons from the file rather than importing new ones when a semantically-equivalent one already exists. Match stroke width, corner style, and sizing (`--ic-xs..xl`).
- **Fonts / sizing:** use the `--fs-*` font-size tokens and existing label classes (`.lbl`, `.help-sec-title`, etc.). Don't hardcode pixel font sizes.
- **Button hover colors:** destructive → `--danger`; accept/download → `--success`; neutral → `--text2`; active/selected → `--accent2`. Match the existing `.btn-icon` / `.pill` primitives.
- **Modal dialogs:** model after `#help-panel` / `#confirm-panel` / `#fb-modal` (rounded corners, `var(--bg2)` surface, header with × close, footer with right-aligned action buttons, Esc-to-close, click-outside-to-close).
- **Dropdowns:** model after `.gv-dd` (custom dropdown with popover list) — don't introduce native `<select>` popups unless the existing custom ones are genuinely unsuitable.
- **Dropdown / chevron alignment:** the dropdown chevron glyph (`▾`, `▼`, or `&#9660;`) MUST be vertically aligned with the label text it sits next to. Unicode arrows have their visual mass low in their em-box, so a same-font-size chevron next to a label reads LOWER than the label's optical center — verify visually and apply a small `transform:translateY(-Npx)` to the chevron when needed (the plot-popup title chevron uses `-1px`; smaller-font chevrons may need a smaller offset). This applies to every dropdown chevron going forward — not just custom ones. When adding a new dropdown, eyeball-align the chevron with the label before submitting; "close enough" is not — users notice the misalignment immediately.
- **Cluster / cell-type palette: avoid extreme luminances.** When generating or assigning cluster colors, keep entries in a moderate luminance band (roughly 0.18 ≤ L ≤ 0.82, where `L = 0.299r + 0.587g + 0.114b`). Very dark colors disappear against the dark-mode background; very light colors disappear against the light-mode background. The client-side `_clampLuminance` helper inside `cc()` enforces this at read time so every callsite (canvas, legend, tooltips, exports) gets clamped colors automatically. If a generator is producing extremes, clamp toward gray rather than letting them through.
- **Gradients: same width + height + end-roundness + border.** All gradient bars in the program (the GENE-panel single-gene gradient, the Colocalized gradient, the Background-shade slider) share `width: 100%`, `height: 13px`, `border-radius: 5px`, and `border: 1px solid var(--border)` (with `box-sizing: border-box` so the border doesn't subtly shift the inner gradient's width). Future gradient strips MUST follow the same convention.
- **Capitalization: sentence case for all UI text.** Every user-facing label, button text, panel header, dropdown option, tab label, search-palette entry, flash-action toast, and help-section title uses sentence case (first word capitalized, rest lowercase). Acronyms stay uppercase (DEG, GO, AI, UMAP, PCA, FLEEK, Claude). Single-word labels are unaffected. Examples: "Selection groups" (not "Selection Groups"), "Cell types" (not "Cell Types"), "Gene expression" (not "Gene Expression"), "Render plot" (not "Render Plot"). When adding a new label, default to sentence case unless you have a strong reason — Title Case mid-program reads as inconsistent immediately.

If a situation genuinely needs a new primitive (because nothing existing fits), add tokens and classes in the `:root` / top CSS block rather than inlining styles, so the new primitive is reusable by the next feature.

### Gene Variability
Two-stage locally standardized trend residual:
1. **Stage 1**: Fit polynomial to log(mean) vs log(Fano) trend
2. **Stage 2**: Fit scatter envelope |residual| vs expression to get local sigma
3. **Score**: z = raw_residual / max(local_sigma, 0.1)

Red-flagged genes: bitmask (bit 0: no expression, bit 1: low mean, bit 2: high mean, bit 3: low fraction)

### Plots: every plot must ship with a download button
**Every plot in FLEEK — current and future — must have a download icon next to it that exports the on-screen rendering to a portable image file.** Use the same SVG download glyph used everywhere else (the `<polyline points="7 10 12 15 17 10"/>...` "tray + down arrow" path) on a `.icon-btn.download` (or `.btn-icon.download`) button.

**Format choices:**
- SVG-rendered plots (raincloud, violin, dot plot, anything assembled into the `#plot-canvas` popup): export as **standalone SVG** via `_plotSave()`. The function combines facet panels into a single outer `<svg>` and bakes a font-family on the root so the file renders correctly outside the browser.
- Canvas-rendered plots (volcano, anything `<canvas>`): export as **PNG** via `canvas.toDataURL("image/png")` then trigger a hidden `<a download>` click — see `downloadVolcanoPNG()` for the pattern.

**Filename convention:** every download in FLEEK uses one template — `<dataset-stem>_<what>_<detail>_<YYYY-MM-DD_HH-MM-SS>.<ext>` — built via the shared helpers `_dlDatasetStem()`, `_dlSafe()`, `_dlTimestamp()` (defined near `_plotSave` in fleek.html). The dataset stem comes first so a downloads folder groups by dataset; the timestamp guards against silent overwrites on repeated exports. New code MUST use the helpers, not roll its own filename string. Examples: `colon_atlas_raincloud_FOXP3_2026-04-27_14-30-22.svg`, `colon_atlas_dotplot_FOXP3-CD4-IL2_2026-04-27_14-30-22.svg`, `colon_atlas_volcano_singlecell_TumorA_vs_TumorB_2026-04-27_14-30-22.png`, `colon_atlas_DEG_pseudobulk_A_vs_B_2026-04-27_14-30-22.tsv`. Server-side h5ad subsets follow the same template (`<parent-stem>_subset_<group>_<N>cells_<ts>.h5ad` from `export_h5ad_subset`). When you add a new download path, use the helpers; don't introduce a new convention.

**When you add a new plot:**
1. Place a `<button class="icon-btn download">` (or `btn-icon download`) in the plot's header / toolbar with the standard download SVG glyph and a `title="Download <plot type> (<format>)"`.
2. Wire its `onclick` to a `download<PlotType><Format>()` helper (e.g. `downloadVolcanoPNG`, `_plotSave`).
3. Add a one-line bullet in the Help panel's plot section saying the download exists.
4. If the plot has multiple panels (facets / splits), serialize **all panels** into one image, not just the first — see `_plotSave`'s flex-row → translated `<g>` pattern.

### Search palette: every named feature must register
**Every named, user-visible feature in FLEEK — panels, modes, toggles, actions, plots, picker rows — MUST have a corresponding entry in `_searchRegistry` (single array near the top of fleek.html's script section, search for `var _searchRegistry`).** The palette (triggered by `/` or `Cmd/Ctrl+K`) is the keyboard-driven jump-anywhere index; if a feature isn't in the registry, the user can't find it that way and the registry's promise as the single source of truth breaks.

**Applies to:**
- New panels (jump to + expand if collapsed).
- New modes (color modes, plot types, DEG modes, pathway sources).
- New toggles, switches, threshold sliders, gradient pickers, scheme pickers.
- New actions (anything bound to a button or shortcut: render plot, run DEG, save session, send feedback, reset camera, etc.).
- New keyboard shortcuts (the entry's `shortcut` field surfaces them in the palette).

**Doesn't apply to:**
- Pure refactors / internal plumbing.
- Per-row dynamic UI (cluster rows, gene rows, gene chips, pathway rows in a results list — these are content, not features).
- Dialogs that are reachable only as a side-effect of another action (Confirm, Prompt — not separately searchable).

**Entry shape:**
```js
{
  id:        "unique-stable-id",      // used by the recent-actions list, NEVER reuse / rename
  label:     "Render Raincloud",      // primary text shown in the palette
  secondary: "Plots",                 // panel/group context, faint
  keywords:  "distribution density",  // extra matchable terms (no need to repeat label words)
  shortcut:  "G",                     // optional — surfaces the existing key binding in the palette
  available: function(){ return D!==null; },  // optional — return false to gray out (still findable)
  panel:     "plot-sec",              // optional — section id to expand before navigating
  target:    "plot-render-btn",       // optional — element id to spotlight after the action
  action:    function(){ ... }        // what Enter does (toggle, mode change, focus, click, etc.)
}
```

**How to apply:**
1. When adding a feature, write the entry alongside the markup/handler. Use a short, lowercase `id` (e.g. `mode-coloc`, `plot-render`, `settings-api-key`) — the id is referenced by the recent-actions storage and renaming it loses recents.
2. Match-target hierarchy: `target` (specific element) > `panel` (whole panel). If both are present the panel is expanded and the element is spotlit.
3. For toggles, `action` should perform the toggle (e.g. `document.getElementById('show-non-expr').click()`) so the palette's behaviour matches the user clicking the control.
4. For `available()`, only gate things that genuinely don't work in the current state (e.g. coloc R slot unavailable in single-gene mode). Don't gate by "panel collapsed" — the palette can expand panels.
5. After adding the entry, smoke-test: open the palette, type a few characters of the label, hit Enter, and verify the spotlight lands on the right thing.

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
python3 -c "import ast; ast.parse(open('rna_fleek/server.py').read()); print('Server syntax OK')"
```

After editing `fleek.html`, always verify:
1. Script tags balanced: `grep -c '<script' fleek.html` should equal `grep -c '</script>' fleek.html`
2. No unescaped closing tags in script (see SVG section above)
3. Server syntax: `python3 -c "import ast; ast.parse(open('rna_fleek/server.py').read())"`

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

#### Color modes (current set)
- `mode="cluster"` — per-cluster palette (default). No per-cell scalar.
- `mode="gene"` with `_geneMode="single"` — per-cell expression of `selGene` (Float32Array in `geneExprCache[selGene]`, normalized to [0,1]).
- `mode="gene"` with `_geneMode="coloc"` — two-channel composition of `_colocR` + `_colocG` (red/green) via `_colocCompose`.
- `mode="pathway"` — per-cell pathway score = mean of overlap-gene normalized expressions for `_pathwayColorId`. Score arrays live in `_pathwayScoreCache[pid]` (Float32Array, length n_cells, [0,1]). Computed on demand by `_pathwayComputeScore` from already-cached genes; recomputed lazily on session restore. Per-view snapshot fields: `pathwayColorId`, `pathwayColorName`. Renderer treats it identically to single-gene (same gradient + ghost + size→expression pipeline) — `useGene` in `recolor()` / `_recolorView` is overloaded to mean "use scalar coloring (gene OR pathway)" so the inner loop didn't grow new branches.
- HUD info overlay shows current color mode in a `Color` row per half (`_updateHudColor` reads each split's snapshot, NOT globals — globals are scratch state during the cross-half render loop).

#### Adding a new way to color the main viewer
**Any new "color mode" for the main viewer (gene, coloc, pathway score, signature score, density, …) MUST be wired into split-screen from day one. Half-finished modes that work in single-view but break in split-screen are a recurring failure mode here.** Concretely, when adding a new mode:

1. **Pick the mode key.** New top-level modes extend `mode` (currently `"cluster"` / `"gene"`); sub-modes nest like `_geneMode` (`"single"` / `"coloc"`). Pick the level that matches how the user thinks about it — pathway scoring is a top-level alternative to gene/cluster, so it should probably get its own `mode` value rather than hide inside `_geneMode`.
2. **Add EVERY mode-defining field to the per-view snapshot.** That means `_views[i]` (around the `_views.push({…})` site) plus `_viewSwitch` (which writes globals back from a snapshot on activation). For pathway coloring this would be e.g. `pathwayId`, `pathwayScoreMethod`, `pathwayScoreCache`. If a field is consulted by the renderer and is NOT in the snapshot, the other split half will silently render with the wrong state.
3. **Teach `_recolorView(vIdx)` the new mode.** The `recolor()` driver already swaps globals around the per-view loop; `_recolorView` just needs the same `if(mode === "...")` branch as the on-screen `recolor()` path. If the renderer needs per-cell data (gene expression, score array), make sure the data is fetched / cached for the half being rendered BEFORE `_recolorView` runs. The coloc save+restore pattern in `recolor()` is the model — both channels are pre-fetched, then each view renders from cache.
4. **Sync UI on active-half switch.** `_viewSwitch` calls `_syncColorModeUIToGlobals` (or equivalent) so panels, indicators, and badges reflect the active half. Any new UI affordance for the new mode (a gene-panel-style picker, a mode toggle, a chip in the info overlay) needs a sync line here.
5. **Update the info-overlay color-mode chip.** It must read from the active half's snapshot in split mode, NOT from globals — globals are scratch during the cross-half render loop.
6. **SVG export.** `_exportComputeRenderState` reads per-view snapshots; the side-by-side path runs each half through it independently. New mode → branch in `_exportComputeRenderState`, branch in the per-view SVG renderer, and tag the filename / title with whatever the user is looking at.
7. **Session restore.** `restoreState` (sessionStorage) and `_applySessionState` (server-side auto-session) must serialize and restore the new fields, otherwise a reload silently drops color state. Each `_views[i]` round-trips through these, so whatever field you added in step 2 is the field to add here too.

A reasonable sanity test before marking the mode "done": split the screen, set the mode differently in each half, switch active halves a few times, hit Recolor, reload the page, then SVG-export the side-by-side. If any step shows the wrong colors or loses state, you're missing one of the seven points above.

### Auto-quit (heartbeat)
- Renamed from "auto-unload" but variables/settings keep the `auto_unload`
  name for session-storage compatibility. Default: OFF (opt-in), 20-second
  foreground timeout and 180-second hidden-tab timeout.
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
- GSEA preranked (ORA is shipped — Pathways panel; GSEA reuses the same panel, server endpoint stub still TODO)
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
