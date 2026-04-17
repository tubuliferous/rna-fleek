# RNA-FLEEK

**Fast Lightweight Expression Exploration Kit**

A browser-based single-cell RNA-seq visualizer. Two files — a Python server and an HTML frontend — give you an interactive 3D/2D UMAP explorer with gene expression overlays, selection groups, differential expression analysis, and more. No build tools, no npm, no Docker.

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

## Download

Pre-built desktop apps — no Python installation required:

| Platform | Download |
|----------|----------|
| **Windows** | [RNA-FLEEK-Windows.zip](https://github.com/tubuliferous/rna-fleek/releases/latest/download/RNA-FLEEK-Windows.zip) |
| **macOS** | [RNA-FLEEK-Mac.zip](https://github.com/tubuliferous/rna-fleek/releases/latest/download/RNA-FLEEK-Mac.zip) |
| **Linux** | [RNA-FLEEK-Linux.tar.gz](https://github.com/tubuliferous/rna-fleek/releases/latest/download/RNA-FLEEK-Linux.tar.gz) |

> **Note:** These links point to the latest release. If no release has been published yet, download from the [Actions artifacts](https://github.com/tubuliferous/rna-fleek/actions) page instead.

### macOS: first-launch security warning

The Mac build isn't code-signed or notarized, so Gatekeeper will block it on first launch with one of two messages:

- *"RNA-FLEEK can't be opened because Apple cannot check it for malicious software"* — right-click (or Control-click) the app and choose **Open**, then click **Open** in the confirmation dialog. This only needs to be done once.
- *"RNA-FLEEK is damaged and can't be opened"* — this is the quarantine flag added by Safari/Chrome. Remove it with:

  ```bash
  xattr -d com.apple.quarantine /Applications/RNA-FLEEK.app
  ```

  (Adjust the path to wherever you placed the app.)

If the right-click trick doesn't produce an "Open" option, try **System Settings → Privacy & Security**, scroll to the bottom, and click **Open Anyway** next to the blocked-app notice.

Or install from source (see Quick Start below).

---

## Quick Start

```bash
pip install scanpy anndata numpy scipy

python fleek_server.py my_dataset.h5ad
```

Open [http://localhost:8080](http://localhost:8080) in your browser. That's it.

To start without a dataset and upload via the browser:

```bash
python fleek_server.py
```

## Features

**Visualization**
- 3D and 2D UMAP with WebGL point cloud (Three.js)
- Perspective depth in 3D — near points appear larger
- Fly-through navigation: orbit, pan, zoom, and translate through the point cloud
- Orbit pivot cross for spatial orientation (depth-tested against point cloud)
- Intro orbit animation on first load
- Dark and light themes with theme-appropriate colormaps

**Gene Expression**
- On-demand gene loading from the server, cached client-side
- Viridis colormap (dark mode) / Blues sequential (light mode)
- Top Variable Genes shown as quick-access chips
- Background CSC index build for instant gene lookups on large datasets
- Gene search across all genes in the dataset

**Selection & Analysis**
- Lasso selection, single-cell pick, cell-type bulk selection
- Multiple named selection groups (non-overlapping, color-coded)
- DEG analysis: Wilcoxon rank-sum or Welch's t-test
- Volcano plot with hover stats and click-to-load-gene
- Cohen's d effect size computed directly on sparse matrices
- Results cached per test type + group combination
- Export cell subsets as .h5ad files

**Cell Types**
- Searchable legend with per-type cell counts
- Sort alphabetically or by count
- Click to toggle, Option-click to solo/unsolo
- Shift-click to select entire type into active group

**Performance**
- Handles 1.6M+ cells on a laptop
- QuickMap: subsample-and-project UMAP for large datasets
- Disk-backed mode for datasets larger than RAM
- Batch PCA projection avoids materializing dense 60k-gene matrices
- DEG pre-filters to expressed genes (3–6× speedup)
- Background threading — visualization is interactive immediately

## Requirements

**Server:**

```
numpy
scipy
scanpy
anndata
```

Optional, for CELLxGENE Census downloads:

```
cellxgene-census
```

Install everything:

```bash
pip install scanpy anndata numpy scipy
```

**Client:** Any modern browser with WebGL (Chrome, Firefox, Safari, Edge).

## Usage

### Server

```bash
# Load a specific file
python fleek_server.py my_data.h5ad

# Custom port
python fleek_server.py my_data.h5ad --port 9090

# QuickMap for large datasets (fit UMAP on 50k subsample, project the rest)
python fleek_server.py huge_dataset.h5ad --fast-umap --fast-umap-n 50000

# Disk-backed mode (keeps expression matrix on disk)
python fleek_server.py huge_dataset.h5ad --backed on

# Auto disk-backed (uses backed mode if file > 40% of available RAM)
python fleek_server.py huge_dataset.h5ad --backed auto

# Start empty, upload via browser
python fleek_server.py
```

### Controls

#### 3D View

| Action | Input |
|--------|-------|
| Rotate | Left-drag |
| Pan | Ctrl+drag or right-drag |
| Fly through | Option+drag or middle-drag |
| Zoom | Scroll wheel |

#### 2D View

| Action | Input |
|--------|-------|
| Pan | Drag |
| Zoom | Scroll wheel |

#### Selection

| Action | Input |
|--------|-------|
| Lasso select | Shift+drag |
| Lasso deselect | Shift+Option+drag |
| Toggle single cell | Shift+click |
| Toggle cell type | Shift+double-click on UMAP, or Shift+click in Cell Types |
| Solo visibility | Option+click dot (click again to restore) |
| Undo / Redo | Cmd+[ / Cmd+] |

### Downloading Data from CELLxGENE Census

```bash
pip install cellxgene-census

# Download blood cells (chunked, resumable)
python download_census.py --tissue blood

# Limit to 500k cells
python download_census.py --tissue liver --max-cells 500000

# Custom chunk size
python download_census.py --tissue brain --chunk-size 50000
```

Available tissues: blood, lung, brain, liver, heart, kidney, and others from the Census `tissue_general` field.

## Architecture

```
fleek_server.py    Python HTTP server — loads h5ad, computes UMAP, serves API
fleek.html         Single-file frontend — Three.js WebGL, all JS/CSS inline
download_census.py Chunked CELLxGENE Census downloader with resume support
preprocess_census.py  Preprocessor for standalone .bin files
```

**Server endpoints:**

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Serves fleek.html |
| `/api/init` | GET | Binary SCRN payload (coordinates, clusters, metadata) |
| `/api/gene?name=X` | GET | Float32 expression array for one gene |
| `/api/search?q=X` | GET | Gene name autocomplete |
| `/api/deg` | POST | Differential expression analysis |
| `/api/export` | POST | Download cell subset as .h5ad |
| `/api/upload` | POST | Upload new dataset |
| `/api/progress` | GET | Background task progress |
| `/api/load` | POST | Load server-side file by path |

**Data format:** The init payload uses a custom binary format ("SCRN") — 4-byte magic, uint32 header length, JSON header (padded to 4-byte boundary), followed by float32 coordinate arrays and int32 cluster IDs. Gene expression is sent as raw float32 arrays on demand.

**UMAP caching:** Computed embeddings are saved to `<filename>.umap_cache.npz` alongside the h5ad file. Subsequent loads skip UMAP computation entirely.

## Methods

**UMAP** — Computed via scanpy with default parameters (n_neighbors=15, min_dist=0.5). Both 2D and 3D embeddings are generated. QuickMap fits on a stratified subsample preserving cluster proportions, then projects remaining cells via batch PCA and UMAP transform.

**Differential expression** — Uses scanpy's `rank_genes_groups`. Genes expressed in <3 cells across both groups are filtered before testing. P-values adjusted via Benjamini-Hochberg. If the data appears to contain raw counts, `normalize_total` and `log1p` are auto-applied to the test subset.

**Cohen's d** — Difference in means divided by pooled standard deviation with Bessel correction. Computed directly on sparse matrices using E[X²]−E[X]² to avoid materializing dense arrays.

**Expression normalization** — For visualization, expression is normalized to [0, 1] using the 98th percentile of expressing cells as the ceiling. This prevents outliers from washing out the color range.

**Pseudoreplication caveat** — Cells from the same donor are not independent. P-values from single-cell DEG tests are anti-conservative. For publication results, use pseudobulk approaches.

## Performance Notes

| Dataset size | RAM usage | UMAP time | CSC index |
|---|---|---|---|
| ~100k cells | ~2 GB | ~30s | ~5s |
| ~500k cells | ~6 GB | ~3 min | ~30s |
| ~1.6M cells | ~9 GB (batch PCA) | ~15 min | ~5 min |

QuickMap at 50k fit cells reduces UMAP time to ~2 minutes regardless of dataset size.

## State Persistence

All visualization state is saved to sessionStorage and survives page refreshes and browser zoom changes: camera position, selection groups, gene expression mode, visibility toggles, DEG results, and settings. Undo/redo history (up to 100 steps) is maintained within a session.

## License

MIT
