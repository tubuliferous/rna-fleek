#!/usr/bin/env python3
"""
RNA-FLEEK Server — Load h5ad and serve to browser visualizer.

Requirements:
    pip install scanpy anndata numpy scipy

Usage:
    python fleek_server.py liver_census.h5ad
    python fleek_server.py liver_census.h5ad --port 8080 --max-cells 500000

Then open: http://localhost:8080
"""

import argparse
import io
import json
import math
import os
import struct
import sys
import time
import threading
import warnings
from http.server import HTTPServer, SimpleHTTPRequestHandler
from socketserver import ThreadingMixIn
from pathlib import Path

# Suppress noisy pandas/scanpy warnings
warnings.filterwarnings("ignore", message="DataFrame is highly fragmented")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Suffix used.*to deduplicate")
from urllib.parse import urlparse, parse_qs

import numpy as np

# PyInstaller bundles data files into sys._MEIPASS; normal runs use script dir
_BUNDLE_DIR = Path(getattr(sys, '_MEIPASS', Path(__file__).parent))

# Globals set at startup
ADATA = None
UMAP_2D = None
UMAP_3D = None
# Per-embedding provenance: {"umap_2d": "obsm" | "cache_quick" | "cache_full" |
# "computed" | "derived_from_3d" | "unavailable", ...}.
EMBEDDING_SOURCES = {}
PCA_2D = None
PCA_3D = None
PACMAP_2D = None
PACMAP_3D = None
PACMAP_COMPUTING = False
CLUSTER_IDS = None
CLUSTER_NAMES = None
CLUSTER_COL = None    # obs column name used for clustering
CLUSTER_COLORS = []   # hex color per cluster, matching client palette
GENE_NAMES_LIST = None
N_CELLS = 0
HTML_DIR = None
PROGRESS = {"status": "idle", "message": "", "pct": 0}  # for client polling
PROCESSING_LOCK = threading.Lock()
ABORT_REQUESTED = False  # Set True to cancel in-progress loading
ANNOT_STATUS = {"active": False, "method": "", "step": "", "pct": 0}  # For client polling
DEG_SETTINGS = {"fast": True, "max_cells": 50000}  # QuickDEG defaults
PACMAP_SETTINGS = {"fast": True, "max_cells": 50000}  # QuickPaCMAP defaults
PACMAP_FAST = True  # QuickMap (PaCMAP) — subsample-and-project
BACKED = False  # True if adata.X is on disk (backed mode)
X_CSC = None  # CSC copy of expression matrix for fast column access
HVG_NAMES = []  # ordered list of HVG names for the client
HVG_VAR = {}    # gene_name -> variance (for tooltip/list display)
ALL_GENE_VAR = None  # numpy array of variance for all genes (computed in background)
ALL_GENE_VAR_METRIC = "dispersions_norm"  # metric name for client display
GENE_CUTOFF_FLAG = None  # int8 bitmask: 1=not expressed, 2=low mean, 4=high mean, 8=low fraction
GENE_INDEX = {}  # gene_name -> column index
CSC_BUILDING = False  # True while background CSC build is running
CSC_CACHED = False    # True if CSC was loaded from cache
CSC_TIME = 0          # seconds taken to build/load CSC index
LOAD_SETTINGS = {}  # Sent to client so it knows what options were active
LOADED_PATH = ""    # path to currently loaded h5ad, for annotation caching
MARKER_DB = None    # {cell_type: [genes]} — loaded from file or built-in
GO_DB = None        # Gene Ontology lookup: {terms, gene_to_terms, synonyms, ...}
GO_DATASET_GENES = None  # Set of gene names in current dataset (for intersection)
GO_SYNONYM_MAP = None    # synonym -> canonical gene symbol (for dataset matching)
COUNTS_MATRIX = None  # raw counts matrix reference (for pseudo-bulk DEG)
COUNTS_LABEL = None   # where counts were found (e.g. "X", "layers['counts']")
COUNTS_INFO = None    # {"label","median_total","max","fano_median","looks_smoothed","pb_ok","reason"}
try:
    import pydeseq2 as _pydeseq2_probe  # noqa: F401
    PYDESEQ2_AVAILABLE = True
except Exception:
    PYDESEQ2_AVAILABLE = False
OBS_COLS = {}       # {col_name: {"categories": [...], "counts": [...], "codes": np.int16 array}}
_DETECTED_ORGANISM = None  # (name, reason) tuple from _detect_organism()
AUTO_UNLOAD = False   # If True, unload dataset when no client heartbeat
AUTO_UNLOAD_TIMEOUT = 1.0  # seconds without heartbeat before unload
_LAST_HEARTBEAT = 0.0  # time.time() of last heartbeat
_HEARTBEAT_TIMER = None  # threading.Timer for delayed unload

# ── Cache location management ──
# CACHE_MODE: "dataset" = store caches alongside .h5ad file (default)
#             "server"  = store caches in a server-local directory
# Falls back to server dir when dataset directory is not writable.
# Reads check both locations (dataset dir first, then server dir).
CACHE_MODE = "dataset"
CACHE_DIR = None  # Path object, set at startup (defaults to ~/.fleek_cache)

# Matches the JS palette P[] in fleek.html exactly — used to bake colors into exports
_FLEEK_PALETTE = [
    (.91,.30,.24),(.20,.60,.86),(.15,.68,.38),(.95,.61,.16),(.56,.27,.68),
    (.90,.49,.67),(.10,.74,.61),(.93,.84,.19),(.40,.40,.40),(1.,.42,.42),
    (.36,.42,.75),(.30,.82,.22),(.80,.47,.20),(.60,.35,.71),(.16,.50,.73),
    (.85,.17,.36),(.47,.75,.86),(.72,.67,.28),(.55,.22,.32),(.22,.70,.50),
    (.78,.56,.73),(.60,.80,.42),(.35,.30,.55),(.95,.72,.55),
]
def _palette_hex(idx):
    r, g, b = _FLEEK_PALETTE[idx % len(_FLEEK_PALETTE)]
    return f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"

def _json_safe(obj):
    """Recursively convert NaN / ±Inf floats to None so json.dumps produces
    strict JSON that the browser can parse. Numpy scalars are coerced to
    builtins first so isnan/isinf work on them."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    # numpy scalar / array element
    if hasattr(obj, "item") and not isinstance(obj, (str, bytes)):
        try:
            v = obj.item()
            if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                return None
            return v
        except Exception:
            pass
    return obj

def _init_cache_dir():
    """Initialize the server-side cache directory."""
    global CACHE_DIR
    if CACHE_DIR is None:
        CACHE_DIR = Path.home() / ".fleek_cache"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

def _dataset_cache_path(suffix):
    """Cache path alongside the dataset file."""
    if not LOADED_PATH:
        return None
    return Path(LOADED_PATH).with_suffix(suffix)

def _server_cache_path(suffix):
    """Cache path in the server cache directory, namespaced by dataset."""
    if not LOADED_PATH or CACHE_DIR is None:
        return None
    ds_name = Path(LOADED_PATH).stem
    return CACHE_DIR / f"{ds_name}{suffix}"

def _is_writable(path):
    """Check if we can write to the directory containing path."""
    try:
        d = path.parent
        if not d.exists():
            return False
        test = d / f".fleek_write_test_{os.getpid()}"
        test.touch()
        test.unlink()
        return True
    except (OSError, PermissionError):
        return False

def _cache_read_path(suffix):
    """Find an existing cache file. Checks dataset dir first, then server dir.
    Returns (path, writable) or (None, False) if no cache exists."""
    dp = _dataset_cache_path(suffix)
    if dp and dp.exists():
        return dp, _is_writable(dp)
    sp = _server_cache_path(suffix)
    if sp and sp.exists():
        return sp, True
    return None, False

def _cache_write_path(suffix):
    """Determine where to write a new cache file.
    In 'dataset' mode, tries dataset dir first, falls back to server dir.
    In 'server' mode, always uses server dir."""
    _init_cache_dir()
    if CACHE_MODE == "dataset":
        dp = _dataset_cache_path(suffix)
        if dp and _is_writable(dp):
            return dp
        print(f"  Cache: dataset dir not writable, using server cache dir")
    return _server_cache_path(suffix)

# ── Built-in compact marker set for common human cell types ──
_BUILTIN_MARKERS = {
    "T cell": ["CD3D","CD3E","CD3G","CD2","TRAC","CD28"],
    "CD4+ T cell": ["CD3D","CD3E","CD4","IL7R","MAL","LEF1"],
    "CD8+ T cell": ["CD3D","CD3E","CD8A","CD8B","GZMK","NKG7"],
    "Regulatory T cell": ["FOXP3","IL2RA","CTLA4","TIGIT","IKZF2"],
    "NK cell": ["NKG7","GNLY","KLRD1","KLRB1","NCAM1","GZMB","PRF1"],
    "B cell": ["CD79A","CD79B","MS4A1","CD19","PAX5","BANK1"],
    "Plasma cell": ["JCHAIN","MZB1","SDC1","IGHG1","XBP1","PRDM1"],
    "Monocyte": ["CD14","LYZ","S100A9","S100A8","VCAN","FCN1"],
    "Classical monocyte": ["CD14","LYZ","S100A9","S100A8","VCAN"],
    "Non-classical monocyte": ["FCGR3A","MS4A7","LST1","LILRB2","IFITM3"],
    "Macrophage": ["CD68","CD163","MRC1","MSR1","MARCO","LGMN","C1QA"],
    "Dendritic cell": ["FCER1A","CLEC10A","CD1C","HLA-DQA1","ITGAX"],
    "Plasmacytoid dendritic cell": ["CLEC4C","IL3RA","NRP1","IRF7","TCF4"],
    "Mast cell": ["TPSAB1","TPSB2","CPA3","KIT","HPGDS","HDC"],
    "Neutrophil": ["CSF3R","FCGR3B","CXCR2","S100A12","MMP9"],
    "Eosinophil": ["CLC","RNASE2","PRG2","EPX","CCR3"],
    "Basophil": ["CLC","HDC","GATA2","CPA3","MS4A2"],
    "Platelet": ["PF4","PPBP","GP9","ITGA2B","TUBB1"],
    "Erythrocyte": ["HBA1","HBA2","HBB","HBD","ALAS2","SLC4A1"],
    "Erythroid progenitor": ["GYPA","GYPB","KLF1","TFRC","EPOR"],
    "HSC": ["CD34","THY1","KIT","CRHBP","HLF","AVP","MLLT3"],
    "Megakaryocyte": ["PF4","PPBP","GP9","ITGA2B","TUBB1","ITGB3"],
    "Fibroblast": ["DCN","COL1A1","COL1A2","LUM","PDGFRA","THY1"],
    "Myofibroblast": ["ACTA2","TAGLN","MYL9","COL1A1","COL3A1"],
    "Smooth muscle cell": ["ACTA2","MYH11","TAGLN","CNN1","DES"],
    "Endothelial cell": ["PECAM1","VWF","CDH5","ERG","FLT1","KDR"],
    "Lymphatic endothelial cell": ["PROX1","LYVE1","FLT4","PDPN","CCL21"],
    "Epithelial cell": ["EPCAM","KRT18","KRT19","KRT8","CDH1"],
    "Basal cell": ["KRT5","KRT14","TP63","ITGA6","ITGB4"],
    "Alveolar type 1": ["AGER","PDPN","CAV1","EMP2","HOPX"],
    "Alveolar type 2": ["SFTPC","SFTPB","SFTPA1","ABCA3","NAPSA"],
    "Club cell": ["SCGB1A1","SCGB3A1","CYP2F1","BPIFB1"],
    "Ciliated cell": ["FOXJ1","PIFO","TPPP3","SNTN","CAPS"],
    "Goblet cell": ["MUC5AC","MUC5B","SPDEF","AGR2","TFF3"],
    "Hepatocyte": ["ALB","APOB","HP","TF","TTR","PCK1","CYP3A4"],
    "Cholangiocyte": ["KRT19","KRT7","EPCAM","SOX9","SPP1"],
    "Stellate cell": ["RGS5","ACTA2","DES","PDGFRB","LRAT"],
    "Cardiomyocyte": ["TNNT2","MYH7","MYH6","ACTC1","MYL2"],
    "Pericyte": ["RGS5","PDGFRB","NOTCH3","KCNJ8","ABCC9"],
    "Adipocyte": ["ADIPOQ","LEP","PLIN1","FABP4","LPL"],
    "Neuron": ["RBFOX3","MAP2","SYT1","SNAP25","NEFL","NEFM"],
    "Astrocyte": ["GFAP","AQP4","SLC1A3","SLC1A2","ALDH1L1"],
    "Oligodendrocyte": ["MBP","PLP1","MOG","MAG","MOBP"],
    "Microglia": ["CX3CR1","P2RY12","TMEM119","CSF1R","AIF1"],
    "Schwann cell": ["MPZ","PRX","PMP22","SOX10","S100B"],
    "Podocyte": ["NPHS1","NPHS2","PODXL","WT1","SYNPO"],
    "Proximal tubule cell": ["SLC34A1","LRP2","CUBN","SLC22A6","ALDOB"],
    "Melanocyte": ["PMEL","MLANA","TYR","TYRP1","DCT","MITF"],
    "Keratinocyte": ["KRT1","KRT10","KRT14","KRT5","IVL","LOR"],
    "Chondrocyte": ["COL2A1","ACAN","SOX9","COL9A1","COMP"],
    "Osteoblast": ["BGLAP","RUNX2","SP7","COL1A1","ALPL"],
    "Osteoclast": ["ACP5","CTSK","MMP9","OSCAR","DCSTAMP"],
}

def load_marker_db(marker_path=None):
    """Load marker database from JSON file or use built-in defaults."""
    global MARKER_DB
    
    # Try loading full CellMarker2 database
    paths_to_try = []
    if marker_path:
        paths_to_try.append(marker_path)
    # Look in bundle dir, next to server script, and in working directory
    script_dir = Path(__file__).parent
    paths_to_try.extend([
        _BUNDLE_DIR / "rna_fleek" / "cell_markers.json",
        _BUNDLE_DIR / "cell_markers.json",
        script_dir / "data" / "cell_markers.json",
        script_dir / "cell_markers.json",
        Path("data") / "cell_markers.json",
        Path("cell_markers.json"),
    ])
    
    for p in paths_to_try:
        if Path(p).exists():
            try:
                with open(p) as f:
                    data = json.load(f)
                # Use the global (cross-tissue) human markers
                if "_global" in data and "Human" in data["_global"]:
                    MARKER_DB = data["_global"]["Human"]
                    print(f"  Loaded CellMarker2 database: {len(MARKER_DB)} cell types from {p}")
                    return
            except Exception as e:
                print(f"  Warning: Failed to load {p}: {e}")
    
    # Fall back to built-in
    MARKER_DB = dict(_BUILTIN_MARKERS)
    print(f"  Using built-in marker database: {len(MARKER_DB)} cell types")
    print(f"  (Run download_markers.py for the full CellMarker2 database)")


def _detect_organism():
    """Detect organism from gene names using marker genes + casing heuristic.
    Returns (organism, reason) tuple."""
    if GENE_NAMES_LIST is None or len(GENE_NAMES_LIST) == 0:
        return "human", "default (no gene names)"
    genes = set(GENE_NAMES_LIST)
    # Check for species-specific marker genes (case-sensitive)
    human_markers = {"CD3D", "CD3E", "CD8A", "CD4", "MS4A1", "NKG7", "GZMB",
                     "COL1A1", "EPCAM", "PECAM1", "PTPRC", "HBA1", "GAPDH",
                     "ACTB", "MALAT1", "TMSB4X", "FTL", "FTH1", "S100A9"}
    mouse_markers = {"Cd3d", "Cd3e", "Cd8a", "Cd4", "Ms4a1", "Nkg7", "Gzmb",
                     "Col1a1", "Epcam", "Pecam1", "Ptprc", "Hba-a1", "Gapdh",
                     "Actb", "Malat1", "Tmsb4x", "Ftl1", "Fth1", "S100a9"}
    h_hits = genes & human_markers
    m_hits = genes & mouse_markers
    if len(h_hits) > len(m_hits) and len(h_hits) >= 3:
        return "human", f"matched {len(h_hits)} human marker genes: {', '.join(sorted(h_hits)[:6])}"
    if len(m_hits) > len(h_hits) and len(m_hits) >= 3:
        return "mouse", f"matched {len(m_hits)} mouse marker genes: {', '.join(sorted(m_hits)[:6])}"
    # Fallback: casing heuristic on alpha-only gene names
    sample = [g for g in GENE_NAMES_LIST[:200] if g.isalpha() and len(g) > 1]
    if not sample:
        return "human", "default (no alpha gene names to analyze)"
    upper = sum(1 for g in sample if g == g.upper())
    title = sum(1 for g in sample if g[0].isupper() and g[1:] == g[1:].lower())
    n = len(sample)
    if upper > n * 0.7:
        return "human", f"gene name casing: {upper}/{n} UPPERCASE (human convention)"
    if title > n * 0.5:
        return "mouse", f"gene name casing: {title}/{n} Title case (mouse convention)"
    return "human", f"default (ambiguous casing: {upper}/{n} upper, {title}/{n} title)"


def load_go_db():
    """Load Gene Ontology database if available."""
    global GO_DB, GO_DATASET_GENES, GO_SYNONYM_MAP
    GO_DB = None
    GO_DATASET_GENES = None
    GO_SYNONYM_MAP = None

    global _DETECTED_ORGANISM
    organism, _org_reason = _detect_organism()
    _DETECTED_ORGANISM = (organism.capitalize(), _org_reason)
    print(f"  Organism: {_DETECTED_ORGANISM[0]} ({_org_reason})")
    script_dir = Path(__file__).parent
    paths_to_try = [
        _BUNDLE_DIR / "rna_fleek" / f"go_{organism}.json",
        _BUNDLE_DIR / f"go_{organism}.json",
        script_dir / "data" / f"go_{organism}.json",
        script_dir / f"go_{organism}.json",
        Path("data") / f"go_{organism}.json",
        Path(f"go_{organism}.json"),
    ]

    for p in paths_to_try:
        if Path(p).exists():
            try:
                with open(p) as f:
                    GO_DB = json.load(f)
                print(f"  Loaded GO database: {GO_DB.get('n_terms', '?')} terms, {GO_DB.get('n_genes', '?')} genes ({organism}) from {p}")
                break
            except Exception as e:
                print(f"  Warning: Failed to load GO from {p}: {e}")

    if GO_DB is None:
        print(f"  No GO database found for {organism} (run: python download_go.py --organism {organism})")
        return

    # Build dataset gene set and case-insensitive lookup for robust matching
    if GENE_NAMES_LIST:
        GO_DATASET_GENES = set(GENE_NAMES_LIST)
        # Case-insensitive map: uppercase -> actual dataset gene name
        _gene_upper = {}
        for g in GENE_NAMES_LIST:
            _gene_upper[g.upper()] = g
        # GO->dataset resolver: maps GO gene symbol to dataset gene name
        GO_SYNONYM_MAP = {}
        # First: direct GO gene names -> dataset (case-insensitive)
        go_gene_to_terms = GO_DB.get("gene_to_terms", {})
        for go_gene in go_gene_to_terms:
            if go_gene in GO_DATASET_GENES:
                GO_SYNONYM_MAP[go_gene] = go_gene
            elif go_gene.upper() in _gene_upper:
                GO_SYNONYM_MAP[go_gene] = _gene_upper[go_gene.upper()]
        # Second: synonyms -> dataset (case-insensitive)
        synonyms = GO_DB.get("synonyms", {})
        for syn, canonical in synonyms.items():
            if canonical in GO_SYNONYM_MAP:
                # canonical already resolved — map synonym to same dataset gene
                if syn not in GO_SYNONYM_MAP:
                    GO_SYNONYM_MAP[syn] = GO_SYNONYM_MAP[canonical]
            elif syn.upper() in _gene_upper:
                GO_SYNONYM_MAP[syn] = _gene_upper[syn.upper()]
        n_resolved = sum(1 for v in GO_SYNONYM_MAP.values() if v in GO_DATASET_GENES)
        print(f"  GO matching: {n_resolved} GO genes resolved to dataset genes")


MATRIX_REGISTRY = []
"""List of dicts classifying every matrix present in the loaded h5ad. Populated by
_classify_all_matrices(). Each entry contains:
    path          — how to reference it from a client ("X", "raw.X", "layers[NAME]")
    label         — human-friendly name
    shape         — (n_cells, n_vars) tuple for display
    dtype         — numpy dtype string
    is_integer    — bool (sample appears to be integers)
    max           — max value in 500x500 sample
    median_total  — median per-cell total over ~200 cells
    fano_median   — median Fano across moderately-expressed genes, or None
    kind          — "raw_counts" | "denoised_counts" | "log_normalized" | "scaled" | "unknown"
    notes         — short free-text classification reasoning
"""

MATRIX_ROUTING = {}
"""Per-task matrix routing decisions. Populated by _build_matrix_routing().
Keys are task names; values are:
    path      — matrix path string matching MATRIX_REGISTRY entry
    transform — "none" | "normalize+log1p"
    label     — human-friendly task label
    caveat    — short warning string, or "" if no caveat
"""


def _classify_all_matrices():
    """Walk .X, .raw.X, and all layers; classify each. Populates MATRIX_REGISTRY."""
    global MATRIX_REGISTRY
    MATRIX_REGISTRY = []
    if ADATA is None:
        return

    import scipy.sparse as sp
    candidates = [("X", "adata.X", ADATA.X)]
    if ADATA.raw is not None:
        candidates.append(("raw.X", "adata.raw.X", ADATA.raw.X))
    for name in ADATA.layers:
        candidates.append((f"layers[{name!r}]", f"adata.layers['{name}']", ADATA.layers[name]))

    for path, label, mat in candidates:
        try:
            shape = tuple(mat.shape)
            dtype_s = str(getattr(mat, "dtype", ""))
            sample = mat[:500, :500]
            if sp.issparse(sample):
                sample = sample.toarray()
            elif hasattr(sample, "toarray"):
                sample = sample.toarray()
            sample = np.asarray(sample, dtype=np.float64)
            any_neg = bool((sample < 0).any())
            max_val = float(sample.max()) if sample.size else 0.0
            min_nz = float(sample[sample > 0].min()) if (sample > 0).any() else 0.0
            # Integer check on sample
            is_int = bool(np.allclose(sample[sample >= 0], np.round(sample[sample >= 0])))
            # Median per-cell total across first 200 rows
            row_sums = np.array(mat[:200, :].sum(axis=1)).flatten()
            median_total = float(np.median(row_sums)) if row_sums.size else 0.0

            # Fano on moderately-expressed genes, on a small subsample of cells.
            fano_median = None
            try:
                n_cells = mat.shape[0]
                n_sub = min(3000, n_cells)
                if n_sub > 0:
                    stride = max(1, n_cells // n_sub)
                    row_idx = np.arange(0, n_cells, stride)[:n_sub]
                    sub = mat[row_idx, :]
                    if sp.issparse(sub):
                        mean_arr = np.asarray(sub.mean(axis=0)).ravel()
                        sq_mean = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
                        gene_sum = np.asarray(sub.sum(axis=0)).ravel()
                    else:
                        sub_d = np.asarray(sub, dtype=np.float64)
                        mean_arr = sub_d.mean(axis=0)
                        sq_mean = (sub_d ** 2).mean(axis=0)
                        gene_sum = sub_d.sum(axis=0)
                    var_arr = sq_mean - mean_arr ** 2
                    moderate = (mean_arr > 0) & (var_arr > 0) & (gene_sum > 50)
                    if moderate.sum() >= 50:
                        fano_median = float(np.median(var_arr[moderate] / mean_arr[moderate]))
            except Exception:
                pass

            # Classify
            kind = "unknown"
            notes = ""
            if any_neg:
                kind = "scaled"
                notes = "has negative values — likely z-scored / sc.pp.scale output"
            elif is_int and max_val >= 3 and median_total >= 50:
                if fano_median is not None and fano_median < 1.5:
                    kind = "denoised_counts"
                    notes = (f"integer + Fano≈{fano_median:.2f} (Poisson-like) — denoiser output "
                             f"(e.g. scVI denoised_counts)")
                else:
                    kind = "raw_counts"
                    if fano_median is not None:
                        notes = f"integer UMI counts (Fano≈{fano_median:.2f})"
                    else:
                        notes = "integer counts"
            elif (not is_int) and min_nz > 0 and max_val < 100:
                kind = "log_normalized"
                if max_val < 15:
                    notes = f"decimal values, max≈{max_val:.2f} — looks log-normalized"
                else:
                    notes = f"decimal values, max≈{max_val:.2f} — probably normalized (possibly log)"
            elif (not is_int):
                kind = "log_normalized"
                notes = f"decimal values, max≈{max_val:.2f} — probably normalized"
            else:
                notes = f"integer but unusual (max={max_val:.0f}, median/cell={median_total:.0f})"

            # Source-hint heuristics: layer names often encode provenance that the
            # numerical shape alone can't reveal. E.g. "denoised_norm" looks just
            # like any other log-normalized matrix numerically, but the name tells
            # us the input was model-smoothed data.
            source_hint = ""
            label_lower = label.lower()
            denoised_tags = ["denois", "smooth", "magic", "imput", "scvi", "saver", "dca"]
            raw_tags = ["raw", "umi", "counts", "spliced"]
            if any(tag in label_lower for tag in denoised_tags):
                source_hint = "denoised"
            elif any(tag in label_lower for tag in raw_tags):
                source_hint = "raw"

            # Amend notes / kind labeling when a log-normalized matrix appears
            # to be derived from denoised counts.
            if kind == "log_normalized" and source_hint == "denoised":
                notes = (notes + " — name suggests it's been built from denoised counts, "
                                  "so values are model-smoothed, not raw observations").strip()

            MATRIX_REGISTRY.append({
                "path": path,
                "label": label,
                "shape": list(shape),
                "dtype": dtype_s,
                "is_integer": is_int,
                "max": max_val,
                "median_total": median_total,
                "fano_median": fano_median,
                "kind": kind,
                "source_hint": source_hint,
                "notes": notes,
            })
        except Exception as e:
            MATRIX_REGISTRY.append({
                "path": path, "label": label, "shape": None, "dtype": "?",
                "kind": "unknown", "notes": f"(classification failed: {e})",
            })


def _get_matrix_by_path(path):
    """Resolve a MATRIX_REGISTRY path string back to an actual matrix reference."""
    if ADATA is None:
        return None
    if path == "X":
        return ADATA.X
    if path == "raw.X":
        return ADATA.raw.X if ADATA.raw is not None else None
    if path.startswith("layers[") and path.endswith("]"):
        name = path[len("layers["):-1]
        # strip surrounding quotes if present
        if (name.startswith("'") and name.endswith("'")) or (name.startswith('"') and name.endswith('"')):
            name = name[1:-1]
        if name in ADATA.layers:
            return ADATA.layers[name]
    return None


def _build_matrix_routing():
    """Assign the best matrix for each analysis task based on MATRIX_REGISTRY.

    Preference hierarchy:
      * variability (Fano-based):      raw_counts  →  fallback: log_normalized (different interpretation)
      * single-cell DEG:               raw_counts + normalize+log1p  →  log_normalized as-is  →  denoised as-is (caveat)
      * pseudo-bulk:                   raw_counts only (no fallback)
      * gene expression (visualization): log_normalized  →  raw_counts + normalize+log1p on-the-fly  →  denoised as-is (caveat)
      * embedding fallback (PCA):      log_normalized  →  raw_counts + normalize+log1p  →  denoised (caveat)
    """
    global MATRIX_ROUTING
    MATRIX_ROUTING = {}
    if not MATRIX_REGISTRY:
        return

    # Indexes by kind
    by_kind = {}
    for entry in MATRIX_REGISTRY:
        by_kind.setdefault(entry["kind"], []).append(entry)

    def pick_first(kinds):
        for k in kinds:
            if k in by_kind and by_kind[k]:
                return by_kind[k][0]
        return None

    # When a log-normalized matrix has a denoised-source hint, swap in the
    # task-specific caveat that would apply to denoised_counts directly — the
    # log1p step doesn't rescue the underlying problem (biological overdispersion
    # has already been regressed out by the upstream denoiser).
    denoised_via_lognorm_caveat = {
        "variability": ("Running on log-normalized DENOISED data — the upstream denoiser regressed out "
                        "biological overdispersion, so Fano is artificially tight and trend residuals "
                        "are dominated by residual numerical noise, not biology. Rankings are unreliable."),
        "deg_sc": ("Running on log-normalized DENOISED data — within-group variance has been artificially "
                   "compressed by the upstream denoiser, so Wilcoxon/t-test p-values and Cohen's d are "
                   "anti-conservative. Gene rankings still order roughly correctly; statistical magnitudes "
                   "shouldn't be trusted."),
        "gene_display": ("Display values come from log-normalized DENOISED counts — gradients are smooth "
                         "because the denoiser filled in zeros and regressed away noise. You're seeing the "
                         "model's inference, not raw observations."),
        "embedding": ("PCA built on log-normalized DENOISED data — cluster structure preserved but "
                      "amplitudes / dispersions are compressed."),
    }

    def route(task, task_label, preferred_chain, transforms_per_kind, caveat_per_kind):
        pick = None
        for kind in preferred_chain:
            if kind in by_kind and by_kind[kind]:
                pick = by_kind[kind][0]
                pick_kind = kind
                break
        else:
            pick_kind = None
        if pick is None:
            MATRIX_ROUTING[task] = {
                "path": None, "matrix_label": None, "matrix_kind": None,
                "transform": "none", "label": task_label,
                "caveat": "No suitable matrix available — task disabled.",
            }
            return
        caveat = caveat_per_kind.get(pick_kind, "")
        # If a log-normalized pick is actually built from denoised counts, swap in
        # the harder task-specific caveat rather than just appending a soft note —
        # the log1p step doesn't restore the biological variability the denoiser ate.
        if pick_kind == "log_normalized" and pick.get("source_hint") == "denoised":
            harder = denoised_via_lognorm_caveat.get(task)
            if harder:
                caveat = harder
            else:
                extra = ("Underlying values were denoised before normalization; results inherit "
                         "denoised-data caveats.")
                caveat = (caveat + " " + extra).strip() if caveat else extra
        MATRIX_ROUTING[task] = {
            "path": pick["path"],
            "matrix_label": pick["label"],
            "matrix_kind": pick_kind,
            "matrix_source_hint": pick.get("source_hint", ""),
            "transform": transforms_per_kind.get(pick_kind, "none"),
            "label": task_label,
            "caveat": caveat,
        }

    # Variability (Fano): raw only ideal; log_normalized works approximately
    route("variability", "Gene variability",
          ["raw_counts", "log_normalized", "denoised_counts", "unknown", "scaled"],
          {"raw_counts": "none", "log_normalized": "none", "denoised_counts": "none", "unknown": "none", "scaled": "none"},
          {"log_normalized": "Running on log-normalized data: trend-residual method still works, but Fano interpretation is approximate.",
           "denoised_counts": "Running on denoised (Poisson) data — rankings dominated by residual numerical noise, not biology. Results are unreliable.",
           "scaled": "Running on scaled/z-scored data — Fano math is invalid. Results are meaningless."})
    # Single-cell DEG: any reasonable matrix works; raw gets normalized
    route("deg_sc", "Single-cell DEG",
          ["raw_counts", "log_normalized", "denoised_counts", "unknown", "scaled"],
          {"raw_counts": "normalize+log1p", "log_normalized": "none", "denoised_counts": "normalize+log1p", "unknown": "none", "scaled": "none"},
          {"denoised_counts": "Running on denoised counts — p-values and Cohen's d are anti-conservative (within-group variance artificially low). Rankings are still meaningful.",
           "scaled": "Running on scaled data — fold changes and Cohen's d may be misleading.",
           "unknown": "Matrix type uncertain — interpret p-values cautiously."})
    # Pseudo-bulk: raw only
    route("deg_pb", "Pseudo-bulk DEG",
          ["raw_counts"],
          {"raw_counts": "none"},
          {})
    # Gene expression display: log_normalized is ideal
    route("gene_display", "Gene expression display",
          ["log_normalized", "raw_counts", "denoised_counts", "unknown", "scaled"],
          {"raw_counts": "normalize+log1p", "log_normalized": "none", "denoised_counts": "normalize+log1p", "unknown": "none", "scaled": "none"},
          {"raw_counts": "Normalized + log1p on-the-fly for display.",
           "denoised_counts": "Denoised integer counts normalized + log1p'd on-the-fly. Gradients may be artificially smooth.",
           "scaled": "Scaled matrix may show negative values — colormap interpretation awkward."})
    # Embedding input (used only when we need to compute PCA/UMAP from scratch)
    route("embedding", "Embedding (PCA fallback)",
          ["log_normalized", "raw_counts", "denoised_counts", "unknown", "scaled"],
          {"raw_counts": "normalize+log1p", "log_normalized": "none", "denoised_counts": "normalize+log1p", "unknown": "none", "scaled": "none"},
          {"raw_counts": "Normalized + log1p before PCA.",
           "denoised_counts": "PCA on denoised data — cluster structure preserved but amplitudes off.",
           "scaled": "PCA on already-scaled data may over-weight low-variance genes."})


def detect_counts():
    """Classify every matrix in the loaded adata, then pick the best raw-counts
    entry for pseudo-bulk. MATRIX_REGISTRY / MATRIX_ROUTING get populated as a
    side-effect so every analysis task can share the same hierarchy.
    """
    global COUNTS_MATRIX, COUNTS_LABEL, COUNTS_INFO
    COUNTS_MATRIX = None
    COUNTS_LABEL = None
    COUNTS_INFO = None
    if ADATA is None:
        return

    import scipy.sparse as sp

    # Full classification + task routing in one pass.
    _classify_all_matrices()
    _build_matrix_routing()

    # Prefer matrices classified as "raw_counts"; among those, lean on layer-name
    # hints (counts / raw_counts / umi / spliced) and then the highest Fano.
    hint_names = {"counts", "raw_counts", "raw", "umi", "spliced"}
    raw_candidates = [e for e in MATRIX_REGISTRY if e["kind"] == "raw_counts"]
    def _score(entry):
        s = 0.0
        label_low = entry["label"].lower()
        for hint in hint_names:
            if hint in label_low:
                s += 10.0
                break
        if entry.get("fano_median") is not None:
            s += float(entry["fano_median"])
        s += 0.1 * np.log1p(entry.get("max", 0))
        return s
    raw_candidates.sort(key=_score, reverse=True)

    best = None
    if raw_candidates:
        best = (raw_candidates[0]["path"], _get_matrix_by_path(raw_candidates[0]["path"]))

    # Fallback: honor the old heuristic (integer + reasonable depth + no Fano gate)
    # in case the classifier drops everything into denoised by mistake.
    if best is None:
        for entry in MATRIX_REGISTRY:
            if entry["kind"] in ("denoised_counts",) and entry.get("is_integer") and entry.get("max", 0) >= 3 and entry.get("median_total", 0) >= 50:
                best = (entry["path"], _get_matrix_by_path(entry["path"]))
                break

    if best:
        COUNTS_LABEL, COUNTS_MATRIX = best
        # Quick stats
        row_sums = np.array(COUNTS_MATRIX[:1000, :].sum(axis=1)).flatten()
        median_total = float(np.median(row_sums))
        max_sample = COUNTS_MATRIX[:500, :500]
        max_val = float(np.asarray(max_sample.toarray() if hasattr(max_sample, 'toarray') else max_sample).max())

        # Fano-factor check to detect smoothed/denoised integer counts.
        # Sample up to 10k cells × 500 well-expressed genes, compute var/mean per gene.
        fano_median = None
        looks_smoothed = False
        try:
            n_cells_full = COUNTS_MATRIX.shape[0]
            n_sample = min(10000, n_cells_full)
            # deterministic stride-based sample (no RNG dependency)
            stride = max(1, n_cells_full // n_sample)
            row_idx = np.arange(0, n_cells_full, stride)[:n_sample]
            sub = COUNTS_MATRIX[row_idx, :]
            gene_sums = np.asarray(sub.sum(axis=0)).flatten()
            good = np.where(gene_sums > 50)[0]
            if len(good) >= 50:
                # limit to 500 genes to bound cost
                if len(good) > 500:
                    g_stride = max(1, len(good) // 500)
                    good = good[::g_stride][:500]
                block = sub[:, good]
                if sp.issparse(block):
                    block = block.toarray()
                block = np.asarray(block, dtype=np.float64)
                means = block.mean(axis=0)
                vars_ = block.var(axis=0)
                fano = vars_ / np.maximum(means, 1e-9)
                fano_median = float(np.median(fano))
                # Empirical threshold: real UMI scRNA-seq data has median Fano
                # well above 1.5 because of biological overdispersion + cell
                # type mixing. Poisson-sampled outputs sit near 1.0.
                looks_smoothed = fano_median < 1.5
        except Exception as e:
            print(f"    (Fano check skipped: {e})")

        pb_ok = not looks_smoothed
        reason = ""
        if looks_smoothed:
            reason = (f"median Fano={fano_median:.2f} suggests smoothed/denoised integer counts "
                      f"(e.g. scVI denoised_counts, MAGIC). Pseudo-bulk disabled — NB dispersion "
                      f"estimates collapse on Poisson data, yielding anti-conservative p-values.")

        COUNTS_INFO = {
            "label": COUNTS_LABEL,
            "median_total": median_total,
            "max": max_val,
            "fano_median": fano_median,
            "looks_smoothed": looks_smoothed,
            "pb_ok": pb_ok,
            "reason": reason,
        }

        print(f"  Counts detected in: adata.{COUNTS_LABEL}")
        print(f"    Median counts/cell: {median_total:.0f}, max in sample: {max_val:.0f}"
              f"{', median Fano: ' + f'{fano_median:.2f}' if fano_median is not None else ''}")
        if looks_smoothed:
            print(f"    ⚠ {reason}")
            # Matrix stays available for display, but PB is gated on pb_ok.
    else:
        print("  ⚠ No raw counts detected — pseudo-bulk DEG will be unavailable")
        COUNTS_INFO = None


def run_pseudobulk_deg(group_a_cells, group_b_cells, replicate_col, group_a_name="Group A", group_b_name="Group B", method_pref="auto"):
    """Run pseudo-bulk DEG using pyDESeq2 or fallback to scipy t-test.

    Conditions are defined by selection groups (cell index arrays).
    Replicates are defined by an obs column (e.g. 'donor', 'batch').
    Counts are aggregated per (group × replicate) combination.

    Args:
        group_a_cells: list of cell indices for condition A
        group_b_cells: list of cell indices for condition B
        replicate_col: obs column defining biological replicates
        group_a_name: display name for group A
        group_b_name: display name for group B
        method_pref: "auto" (DESeq2 then fall back), "pydeseq2" (require DESeq2,
                     error if unavailable), or "ttest" (skip DESeq2 entirely and
                     use Welch's t-test on log2-CPM).
    Returns:
        dict with: genes, log2fc, padj, method, n_samples_a, n_samples_b, ...
    """
    import scipy.sparse as sp
    global ABORT_REQUESTED
    # Reset any stale abort flag from a prior cancel so this run can check it fresh.
    ABORT_REQUESTED = False

    if ADATA is None or COUNTS_MATRIX is None:
        raise ValueError("No dataset or counts matrix available")
    if COUNTS_INFO and not COUNTS_INFO.get("pb_ok", True):
        raise ValueError(
            "Pseudo-bulk disabled: " + (COUNTS_INFO.get("reason")
            or "detected counts matrix does not look like raw UMI data"))
    if len(group_a_cells) == 0:
        raise ValueError(f"{group_a_name} has no cells")
    if len(group_b_cells) == 0:
        raise ValueError(f"{group_b_name} has no cells")

    obs = ADATA.obs
    if replicate_col not in obs.columns:
        raise ValueError(f"Column '{replicate_col}' not found in obs")

    # Get replicate labels for each group's cells
    idx_a = np.array(group_a_cells, dtype=int)
    idx_b = np.array(group_b_cells, dtype=int)
    reps_a = obs.iloc[idx_a][replicate_col].astype(str).values
    reps_b = obs.iloc[idx_b][replicate_col].astype(str).values
    unique_a = sorted(set(reps_a))
    unique_b = sorted(set(reps_b))

    print(f"  Pseudo-bulk: {group_a_name} ({len(unique_a)} replicates, {len(idx_a)} cells) vs {group_b_name} ({len(unique_b)} replicates, {len(idx_b)} cells)")

    if len(unique_a) < 2:
        raise ValueError(f"Need ≥2 replicates per group. {group_a_name} has cells in {len(unique_a)} replicate(s) of '{replicate_col}'.")
    if len(unique_b) < 2:
        raise ValueError(f"Need ≥2 replicates per group. {group_b_name} has cells in {len(unique_b)} replicate(s) of '{replicate_col}'.")

    # Aggregate counts per (group × replicate)
    n_genes = ADATA.n_vars
    all_samples = []
    all_conditions = []
    COND_A = "A"
    COND_B = "B"

    for rep_id in unique_a:
        _check_abort()
        cell_idx = idx_a[reps_a == rep_id]
        if len(cell_idx) == 0:
            continue
        counts = COUNTS_MATRIX[cell_idx, :]
        agg = np.asarray(counts.sum(axis=0)).flatten()
        all_samples.append(agg)
        all_conditions.append(COND_A)

    for rep_id in unique_b:
        _check_abort()
        cell_idx = idx_b[reps_b == rep_id]
        if len(cell_idx) == 0:
            continue
        counts = COUNTS_MATRIX[cell_idx, :]
        agg = np.asarray(counts.sum(axis=0)).flatten()
        all_samples.append(agg)
        all_conditions.append(COND_B)

    _check_abort()

    count_matrix = np.vstack(all_samples)  # samples × genes
    conditions = np.array(all_conditions)
    n_a = sum(1 for c in conditions if c == COND_A)
    n_b = sum(1 for c in conditions if c == COND_B)

    print(f"  Pseudo-bulk matrix: {count_matrix.shape[0]} samples × {count_matrix.shape[1]} genes")

    # Filter low-count genes (at least 10 total counts across all samples)
    gene_totals = count_matrix.sum(axis=0)
    keep = gene_totals >= 10
    n_kept = keep.sum()
    print(f"  {n_kept} genes pass minimum count filter (≥10 total)")

    gene_names = np.array(ADATA.var_names)

    # Try pyDESeq2 first, fallback to scipy if it's unavailable or fails at runtime.
    # Runtime failures we've seen:
    #   - joblib/loky workers killed (exit code 2) inside PyInstaller-bundled apps
    #     because the frozen executable can't relaunch itself as a worker. Passing
    #     n_cpus=1 keeps everything in-process and sidesteps that.
    #   - Rare numerical corner cases (singular design, all-zero genes after filter)
    #     that raise from DESeq2's NB fit. We fall back to the t-test on CPM rather
    #     than giving the user a hard error.
    method = "pydeseq2"
    _pd_err = None
    # User explicitly asked for the t-test path — skip DESeq2.
    if method_pref == "ttest":
        _pd_err = RuntimeError("user selected t-test (CPM)")
        print("  Pseudo-bulk: user requested t-test (CPM) — bypassing DESeq2")
    try:
        if method_pref == "ttest":
            raise _pd_err  # short-circuit into the fallback block below
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        import pandas as pd

        count_df = pd.DataFrame(count_matrix[:, keep], columns=gene_names[keep])
        meta_df = pd.DataFrame({"condition": conditions}, index=[f"sample_{i}" for i in range(len(conditions))])
        count_df.index = meta_df.index

        _check_abort()  # last chance before entering the long pyDESeq2 fit
        try:
            dds = DeseqDataSet(counts=count_df, metadata=meta_df, design="~condition", n_cpus=1)
        except TypeError:
            # Older pyDESeq2 versions used a different kwarg name (`inference`/`n_processes`).
            dds = DeseqDataSet(counts=count_df, metadata=meta_df, design="~condition")
        dds.deseq2()

        _check_abort()
        stat = DeseqStats(dds, contrast=["condition", COND_A, COND_B])
        stat.summary()

        results_df = stat.results_df
        genes_out = results_df.index.tolist()
        log2fc = results_df["log2FoldChange"].fillna(0).tolist()
        padj = results_df["padj"].fillna(1).tolist()

        print(f"  pyDESeq2 complete: {len(genes_out)} genes tested")

    except AbortError:
        # Cancellation — let the handler see it, not the fallback.
        raise
    except (ImportError, Exception) as _exc:
        _pd_err = _exc

    if _pd_err is not None:
        # Respect the "require DESeq2" preference: surface the real error
        # instead of silently falling back. (user-chosen ttest is exempted —
        # it's expected to follow the fallback code path.)
        if method_pref == "pydeseq2" and not isinstance(_pd_err, RuntimeError):
            raise ValueError(f"pyDESeq2 required but failed: {_pd_err}")
        method = "ttest"
        if method_pref == "ttest":
            pass  # already printed above
        elif isinstance(_pd_err, ImportError):
            # ImportError can mean either "pyDESeq2 not installed" OR "installed
            # but a dep failed to import" (e.g. formulaic / anndata2ri / a
            # version-mismatched numpy). Print the underlying reason so users
            # can diagnose without guessing.
            print(f"  pyDESeq2 import failed ({_pd_err}) — falling back to t-test on CPM")
            print(f"    Check: python -c 'from pydeseq2.dds import DeseqDataSet'")
        else:
            print(f"  pyDESeq2 failed ({type(_pd_err).__name__}: {_pd_err}) — falling back to t-test on CPM")
        try:
            from scipy import stats

            lib_sizes = count_matrix.sum(axis=1, keepdims=True)
            lib_sizes[lib_sizes == 0] = 1
            cpm = count_matrix / lib_sizes * 1e6
            log_cpm = np.log2(cpm + 1)

            mask_cond_a = conditions == COND_A
            mask_cond_b = conditions == COND_B

            genes_out = []
            log2fc = []
            pvals = []

            for gi in range(n_genes):
                if not keep[gi]:
                    continue
                # Check abort every 500 genes — cheap enough to skip every iteration.
                if (gi & 0x1FF) == 0:
                    _check_abort()
                vals_a = log_cpm[mask_cond_a, gi]
                vals_b = log_cpm[mask_cond_b, gi]
                fc = vals_a.mean() - vals_b.mean()
                try:
                    _, pv = stats.ttest_ind(vals_a, vals_b, equal_var=False)
                except Exception:
                    pv = 1.0
                if np.isnan(pv):
                    pv = 1.0
                genes_out.append(gene_names[gi])
                log2fc.append(float(fc))
                pvals.append(float(pv))

            from statsmodels.stats.multitest import multipletests
            _, padj_arr, _, _ = multipletests(pvals, method="fdr_bh")
            padj = padj_arr.tolist()

            print(f"  t-test on CPM complete: {len(genes_out)} genes tested")

        except AbortError:
            # Cancellation — propagate unwrapped.
            raise
        except Exception as e:
            raise ValueError(f"DEG fallback (t-test on CPM) failed: {e}")

    return {
        "genes": genes_out,
        "log2fc": log2fc,
        "padj": padj,
        "method": method,
        "n_samples_a": n_a,
        "n_samples_b": n_b,
        "group_a_name": group_a_name,
        "group_b_name": group_b_name,
        "replicate_col": replicate_col,
        "n_genes_tested": len(genes_out),
    }


def _remap_cached_annotations(cached):
    """Remap a cached annotation payload so its results[] use the CURRENTLY loaded
    dataset's cluster ordering. When a subset is loaded from a parent's cache, the
    subset's cluster_ids may shift (because unused categories were dropped) even
    though the cluster names still match. Without remapping, the client would
    assign each cluster's predictions to the wrong position.

    Returns a new dict safe to send to the client. No-op if cache already matches.
    """
    if not cached or not isinstance(cached, dict):
        return cached
    cached_names = cached.get("cluster_names")
    if cached_names == CLUSTER_NAMES and cached.get("n_cells") == N_CELLS:
        return cached  # exact match — nothing to do
    name_to_result = {}
    for r in cached.get("results", []) or []:
        if isinstance(r, dict):
            cn = r.get("cluster_name")
            if cn is not None:
                name_to_result[cn] = r
    if not name_to_result:
        return cached  # unexpected format; leave alone
    new_results = []
    for cid, cname in enumerate(CLUSTER_NAMES or []):
        src = name_to_result.get(cname)
        if src:
            new_r = dict(src)
            new_r["cluster_id"] = cid
            new_r["cluster_name"] = cname
            try:
                new_r["n_cells"] = int(np.sum(CLUSTER_IDS == cid)) if CLUSTER_IDS is not None else src.get("n_cells", 0)
            except Exception:
                pass
            new_results.append(new_r)
        else:
            new_results.append({
                "cluster_id": cid, "cluster_name": cname,
                "n_cells": int(np.sum(CLUSTER_IDS == cid)) if CLUSTER_IDS is not None else 0,
                "top_genes": [], "predictions": []
            })
    out = dict(cached)
    out["results"] = new_results
    out["cluster_names"] = list(CLUSTER_NAMES) if CLUSTER_NAMES else []
    out["n_cells"] = N_CELLS
    out["remapped_from_parent"] = True
    return out


def annotate_clusters(top_n=50, test="wilcoxon"):
    """Score clusters against marker database using shared DEG results."""
    from scipy.stats import fisher_exact

    if ADATA is None or MARKER_DB is None:
        return {"error": "No data or marker database loaded"}
    
    if not MARKER_DB:
        load_marker_db()
    if GO_DB is None:
        load_go_db()
    detect_counts()

    # Check annotation cache (with parent stem fallback for subsets)
    _parent_stem = ADATA.uns.get("fleek_parent_stem") if ADATA is not None else None
    cache_path, _cw = _cache_read_path(".annot_markers.json")
    if not cache_path and _parent_stem and CACHE_DIR:
        _pp = CACHE_DIR / f"{_parent_stem}.annot_markers.json"
        if _pp.exists():
            cache_path = _pp
    if cache_path:
        try:
            with open(cache_path) as f:
                cached_annot = json.load(f)
            _cn = cached_annot.get("cluster_names")
            _valid = (_cn == CLUSTER_NAMES and cached_annot.get("n_cells") == N_CELLS) or (_cn and set(CLUSTER_NAMES).issubset(set(_cn)))
            if _valid:
                cached_annot = _remap_cached_annotations(cached_annot)
                print(f"  Marker annotation loaded from cache" + (" (remapped for subset)" if cached_annot.get("remapped_from_parent") else ""))
                cached_annot["cached"] = True
                return cached_annot
        except Exception as e:
            print(f"  Marker annotation cache load failed: {e}")
    
    t0 = time.time()
    ANNOT_STATUS["active"] = True
    ANNOT_STATUS["method"] = "markers"

    # Get shared DEG results (cached if already computed by Claude path)
    cluster_genes, small_clusters = _get_shared_deg(top_n=top_n, test=test)
    
    if not cluster_genes:
        ANNOT_STATUS["active"] = False
        return {"error": "DEG computation failed"}

    # Score against marker database
    n_universe = len(GENE_NAMES_LIST)
    ANNOT_STATUS["step"] = f"Scoring {len(CLUSTER_NAMES)} clusters against {len(MARKER_DB)} signatures..."
    ANNOT_STATUS["pct"] = 70
    
    results = []
    for cid, cname in enumerate(CLUSTER_NAMES):
        if cid % 20 == 0:
            pct = 70 + int(25 * cid / max(1, len(CLUSTER_NAMES)))
            ANNOT_STATUS["pct"] = pct
            ANNOT_STATUS["step"] = f"Scoring cluster {cid+1}/{len(CLUSTER_NAMES)}..."
        
        top_genes_list = cluster_genes.get(cname, [])
        top_genes = set(g.upper() for g in top_genes_list)
        n_cells_cluster = int(np.sum(CLUSTER_IDS == cid))
        
        if not top_genes:
            results.append({
                "cluster_id": cid, "cluster_name": cname,
                "n_cells": n_cells_cluster, "top_genes": [], "predictions": []
            })
            continue
        
        predictions = []
        for cell_type, markers in MARKER_DB.items():
            marker_set = set(m.upper() for m in markers)
            if len(marker_set) < 2:
                continue
            overlap = top_genes & marker_set
            if not overlap:
                continue
            score = len(overlap) / len(marker_set)
            a = len(overlap)
            b = len(marker_set) - a
            c = len(top_genes) - a
            d = n_universe - a - b - c
            _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            predictions.append({
                "cell_type": cell_type, "score": round(score, 3),
                "pval": float(f"{pval:.2e}"), "markers_found": sorted(overlap),
                "markers_total": len(marker_set), "n_found": len(overlap)
            })
        
        predictions.sort(key=lambda x: (-x["score"], x["pval"]))
        results.append({
            "cluster_id": cid, "cluster_name": cname,
            "n_cells": n_cells_cluster,
            "top_genes": sorted(top_genes)[:20],
            "predictions": predictions[:10]
        })
    
    elapsed = round(time.time() - t0, 1)
    print(f"  Marker annotation complete ({elapsed}s) — {len(results)} clusters scored")
    ANNOT_STATUS["active"] = False
    ANNOT_STATUS["step"] = ""
    ANNOT_STATUS["pct"] = 0
    
    result = {
        "results": results, "elapsed": elapsed,
        "n_clusters": len(CLUSTER_NAMES), "marker_db_size": len(MARKER_DB),
        "top_n": top_n, "cluster_names": CLUSTER_NAMES, "n_cells": N_CELLS
    }
    
    try:
        wp = _cache_write_path(".annot_markers.json")
        with open(wp, "w") as f:
            json.dump(result, f, separators=(",", ":"))
        print(f"  Marker annotation cached to {wp}")
    except Exception as e:
        print(f"  Marker annotation cache save failed: {e}")
    
    return result




# ── Obs metadata column detection ──

# Priority list for the clustering column — first match wins
_CELL_TYPE_COL_NAMES = [
    "cell_type", "cell_type_subsets", "cell_type_major", "celltype",
    "CellType", "louvain", "leiden",
]

# Columns that are already used for clustering — exclude from slice UI
_OBS_EXCLUDE = set(_CELL_TYPE_COL_NAMES) | {
    "fleek_cluster", "n_genes", "n_genes_by_counts", "total_counts",
    "total_counts_mt", "pct_counts_mt", "log1p_total_counts",
    "log1p_n_genes_by_counts", "n_counts", "n_cells",
}

def _detect_obs_cols(adata, clustering_col=None):
    """Return dict of categorical obs columns suitable for slice filtering.
    Excludes the active clustering column and known non-informative columns."""
    exclude = _OBS_EXCLUDE.copy()
    if clustering_col:
        exclude.add(clustering_col)
    result = {}
    for col in adata.obs.columns:
        if col in exclude or col.startswith("_"):
            continue
        try:
            series = adata.obs[col]
            cat = series.astype("category")
            cats = cat.cat.categories
            n_unique = len(cats)
            if n_unique < 2 or n_unique > 500:
                continue
            # Skip purely numeric columns with many unique values
            import pandas as pd
            if pd.api.types.is_float_dtype(series) and n_unique > 50:
                continue
            codes = cat.cat.codes.to_numpy().astype(np.int16)
            counts = np.bincount(codes[codes >= 0], minlength=n_unique).tolist()
            result[col] = {
                "categories": [str(c) for c in cats],
                "counts": counts,
                "codes": codes,
            }
        except Exception:
            pass
    return result

# ── LLM-based cell type annotation ──

def _load_api_key():
    """Load Anthropic API key from env var or ~/.fleek.env file.
    Env var takes priority. ~/.fleek.env should be chmod 600."""
    key = os.environ.get("ANTHROPIC_API_KEY", "")
    if key:
        return key
    env_path = Path.home() / ".fleek.env"
    if env_path.exists():
        try:
            for line in env_path.read_text().splitlines():
                line = line.strip()
                if line.startswith("ANTHROPIC_API_KEY=") and not line.startswith("#"):
                    return line.split("=", 1)[1].strip().strip('"').strip("'")
        except Exception:
            pass
    return ""

ANTHROPIC_API_KEY = _load_api_key()

# ── Shared DEG cache (used by both marker and LLM annotation) ──
_DEG_CACHE = {"key": None, "cluster_genes": None, "small_clusters": None, "full_rankings": None}
_DEG_LOCK = threading.Lock()

def _get_shared_deg(top_n=50, test="wilcoxon"):
    """Run one-vs-rest DEG and return {cluster_name: [top_genes]} dict.
    
    For large datasets (>50k cells), uses stratified subsampling — marker genes
    are virtually identical whether tested on 50k or 600k cells, but 10-20× faster.
    
    Cached in memory — second call with same dataset+test returns instantly.
    Thread-safe: first caller computes, others wait and get the cache.
    """
    import scanpy as sc
    import scipy.sparse as sp
    from collections import Counter

    DEG_MAX_CELLS = DEG_SETTINGS.get("max_cells", 50000) if DEG_SETTINGS.get("fast", True) else 0

    if ADATA is None:
        return {}, set()

    cache_key = (N_CELLS, len(CLUSTER_NAMES), test, DEG_MAX_CELLS)

    with _DEG_LOCK:
        # Cache hit — return immediately (full_rankings must also be populated if present in cache schema)
        if (_DEG_CACHE["key"] == cache_key and _DEG_CACHE["cluster_genes"] is not None
                and _DEG_CACHE["full_rankings"] is not None):
            print("  DEG results from memory cache")
            return _DEG_CACHE["cluster_genes"], _DEG_CACHE["small_clusters"]

        # ── Disk-cache hit ─────────────────────────────────────────────────
        # The full shared-DEG rankings (used by the Gene-variability cluster
        # dropdown and by both annotation paths) are expensive to recompute.
        # Persist them to .fleek_cluster_genes.json alongside the dataset so a
        # fresh server start / page reload can skip the recomputation.
        try:
            _cr, _ = _cache_read_path(".fleek_cluster_genes.json")
            if _cr:
                with open(_cr) as _f:
                    _disk = json.load(_f)
                _dk = tuple(_disk.get("key", []))
                if list(_dk) == list(cache_key):
                    # Hydrate caches from disk.
                    _DEG_CACHE["key"] = cache_key
                    _DEG_CACHE["cluster_genes"] = _disk.get("cluster_genes", {})
                    _DEG_CACHE["small_clusters"] = set(_disk.get("small_clusters", []))
                    _fr = _disk.get("full_rankings", {})
                    # full_rankings may have encoded `None` cluster entries as null
                    # in JSON — keep them as None in memory.
                    _DEG_CACHE["full_rankings"] = {k: v for k, v in _fr.items()}
                    print("  Shared DEG loaded from disk cache (.fleek_cluster_genes.json)")
                    return _DEG_CACHE["cluster_genes"], _DEG_CACHE["small_clusters"]
        except Exception as _ce:
            print(f"  Shared-DEG disk cache read failed ({_ce}); recomputing")

        t_deg = time.time()
        ANNOT_STATUS["step"] = "Preparing data..."
        ANNOT_STATUS["pct"] = 5

        adata = ADATA
        if "fleek_cluster" not in adata.obs.columns:
            adata.obs["fleek_cluster"] = [CLUSTER_NAMES[cid] for cid in CLUSTER_IDS]

        # Identify valid/small clusters BEFORE any copying
        cluster_counts = Counter(adata.obs["fleek_cluster"])
        for cname in CLUSTER_NAMES:
            if cname not in cluster_counts:
                cluster_counts[cname] = 0
        valid_clusters = {c for c, n in cluster_counts.items() if n >= 5}
        small_clusters = {c for c, n in cluster_counts.items() if n < 5}
        n_empty = sum(1 for c, n in cluster_counts.items() if n == 0)
        if small_clusters:
            print(f"  Skipping {len(small_clusters)} clusters with <5 cells ({n_empty} empty)")

        if len(valid_clusters) < 2:
            print(f"  Only {len(valid_clusters)} valid cluster(s) — need ≥2 for DEG, returning empty")
            return {cname: [] for cname in CLUSTER_NAMES}, small_clusters

        # SUBSAMPLE FIRST — before any copy or normalization.
        # On a 1.6M-cell dataset, copying the full adata for normalization
        # doubles memory usage and can trigger an OOM kill. By subsampling
        # to DEG_MAX_CELLS first (on indices only, no copy), we only ever
        # copy the small subset.
        valid_mask = adata.obs["fleek_cluster"].isin(valid_clusters).values
        valid_idx = np.where(valid_mask)[0]
        n_valid_cells = len(valid_idx)

        if n_valid_cells > DEG_MAX_CELLS:
            ANNOT_STATUS["step"] = f"Subsampling {DEG_MAX_CELLS:,} / {n_valid_cells:,} cells..."
            ANNOT_STATUS["pct"] = 8
            rng = np.random.default_rng(42)
            deg_clusters = adata.obs["fleek_cluster"].values[valid_idx]
            unique_clusters = np.unique(deg_clusters)
            cluster_local_idx = {c: np.where(deg_clusters == c)[0] for c in unique_clusters}
            total_avail = len(valid_idx)
            sub_local = []
            for c in unique_clusters:
                ci = cluster_local_idx[c]
                n_alloc = max(5, int(len(ci) / total_avail * DEG_MAX_CELLS))
                n_alloc = min(n_alloc, len(ci))
                sub_local.extend(rng.choice(ci, n_alloc, replace=False).tolist())
            sub_local = sorted(sub_local)
            use_idx = valid_idx[sub_local]
            print(f"  Subsampled {n_valid_cells:,} → {len(use_idx):,} cells for DEG (stratified)")
        else:
            use_idx = valid_idx

        # Copy only the subsample — never the full dataset
        ANNOT_STATUS["step"] = f"Preparing {len(use_idx):,} cells for DEG..."
        ANNOT_STATUS["pct"] = 10
        adata_deg = adata[use_idx].copy()

        # Route to the best available matrix (raw counts → normalize+log1p if found,
        # otherwise fall through to whatever adata.X is).
        route = MATRIX_ROUTING.get("deg_sc", {}) if MATRIX_ROUTING else {}
        src_path = route.get("path") or "X"
        if src_path != "X":
            try:
                src_full = _get_matrix_by_path(src_path)
                if src_full is not None:
                    src_sub = src_full[use_idx, :]
                    if not sp.issparse(src_sub) and hasattr(src_sub, "toarray"):
                        src_sub = src_sub.toarray()
                    adata_deg.X = src_sub
                    print(f"  Annot DEG: using {route.get('matrix_label') or src_path} as source matrix")
            except Exception as _e:
                print(f"  Annot DEG: failed to switch to {src_path} ({_e}), falling back to adata.X")

        transform = route.get("transform", "none")
        if transform == "normalize+log1p":
            ANNOT_STATUS["step"] = "Normalizing raw counts..."
            ANNOT_STATUS["pct"] = 12
            print("  Annot DEG: applying normalize_total + log1p before testing")
            sc.pp.normalize_total(adata_deg, target_sum=1e4)
            sc.pp.log1p(adata_deg)

        n_valid = len(valid_clusters)
        ANNOT_STATUS["step"] = f"Running DEG ({n_valid} clusters, {adata_deg.n_obs:,} cells)..."
        ANNOT_STATUS["pct"] = 15
        print(f"  Running one-vs-rest {test} for {n_valid} clusters on {adata_deg.n_obs:,} cells...")

        n_all_genes = adata_deg.n_vars
        try:
            sc.tl.rank_genes_groups(adata_deg, groupby="fleek_cluster", method=test,
                                    n_genes=n_all_genes, use_raw=False)
        except Exception as e:
            print(f"  rank_genes_groups failed: {e}, retrying with t-test...")
            try:
                sc.tl.rank_genes_groups(adata_deg, groupby="fleek_cluster", method="t-test",
                                        n_genes=n_all_genes, use_raw=False)
            except Exception as e2:
                print(f"  DEG failed: {e2}")
                return {}, small_clusters

        rgg = adata_deg.uns["rank_genes_groups"]

        # Extract top genes per cluster (for annotation backward compat)
        result = {}
        for cid, cname in enumerate(CLUSTER_NAMES):
            if cname in small_clusters:
                result[cname] = []
                continue
            try:
                names = rgg["names"][cname][:top_n]
                scores = rgg["scores"][cname][:top_n]
                genes = [str(g) for g, s in zip(names, scores) if s > 0]
                result[cname] = genes[:top_n]
            except (KeyError, IndexError):
                result[cname] = []

        # Build full rankings with all metrics + Cohen's d
        full_rankings = {}
        ANNOT_STATUS["step"] = "Computing effect sizes..."
        ANNOT_STATUS["pct"] = 55
        X = adata_deg.X
        cluster_labels = adata_deg.obs["fleek_cluster"].values
        for cname in CLUSTER_NAMES:
            if cname in small_clusters:
                full_rankings[cname] = None
                continue
            try:
                all_names = [str(g) for g in rgg["names"][cname]]
                all_lfc = rgg["logfoldchanges"][cname].astype(np.float64)
                all_padj = rgg["pvals_adj"][cname].astype(np.float64)
                # Cohen's d: vectorized over all genes
                mask_in = cluster_labels == cname
                n1 = int(mask_in.sum())
                n2 = int((~mask_in).sum())
                if n1 >= 2 and n2 >= 2:
                    X_in = X[mask_in]
                    X_out = X[~mask_in]
                    m1 = np.asarray(X_in.mean(axis=0)).ravel()
                    m2 = np.asarray(X_out.mean(axis=0)).ravel()
                    sq1 = np.asarray(X_in.power(2).mean(axis=0)).ravel() if sp.issparse(X_in) else np.asarray((X_in**2).mean(axis=0)).ravel()
                    sq2 = np.asarray(X_out.power(2).mean(axis=0)).ravel() if sp.issparse(X_out) else np.asarray((X_out**2).mean(axis=0)).ravel()
                    v1 = np.maximum(sq1 - m1**2, 0) * n1 / (n1 - 1)
                    v2 = np.maximum(sq2 - m2**2, 0) * n2 / (n2 - 1)
                    pooled = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
                    d_all = (m1 - m2) / np.maximum(pooled, 1e-9)
                    # Map Cohen's d by gene name (d_all is indexed by var position, not by ranking)
                    gene_to_varidx = {g: j for j, g in enumerate(adata_deg.var_names)}
                    all_d = np.array([float(d_all[gene_to_varidx[g]]) if g in gene_to_varidx else 0.0 for g in all_names])
                else:
                    all_d = np.zeros(len(all_names))
                # Replace non-finite values
                all_lfc = np.where(np.isfinite(all_lfc), all_lfc, 0.0)
                all_padj = np.where(np.isfinite(all_padj), all_padj, 1.0)
                all_d = np.where(np.isfinite(all_d), all_d, 0.0)
                full_rankings[cname] = {
                    "genes": all_names,
                    "log2fc": np.round(all_lfc, 4).tolist(),
                    "padj": all_padj.tolist(),
                    "cohens_d": np.round(all_d, 4).tolist(),
                }
            except (KeyError, IndexError):
                full_rankings[cname] = None
        print(f"  Full cluster rankings computed for {sum(1 for v in full_rankings.values() if v)} clusters")

        if adata_deg is not adata:
            del adata_deg

        elapsed_deg = time.time() - t_deg
        print(f"  Shared DEG complete ({elapsed_deg:.1f}s) — cached for reuse")

        _DEG_CACHE["key"] = cache_key
        _DEG_CACHE["cluster_genes"] = result
        _DEG_CACHE["small_clusters"] = small_clusters
        _DEG_CACHE["full_rankings"] = full_rankings

        # Persist to disk so subsequent server starts / page reloads skip
        # this (can be multi-minute) recomputation on large datasets.
        try:
            _wp = _cache_write_path(".fleek_cluster_genes.json")
            if _wp:
                _payload = {
                    "key": list(cache_key),
                    "cluster_genes": result,
                    "small_clusters": list(small_clusters),
                    "full_rankings": full_rankings,
                }
                with open(_wp, "w") as _wf:
                    json.dump(_json_safe(_payload), _wf, separators=(",", ":"))
                print(f"  Shared DEG written to {_wp}")
        except Exception as _we:
            print(f"  Shared-DEG disk cache write failed ({_we}) — keeping in-memory only")

        return result, small_clusters


def annotate_clusters_llm(top_n=30, test="wilcoxon", tissue_hint=""):
    """Use Claude API to infer cell types from top marker genes per cluster.
    Batches clusters into groups to avoid token limits.

    Returns same format as annotate_clusters for UI compatibility.
    """
    import urllib.request
    import urllib.error

    if not ANTHROPIC_API_KEY:
        return {"error": "No API key. Set ANTHROPIC_API_KEY env var or add it to ~/.fleek.env."}

    if ADATA is None:
        return {"error": "No dataset loaded"}

    # Check cache (with parent stem fallback for subsets)
    _parent_stem = ADATA.uns.get("fleek_parent_stem") if ADATA is not None else None
    cache_path, _cw = _cache_read_path(".annot_llm.json")
    if not cache_path and _parent_stem and CACHE_DIR:
        _pp = CACHE_DIR / f"{_parent_stem}.annot_llm.json"
        if _pp.exists():
            cache_path = _pp
    if cache_path:
        try:
            with open(cache_path) as f:
                cached = json.load(f)
            _cn = cached.get("cluster_names")
            _valid = (_cn == CLUSTER_NAMES and cached.get("n_cells") == N_CELLS) or (_cn and set(CLUSTER_NAMES).issubset(set(_cn)))
            if _valid:
                cached = _remap_cached_annotations(cached)
                print(f"  LLM annotation loaded from cache" + (" (remapped for subset)" if cached.get("remapped_from_parent") else ""))
                cached["cached"] = True
                return cached
        except Exception as e:
            print(f"  LLM annotation cache load failed: {e}")

    t0 = time.time()
    ANNOT_STATUS["active"] = True
    ANNOT_STATUS["method"] = "claude"
    ANNOT_STATUS["step"] = f"Computing top genes ({N_CELLS:,} cells)..."
    ANNOT_STATUS["pct"] = 5

    # Step 1: Get top genes per cluster (shared with marker path)
    print("  LLM annotation: getting shared DEG results...")
    cluster_genes, _ = _get_shared_deg(top_n=50, test=test)

    ANNOT_STATUS["step"] = "Preparing API batches..."
    ANNOT_STATUS["pct"] = 40

    # Step 2: Build cluster info list
    cluster_info = []
    for cid, cname in enumerate(CLUSTER_NAMES):
        genes = cluster_genes.get(cname, [])
        n_cells = int(np.sum(CLUSTER_IDS == cid))
        cluster_info.append({"cid": cid, "cname": cname, "n_cells": n_cells, "genes": genes})

    # Step 3: Batch clusters (~50 per call to stay well within token limits)
    BATCH_SIZE = 50
    all_predictions_raw = []
    total_in_tok = 0
    total_out_tok = 0
    model_name = "claude"
    n_batches = (len(cluster_info) + BATCH_SIZE - 1) // BATCH_SIZE

    tissue_ctx = f" from {tissue_hint} tissue" if tissue_hint else ""

    for batch_idx in range(n_batches):
        batch = cluster_info[batch_idx * BATCH_SIZE : (batch_idx + 1) * BATCH_SIZE]
        cluster_lines = []
        for ci in batch:
            if ci["genes"]:
                cluster_lines.append(f'  Cluster "{ci["cname"]}" ({ci["n_cells"]:,} cells): {", ".join(ci["genes"][:25])}')
            else:
                cluster_lines.append(f'  Cluster "{ci["cname"]}" ({ci["n_cells"]:,} cells): [no significant markers]')

        prompt = f"""You are an expert single-cell RNA-seq bioinformatician. Below are clusters from a scRNA-seq dataset{tissue_ctx} with their top differentially expressed marker genes (one-vs-rest, ranked by significance).

For each cluster, identify the most likely cell type. Return ONLY valid JSON — no markdown, no explanation outside the JSON. Use this exact format:

[
  {{"cluster": "cluster_name", "cell_type": "predicted type", "confidence": "high|medium|low", "reasoning": "brief explanation citing key markers"}},
  ...
]

Clusters and their top marker genes:
{chr(10).join(cluster_lines)}

Respond with ONLY the JSON array. No markdown fences. No other text."""

        pct = 40 + int(50 * batch_idx / max(1, n_batches))
        ANNOT_STATUS["step"] = f"Claude API batch {batch_idx+1}/{n_batches}..."
        ANNOT_STATUS["pct"] = pct
        print(f"  LLM annotation: batch {batch_idx+1}/{n_batches} ({len(batch)} clusters)...")
        request_body = json.dumps({
            "model": "claude-sonnet-4-20250514",
            "max_tokens": 8192,
            "messages": [{"role": "user", "content": prompt}]
        }).encode("utf-8")

        req = urllib.request.Request(
            "https://api.anthropic.com/v1/messages",
            data=request_body,
            headers={
                "Content-Type": "application/json",
                "x-api-key": ANTHROPIC_API_KEY,
                "anthropic-version": "2023-06-01"
            },
            method="POST"
        )

        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                resp_data = json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", errors="replace")
            return {"error": f"Claude API error {e.code}: {body[:200]}"}
        except Exception as e:
            return {"error": f"Claude API request failed: {str(e)}"}

        model_name = resp_data.get("model", model_name)
        usage = resp_data.get("usage", {})
        total_in_tok += usage.get("input_tokens", 0)
        total_out_tok += usage.get("output_tokens", 0)

        # Parse response
        text = ""
        for block in resp_data.get("content", []):
            if block.get("type") == "text":
                text += block["text"]

        text = text.strip()
        if text.startswith("```"):
            text = text.split("\n", 1)[1] if "\n" in text else text[3:]
        if text.endswith("```"):
            text = text[:-3].strip()
        if text.startswith("json"):
            text = text[4:].strip()

        try:
            batch_predictions = json.loads(text)
            all_predictions_raw.extend(batch_predictions)
        except json.JSONDecodeError as e:
            print(f"  Warning: batch {batch_idx+1} JSON parse failed: {e}")
            # Try to salvage partial results
            continue

    # Step 4: Convert to standard format
    results = []
    for ci in cluster_info:
        cid = ci["cid"]
        cname = ci["cname"]

        pred = None
        for p in all_predictions_raw:
            if p.get("cluster") == cname:
                pred = p
                break

        predictions = []
        if pred:
            conf = pred.get("confidence", "medium")
            score = {"high": 0.9, "medium": 0.6, "low": 0.3}.get(conf, 0.5)
            predictions.append({
                "cell_type": pred.get("cell_type", "Unknown"),
                "score": score,
                "confidence": conf,
                "reasoning": pred.get("reasoning", ""),
                "n_found": 0,
                "markers_total": len(ci["genes"]),
                "pval": 0
            })

        results.append({
            "cluster_id": cid,
            "cluster_name": cname,
            "n_cells": ci["n_cells"],
            "top_genes": ci["genes"][:20],
            "predictions": predictions
        })

    elapsed = round(time.time() - t0, 1)
    print(f"  LLM annotation complete ({elapsed}s, {model_name}, {total_in_tok}+{total_out_tok} tokens, {n_batches} batches)")
    ANNOT_STATUS["active"] = False
    ANNOT_STATUS["step"] = ""
    ANNOT_STATUS["pct"] = 0

    result = {
        "results": results,
        "elapsed": elapsed,
        "n_clusters": len(CLUSTER_NAMES),
        "method": "llm",
        "model": model_name,
        "tokens": {"input": total_in_tok, "output": total_out_tok},
        "top_n": top_n,
        "cluster_names": CLUSTER_NAMES,
        "n_cells": N_CELLS
    }

    # Save cache (always — cache_path may be None if file didn't exist yet)
    try:
        wp = _cache_write_path(".annot_llm.json")
        with open(wp, "w") as f:
            json.dump(result, f, separators=(",", ":"))
        print(f"  LLM annotation cached to {wp}")
    except Exception as e:
        print(f"  LLM annotation cache save failed: {e}")

    return result

def _progress(msg, pct=None):
    """Update global progress for client polling. Raises AbortError if abort was requested."""
    PROGRESS["message"] = msg
    if pct is not None:
        PROGRESS["pct"] = pct
    PROGRESS["status"] = "processing"
    print(f"  [{pct or '?'}%] {msg}")
    _check_abort()


def _get_available_ram_gb():
    """Estimate available RAM in GB."""
    try:
        import psutil
        return psutil.virtual_memory().available / 1e9
    except ImportError:
        # Assume 16GB if psutil not available
        return 16.0


class AbortError(Exception):
    pass

def _check_abort():
    """Raise AbortError if abort was requested."""
    if ABORT_REQUESTED:
        raise AbortError("Loading aborted by user")

def _reset_all():
    """Free all data globals and reclaim memory.

    On machines near their RAM limit this matters — loading dataset B while
    dataset A is still held would briefly peak at A+B. We drop all top-level
    references, close any backed h5 file handle, and run gc.collect() twice
    (some scipy/numpy buffers are only reclaimed on the second pass because
    the first releases weakref-targets that the second then finalizes).
    """
    global ADATA, UMAP_2D, UMAP_3D, PCA_2D, PCA_3D, PACMAP_2D, PACMAP_3D
    global PACMAP_COMPUTING, CLUSTER_IDS, CLUSTER_NAMES, GENE_NAMES_LIST
    global N_CELLS, BACKED, X_CSC, HVG_NAMES, HVG_VAR, ALL_GENE_VAR, ALL_GENE_VAR_METRIC, GENE_CUTOFF_FLAG, GENE_INDEX, CSC_BUILDING, OBS_COLS
    global CSC_CACHED, CSC_TIME, LOADED_PATH, LOAD_SETTINGS, ABORT_REQUESTED
    global COUNTS_MATRIX, COUNTS_LABEL, COUNTS_INFO
    # Close backed-mode file handle first — anndata keeps an h5py.File open
    # when adata is loaded with backed=True, and the file (plus its page
    # cache) won't be freed until the handle is closed explicitly.
    if ADATA is not None:
        try:
            if getattr(ADATA, "isbacked", False) and getattr(ADATA, "file", None) is not None:
                ADATA.file.close()
        except Exception as _e:
            print(f"  _reset_all: failed to close backed file ({_e})")
    ADATA = None
    UMAP_2D = None; UMAP_3D = None
    PCA_2D = None; PCA_3D = None
    PACMAP_2D = None; PACMAP_3D = None
    global EMBEDDING_SOURCES
    EMBEDDING_SOURCES = {}
    PACMAP_COMPUTING = False
    CLUSTER_IDS = None; CLUSTER_NAMES = None; CLUSTER_COL = None; CLUSTER_COLORS = []
    OBS_COLS.clear()
    GENE_NAMES_LIST = None; N_CELLS = 0
    BACKED = False; X_CSC = None
    HVG_NAMES = []; HVG_VAR = {}; ALL_GENE_VAR = None; GENE_CUTOFF_FLAG = None; GENE_INDEX = {}
    CSC_BUILDING = False; CSC_CACHED = False; CSC_TIME = 0
    LOADED_PATH = ""; LOAD_SETTINGS.clear()
    COUNTS_MATRIX = None; COUNTS_LABEL = None; COUNTS_INFO = None
    global MATRIX_REGISTRY, MATRIX_ROUTING, _CELL_LIB_SIZES
    MATRIX_REGISTRY = []; MATRIX_ROUTING = {}; _CELL_LIB_SIZES = None
    ABORT_REQUESTED = False
    # Clear DEG cache (reinitialize with expected keys, not .clear())
    _DEG_CACHE["key"] = None
    _DEG_CACHE["cluster_genes"] = None
    _DEG_CACHE["small_clusters"] = None
    _DEG_CACHE["full_rankings"] = None
    import gc
    gc.collect()
    gc.collect()

def load_and_prepare(h5ad_path, max_cells=0, n_dims_list=[2, 3], fast_umap=False,
                     fast_umap_subsample=50000, backed="off", native_2d=False):
    """Load h5ad, subsample if needed, compute UMAP(s).
    
    backed: "auto" (use file size vs RAM), "on" (force backed), "off" (force in-memory)
    """
    global ADATA, UMAP_2D, UMAP_3D, PCA_2D, PCA_3D, PACMAP_2D, PACMAP_3D, PACMAP_COMPUTING, CLUSTER_IDS, CLUSTER_NAMES, CLUSTER_COL, CLUSTER_COLORS, GENE_NAMES_LIST, N_CELLS, BACKED, X_CSC, HVG_NAMES, HVG_VAR, ALL_GENE_VAR, ALL_GENE_VAR_METRIC, GENE_CUTOFF_FLAG, GENE_INDEX, CSC_BUILDING, CSC_CACHED, CSC_TIME, LOADED_PATH, ABORT_REQUESTED

    import scanpy as sc
    import scipy.sparse as sp
    import gc

    ABORT_REQUESTED = False  # Reset at start of new load

    # Free old dataset before loading new one
    if ADATA is not None:
        _progress("Freeing previous dataset...", 2)
        _reset_all()

    _check_abort()
    _progress("Loading h5ad file...", 5)
    LOADED_PATH = str(h5ad_path)
    t0 = time.time()

    # Decide backed mode
    file_size_gb = os.path.getsize(h5ad_path) / 1e9
    avail_ram_gb = _get_available_ram_gb()
    use_backed = False
    if backed == "on":
        use_backed = True
    elif backed == "auto":
        # Use backed if file is > 40% of available RAM
        use_backed = file_size_gb > avail_ram_gb * 0.4
    # backed == "off" → use_backed stays False

    if use_backed:
        _progress(f"Loading in backed mode ({file_size_gb:.1f}GB file, {avail_ram_gb:.0f}GB RAM avail)...", 5)
        adata = sc.read_h5ad(h5ad_path, backed="r")
        BACKED = True
    else:
        adata = sc.read_h5ad(h5ad_path)
        BACKED = False
    _progress(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes"
              f" ({'backed' if BACKED else 'in-memory'}, {time.time()-t0:.0f}s)", 10)

    # Subsample
    if max_cells > 0 and adata.n_obs > max_cells:
        _progress(f"Subsampling to {max_cells:,} cells...", 12)
        rng = np.random.default_rng(42)
        idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
        idx.sort()
        if BACKED:
            # Read subset into memory — it's small by definition
            _progress(f"Reading {max_cells:,}-cell subset from disk...", 13)
            adata = adata[idx].to_memory()
            BACKED = False
        else:
            adata = adata[idx].copy()

    N_CELLS = adata.n_obs
    _progress(f"Working with {N_CELLS:,} cells", 15)

    # Load fleek cache (UMAP, PCA, leiden, PaCMAP — all in one file)
    _init_cache_dir()
    _cr, _cr_writable = _cache_read_path(".fleek_cache.npz")
    cache_read = _cr  # may be None
    old_cache = _dataset_cache_path(".umap_cache.npz")
    # Migrate old cache if needed
    if not cache_read and old_cache and old_cache.exists():
        target = _cache_write_path(".fleek_cache.npz")
        try:
            old_cache.rename(target)
            print(f"  Migrated old cache {old_cache.name} → {target.name}")
            cache_read = target
            _cr_writable = True
        except (OSError, PermissionError):
            import shutil
            shutil.copy2(str(old_cache), str(target))
            cache_read = target
            _cr_writable = True
    cached = {}
    _cache_time = 0.0  # accumulate time spent loading caches
    _CACHE_SKIP_VALIDATE = {"leiden_names"}  # byte arrays, not cell-indexed
    if cache_read and cache_read.exists():
        _progress("Checking cache...", 20)
        _ct0 = time.time()
        try:
            cached = dict(np.load(cache_read, allow_pickle=False))
            for k in list(cached.keys()):
                if k in _CACHE_SKIP_VALIDATE:
                    continue
                if cached[k].shape[0] != N_CELLS:
                    _progress("Cache mismatch, will recompute", 20)
                    print(f"  Cache key '{k}' has {cached[k].shape[0]} rows, expected {N_CELLS}")
                    cached = {}
                    break
        except Exception as e:
            _progress(f"Cache load failed, recomputing", 20)
            cached = {}
        _cache_time += time.time() - _ct0

    need_save = False
    neighbors_built = False

    # Cluster IDs — try priority list of known column names, then cache, then compute
    global OBS_COLS
    _cluster_col_used = None
    for _cname in _CELL_TYPE_COL_NAMES:
        if _cname in adata.obs.columns:
            cats = adata.obs[_cname].astype("category")
            CLUSTER_IDS = cats.cat.codes.to_numpy().astype(np.int32)
            CLUSTER_NAMES = cats.cat.categories.tolist()
            if _cname in ("leiden", "louvain"):
                CLUSTER_NAMES = [f"Cluster {c}" for c in cats.cat.categories]
            _cluster_col_used = _cname
            print(f"  Cell types from obs column: '{_cname}' ({len(CLUSTER_NAMES)} types)")
            break
    if _cluster_col_used is None and "leiden_ids" in cached and "leiden_names" in cached:
        # Load cached leiden results
        CLUSTER_IDS = cached["leiden_ids"].astype(np.int32)
        import json as _json
        CLUSTER_NAMES = _json.loads(cached["leiden_names"].tobytes().decode("utf-8"))
        print(f"  Leiden clusters loaded from cache ({len(CLUSTER_NAMES)} clusters)")
    elif _cluster_col_used is None:
        # Run leiden clustering — only when no obs column or cached result was found
        print("  No cell_type or leiden found, computing clusters...")
        _progress("Computing clusters (PCA + neighbors + leiden)...", 15)
        if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
            _ensure_neighbors(adata)
        _progress("Running leiden clustering...", 18)
        sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
        cats = adata.obs["leiden"].astype("category")
        CLUSTER_IDS = cats.cat.codes.to_numpy().astype(np.int32)
        CLUSTER_NAMES = [f"Cluster {c}" for c in cats.cat.categories]
        _cluster_col_used = "leiden"
        # Cache leiden results
        import json as _json
        cached["leiden_ids"] = CLUSTER_IDS
        cached["leiden_names"] = np.frombuffer(_json.dumps(CLUSTER_NAMES).encode("utf-8"), dtype=np.uint8)
        need_save = True

    # ── Cluster colors ──────────────────────────────────────────────────────
    # Priority: fleek_color_map (name→hex) > uns["{col}_colors"] > palette fallback.
    # fleek_color_map is stored as parallel {"names": [...], "colors": [...]}
    # arrays (newer format, slash-safe for HDF5) OR as a plain {name:color} dict
    # (legacy format — still read for backward compatibility, but keys containing
    # "/" will be missing from this format because HDF5 strips them on write).
    CLUSTER_COL = _cluster_col_used
    _color_map_obj = adata.uns.get("fleek_color_map", None)
    _name_to_color = {}
    if isinstance(_color_map_obj, dict) and _color_map_obj:
        _nm = _color_map_obj.get("names")
        _cs = _color_map_obj.get("colors")
        if _nm is not None and _cs is not None:
            # New parallel-array format
            try:
                _nm_list = list(_nm) if not isinstance(_nm, list) else _nm
                _cs_list = list(_cs) if not isinstance(_cs, list) else _cs
                for _k, _v in zip(_nm_list, _cs_list):
                    _name_to_color[str(_k)] = str(_v)
            except Exception:
                pass
        else:
            # Legacy plain-dict format
            for _k, _v in _color_map_obj.items():
                _name_to_color[str(_k)] = str(_v)
    # Also grab the positional {col}_colors array as a secondary fallback for
    # names that fleek_color_map doesn't cover (e.g. legacy subsets where HDF5
    # dropped "/"-containing dict keys on write).
    _uns_colors = []
    if _cluster_col_used:
        _uns_colors = adata.uns.get(f"{_cluster_col_used}_colors", [])
    if isinstance(_uns_colors, np.ndarray):
        _uns_colors = _uns_colors.tolist()
    elif not isinstance(_uns_colors, list):
        _uns_colors = list(_uns_colors) if _uns_colors else []

    def _resolve_color(i, n):
        if n in _name_to_color:
            return str(_name_to_color[n])
        if i < len(_uns_colors) and _uns_colors[i]:
            return str(_uns_colors[i])
        return _palette_hex(i)
    CLUSTER_COLORS = [_resolve_color(i, n) for i, n in enumerate(CLUSTER_NAMES)]

    # Detect slice-able obs metadata columns
    OBS_COLS = _detect_obs_cols(adata, clustering_col=_cluster_col_used)
    if OBS_COLS:
        print(f"  Obs slice columns: {', '.join(OBS_COLS.keys())}")

    # Gene names (use var index, or feature_name if available)
    if "feature_name" in adata.var.columns:
        adata.var.index = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()
        GENE_NAMES_LIST = adata.var_names.tolist()
    else:
        adata.var.index = adata.var.index.astype(str)
        adata.var_names_make_unique()
        GENE_NAMES_LIST = adata.var_names.tolist()

    # Compute UMAP(s) — with caching to avoid recomputation

    # Decide if we should use subsample-and-project
    use_fast = fast_umap and N_CELLS > fast_umap_subsample

    # Figure out which dims need computing
    dims_needed = []
    umap_suffix = "_quick" if use_fast else "_full"
    umap_from_cache = True
    # Pre-existing obsm UMAP (e.g. exported subset or pre-processed h5ad)
    _obsm_umap = adata.obsm.get("X_umap")
    if _obsm_umap is not None and _obsm_umap.shape[0] != N_CELLS:
        _obsm_umap = None  # stale / wrong size — ignore
    for nd in n_dims_list:
        # Check method-specific key first, then legacy key, then obsm
        cache_key = f"umap_{nd}d{umap_suffix}"
        legacy_key = f"umap_{nd}d"
        if cache_key in cached:
            _progress(f"Loaded {nd}D UMAP from cache ({'quick' if use_fast else 'full'})", 25)
            if nd == 2:
                UMAP_2D = cached[cache_key].astype(np.float32)
            else:
                UMAP_3D = cached[cache_key].astype(np.float32)
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "cache_quick" if use_fast else "cache_full"
        elif legacy_key in cached:
            _progress(f"Loaded {nd}D UMAP from cache (legacy)", 25)
            if nd == 2:
                UMAP_2D = cached[legacy_key].astype(np.float32)
            else:
                UMAP_3D = cached[legacy_key].astype(np.float32)
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "cache_legacy"
        elif _obsm_umap is not None and _obsm_umap.shape[1] >= nd:
            # Use pre-existing UMAP from obsm (e.g. subset exported from full dataset)
            umap_arr = _obsm_umap[:, :nd].astype(np.float32)
            if nd == 2:
                UMAP_2D = umap_arr
            else:
                UMAP_3D = umap_arr
            cached[f"umap_{nd}d_full"] = umap_arr  # cache so next load is instant
            need_save = True
            _progress(f"Loaded {nd}D UMAP from obsm (pre-existing)", 25)
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "obsm"
        else:
            dims_needed.append(nd)
            umap_from_cache = False

    # If not native_2d and we need both dims, compute only 3D and derive 2D
    if not native_2d and 2 in dims_needed and 3 in dims_needed:
        dims_needed = [3]  # compute only 3D; 2D will be derived after

    if dims_needed and use_fast:
        # ── Shared setup: subsample + PCA (done ONCE for all dims) ──
        import umap as umap_lib
        import scipy.sparse as sp

        _progress(f"QuickMap: subsampling {fast_umap_subsample:,} cells...", 28)
        rng = np.random.default_rng(42)
        sub_idx = _stratified_subsample(CLUSTER_IDS, fast_umap_subsample, rng)

        sub_adata = adata[sub_idx]
        if BACKED:
            try:
                sub_adata = sub_adata.to_memory()
            except AttributeError:
                sub_adata = sub_adata.copy()
        else:
            sub_adata = sub_adata.copy()

        _progress(f"QuickMap: PCA on {len(sub_idx):,}-cell subsample...", 30)
        sc.pp.normalize_total(sub_adata, target_sum=1e4)
        sc.pp.log1p(sub_adata)
        sc.pp.highly_variable_genes(sub_adata, n_top_genes=2000)
        sub_adata = sub_adata[:, sub_adata.var.highly_variable].copy()

        # Capture scaling params for batch projection
        X_pre = sub_adata.X.toarray() if sp.issparse(sub_adata.X) else np.asarray(sub_adata.X)
        gene_mean = X_pre.mean(axis=0).astype(np.float64)
        gene_std = X_pre.std(axis=0).astype(np.float64)
        gene_std[gene_std == 0] = 1.0
        hvg_names = sub_adata.var_names.tolist()
        del X_pre

        sc.pp.scale(sub_adata, max_value=10)
        sc.tl.pca(sub_adata, n_comps=30)
        pca_comp = sub_adata.varm["PCs"][:, :30].copy()
        X_sub = sub_adata.obsm["X_pca"][:, :30].copy()
        del sub_adata  # free subsample

        # ── Batch PCA projection for all cells (done ONCE) ──
        if "X_pca" not in adata.obsm:
            _progress(f"QuickMap: projecting {N_CELLS:,} cells to PCA space...", 35)
            hvg_set = set(hvg_names)
            hvg_idx = [i for i, g in enumerate(adata.var_names) if g in hvg_set]

            batch_size = 20000
            all_pca = np.empty((N_CELLS, 30), dtype=np.float32)
            is_sparse = sp.issparse(adata.X)
            for bi in range(0, N_CELLS, batch_size):
                be = min(bi + batch_size, N_CELLS)
                chunk_full = adata.X[bi:be, :]
                if is_sparse:
                    row_sums = np.asarray(chunk_full.sum(axis=1)).ravel()
                else:
                    row_sums = np.asarray(chunk_full).sum(axis=1)
                row_sums[row_sums == 0] = 1
                scale_factors = (1e4 / row_sums).reshape(-1, 1)
                chunk_hvg = chunk_full[:, hvg_idx]
                if is_sparse:
                    chunk_hvg = chunk_hvg.toarray()
                chunk_hvg = np.asarray(chunk_hvg, dtype=np.float64)
                chunk_hvg *= scale_factors
                np.log1p(chunk_hvg, out=chunk_hvg)
                chunk_hvg = (chunk_hvg - gene_mean) / gene_std
                np.clip(chunk_hvg, -10, 10, out=chunk_hvg)
                all_pca[bi:be] = (chunk_hvg @ pca_comp).astype(np.float32)
                if N_CELLS > batch_size:
                    pct = 35 + int((be / N_CELLS) * 10)
                    _progress(f"QuickMap: PCA projected {be:,}/{N_CELLS:,} cells", pct)

            adata.obsm["X_pca"] = all_pca
            del all_pca, pca_comp, gene_mean, gene_std

        X_all = adata.obsm["X_pca"][:, :30]
        sub_set = set(sub_idx.tolist())
        rest_idx = np.array([i for i in range(N_CELLS) if i not in sub_set])

        # ── Per-dimension: fit UMAP + project (only the UMAP differs) ──
        for di, nd in enumerate(dims_needed):
            base_pct = 48 + di * 25
            _progress(f"QuickMap {nd}D: fitting UMAP on {len(sub_idx):,} cells...", base_pct)
            t0 = time.time()

            reducer = umap_lib.UMAP(
                n_components=nd, n_neighbors=15, min_dist=0.3,
                low_memory=True, n_jobs=-1
            )
            emb_sub = reducer.fit_transform(X_sub)
            _progress(f"QuickMap {nd}D: projecting {len(rest_idx):,} remaining cells...", base_pct + 12)

            full_emb = np.empty((N_CELLS, nd), dtype=np.float32)
            full_emb[sub_idx] = emb_sub.astype(np.float32)

            if len(rest_idx) > 0:
                batch_size = 50000
                n_rest = len(rest_idx)
                for bi in range(0, n_rest, batch_size):
                    batch_end = min(bi + batch_size, n_rest)
                    batch_idx = rest_idx[bi:batch_end]
                    emb_batch = reducer.transform(X_all[batch_idx])
                    full_emb[batch_idx] = emb_batch.astype(np.float32)
                    proj_pct = base_pct + 12 + int(batch_end / n_rest * 12)
                    _progress(f"QuickMap {nd}D: projected {batch_end:,}/{n_rest:,} cells", proj_pct)

            elapsed = time.time() - t0
            _progress(f"QuickMap {nd}D complete ({elapsed:.0f}s)", base_pct + 24)

            if nd == 2:
                UMAP_2D = full_emb
            else:
                UMAP_3D = full_emb
            cached[f"umap_{nd}d_quick"] = full_emb
            need_save = True
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "computed_quick"

    elif dims_needed:
        # ── Standard full UMAP ──
        neighbors_built = False
        for di, nd in enumerate(dims_needed):
            base_pct = 30 + di * 30
            if not neighbors_built:
                _progress(f"Building neighbor graph for {N_CELLS:,} cells...", base_pct)
                _ensure_neighbors(adata)
                neighbors_built = True

            _progress(f"Computing {nd}D UMAP for {N_CELLS:,} cells...", base_pct + 10)
            t0 = time.time()
            sc.tl.umap(adata, n_components=nd)
            _progress(f"{nd}D UMAP done ({time.time()-t0:.0f}s)", base_pct + 28)

            umap_arr = adata.obsm["X_umap"].astype(np.float32)
            if nd == 2:
                UMAP_2D = umap_arr
            else:
                UMAP_3D = umap_arr
            cached[f"umap_{nd}d_full"] = umap_arr
            need_save = True
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "computed_full"

    # Derive 2D from 3D if not native and 2D wasn't loaded/computed
    if UMAP_2D is None and UMAP_3D is not None:
        UMAP_2D = UMAP_3D[:, :2].copy()
        print(f"  UMAP 2D derived from 3D (slice)")
        EMBEDDING_SOURCES["umap_2d"] = "derived_from_3d"

    if need_save:
        # Also save PCA coords to cache if available
        if "X_pca" in adata.obsm:
            pca_arr = adata.obsm["X_pca"].astype(np.float32)
            cached["pca_2d"] = pca_arr[:, :2].copy()
            cached["pca_3d"] = pca_arr[:, :3].copy()
            cached["pca_full"] = pca_arr.copy()
        _progress("Saving UMAP cache...", 95)
        np.savez(_cache_write_path(".fleek_cache.npz"), **cached)

    ADATA = adata

    # ── Extract PCA coordinates ──
    # Priority: obsm['X_pca'] in the h5ad → cache (fallback for FLEEK-computed
    # cases where the file had no PCA to begin with). When obsm has a PCA that
    # matches the current cell count we always trust it over the cache, because
    # the file is the authoritative source — a cached PCA could be stale if the
    # dataset was regenerated with a different upstream pipeline.
    _obsm_pca = adata.obsm.get("X_pca") if hasattr(adata, "obsm") else None
    if _obsm_pca is not None and getattr(_obsm_pca, "shape", (0,))[0] != N_CELLS:
        _obsm_pca = None  # stale/wrong size — ignore

    if _obsm_pca is not None:
        pca_arr = np.asarray(_obsm_pca, dtype=np.float32)
        PCA_2D = pca_arr[:, :2].copy()
        PCA_3D = pca_arr[:, :3].copy() if pca_arr.shape[1] >= 3 else np.pad(PCA_2D, ((0, 0), (0, 1)))
        print(f"  PCA from obsm ({pca_arr.shape[1]} components) — preferred over cache")
        EMBEDDING_SOURCES["pca_2d"] = "obsm"
        EMBEDDING_SOURCES["pca_3d"] = "obsm"
    elif "pca_full" in cached:
        pca_full = cached["pca_full"].astype(np.float32)
        PCA_2D = pca_full[:, :2].copy()
        PCA_3D = pca_full[:, :3].copy()
        # Restore to obsm so PaCMAP and other tools can use it
        adata.obsm["X_pca"] = pca_full
        print(f"  PCA loaded from cache ({pca_full.shape[1]} components) — obsm had none")
        EMBEDDING_SOURCES["pca_2d"] = "cache_full"
        EMBEDDING_SOURCES["pca_3d"] = "cache_full"
    elif "pca_2d" in cached and "pca_3d" in cached:
        # Legacy cache with only 2D/3D — use for display but PaCMAP won't work
        PCA_2D = cached["pca_2d"].astype(np.float32)
        PCA_3D = cached["pca_3d"].astype(np.float32)
        print(f"  PCA loaded from cache (2D/3D only, no full PCA for PaCMAP) — obsm had none")
        EMBEDDING_SOURCES["pca_2d"] = "cache_legacy"
        EMBEDDING_SOURCES["pca_3d"] = "cache_legacy"
    else:
        PCA_2D = None
        PCA_3D = None
        print(f"  PCA not available (will compute in background)")
        # Compute in background to avoid blocking load
        def _bg_pca():
            global PCA_2D, PCA_3D
            try:
                import scanpy as sc
                sc.tl.pca(adata, n_comps=min(50, adata.n_vars - 1))
                pca_arr = adata.obsm["X_pca"].astype(np.float32)
                PCA_2D = pca_arr[:, :2].copy()
                PCA_3D = pca_arr[:, :3].copy()
                cached["pca_2d"] = PCA_2D
                cached["pca_3d"] = PCA_3D
                cached["pca_full"] = pca_arr.copy()
                try:
                    np.savez(_cache_write_path(".fleek_cache.npz"), **cached)
                except Exception:
                    pass
                print(f"  PCA computed in background ({pca_arr.shape[1]} components)")
            except Exception as e:
                print(f"  Background PCA failed: {e}")
        threading.Thread(target=_bg_pca, daemon=True).start()

    # ── PaCMAP embedding (background subprocess to avoid numba/threading issues) ──
    _pm_suffix = "_quick" if PACMAP_SETTINGS.get("fast", True) else "_full"
    _pm_key2 = f"pacmap_2d{_pm_suffix}"
    _pm_key3 = f"pacmap_3d{_pm_suffix}"
    if _pm_key2 in cached and _pm_key3 in cached:
        PACMAP_2D = cached[_pm_key2].astype(np.float32)
        PACMAP_3D = cached[_pm_key3].astype(np.float32)
        PACMAP_COMPUTING = False
        print(f"  PaCMAP loaded from cache ({'quick' if _pm_suffix == '_quick' else 'full'})")
        src = "cache_quick" if _pm_suffix == "_quick" else "cache_full"
        EMBEDDING_SOURCES["pacmap_2d"] = src
        EMBEDDING_SOURCES["pacmap_3d"] = src
    elif _pm_key3 in cached:
        # 3D cached but no native 2D — derive 2D from 3D
        PACMAP_3D = cached[_pm_key3].astype(np.float32)
        PACMAP_2D = PACMAP_3D[:, :2].copy()
        PACMAP_COMPUTING = False
        print(f"  PaCMAP 3D from cache, 2D derived (slice)")
        EMBEDDING_SOURCES["pacmap_3d"] = "cache_full" if _pm_suffix == "_full" else "cache_quick"
        EMBEDDING_SOURCES["pacmap_2d"] = "derived_from_3d"
    elif "pacmap_2d" in cached and "pacmap_3d" in cached:
        # Legacy cache fallback
        PACMAP_2D = cached["pacmap_2d"].astype(np.float32)
        PACMAP_3D = cached["pacmap_3d"].astype(np.float32)
        PACMAP_COMPUTING = False
        print(f"  PaCMAP loaded from cache (legacy)")
        EMBEDDING_SOURCES["pacmap_2d"] = "cache_legacy"
        EMBEDDING_SOURCES["pacmap_3d"] = "cache_legacy"
    elif "pacmap_3d" in cached:
        PACMAP_3D = cached["pacmap_3d"].astype(np.float32)
        PACMAP_2D = PACMAP_3D[:, :2].copy()
        PACMAP_COMPUTING = False
        print(f"  PaCMAP 3D from cache (legacy), 2D derived (slice)")
        EMBEDDING_SOURCES["pacmap_3d"] = "cache_legacy"
        EMBEDDING_SOURCES["pacmap_2d"] = "derived_from_3d"
    elif "X_pacmap" in adata.obsm and adata.obsm["X_pacmap"].shape[0] == N_CELLS:
        # Pre-existing PaCMAP from obsm (e.g. exported subset)
        _obsm_pm = adata.obsm["X_pacmap"].astype(np.float32)
        if _obsm_pm.shape[1] >= 2:
            PACMAP_2D = _obsm_pm[:, :2]
        if _obsm_pm.shape[1] >= 3:
            PACMAP_3D = _obsm_pm[:, :3]
        elif PACMAP_2D is not None:
            PACMAP_3D = np.pad(PACMAP_2D, ((0, 0), (0, 1)))  # pad with zeros if only 2D available
        PACMAP_COMPUTING = False
        if PACMAP_2D is not None:
            cached["pacmap_2d_full"] = PACMAP_2D
            EMBEDDING_SOURCES["pacmap_2d"] = "obsm"
        if PACMAP_3D is not None:
            cached["pacmap_3d_full"] = PACMAP_3D
            if EMBEDDING_SOURCES.get("pacmap_3d") is None:
                EMBEDDING_SOURCES["pacmap_3d"] = "obsm" if _obsm_pm.shape[1] >= 3 else "padded_from_2d"
        need_save = True
        print(f"  PaCMAP loaded from obsm (pre-existing)")
    elif "X_pca" in adata.obsm:
        def _bg_pacmap_subprocess():
            global PACMAP_2D, PACMAP_3D, PACMAP_COMPUTING
            import tempfile, subprocess
            PACMAP_COMPUTING = True
            try:
                # Check if pacmap is installed
                proc_check = subprocess.run(
                    [sys.executable, "-c", "import pacmap"],
                    capture_output=True, timeout=10
                )
                if proc_check.returncode != 0:
                    print("  PaCMAP not installed (pip install pacmap)")
                    PACMAP_COMPUTING = False
                    return

                pca_input = adata.obsm["X_pca"][:, :30].astype(np.float32)
                n_total = pca_input.shape[0]
                PACMAP_SUBSAMPLE = PACMAP_SETTINGS.get("max_cells", 50000) if PACMAP_SETTINGS.get("fast", True) else 0

                with tempfile.TemporaryDirectory() as tmpdir:
                    input_path = os.path.join(tmpdir, "pca_in.npy")
                    out2_path = os.path.join(tmpdir, "pacmap_2d.npy")
                    out3_path = os.path.join(tmpdir, "pacmap_3d.npy")
                    np.save(input_path, pca_input)

                    if PACMAP_SUBSAMPLE > 0 and n_total > PACMAP_SUBSAMPLE:
                        # Stratified subsample + nearest-neighbor projection
                        sub_idx = _stratified_subsample(CLUSTER_IDS, PACMAP_SUBSAMPLE,
                                                        np.random.default_rng(42))
                        sub_path = os.path.join(tmpdir, "sub_idx.npy")
                        np.save(sub_path, np.array(sub_idx, dtype=np.int32))

                        _do_native_pm = native_2d
                        script = f'''
import numpy as np, pacmap, time
from scipy.spatial import cKDTree
X = np.load("{input_path}")
sub_idx = np.load("{sub_path}")
X_sub = X[sub_idx]
nn = min(10, max(2, len(sub_idx) // 50))
t0 = time.time()
n_total = X.shape[0]
n_sub = len(sub_idx)
do_native_2d = {_do_native_pm}
if do_native_2d:
    print(f"  PaCMAP 2D ({{n_sub:,}} / {{n_total:,}} cells, {{nn}} neighbors)...")
    pm2 = pacmap.PaCMAP(n_components=2, n_neighbors=nn, verbose=False)
    Y2_sub = pm2.fit_transform(X_sub).astype(np.float32)
print(f"  PaCMAP 3D{{'' if do_native_2d else ' (2D will be derived)'}}")
pm3 = pacmap.PaCMAP(n_components=3, n_neighbors=nn, verbose=False)
Y3_sub = pm3.fit_transform(X_sub).astype(np.float32)
if not do_native_2d:
    Y2_sub = Y3_sub[:, :2].copy()
# Project remaining cells via k-NN weighted interpolation in PCA space
print(f"  Projecting {{n_total - n_sub:,}} remaining cells (k=5 weighted)...")
tree = cKDTree(X_sub)
rest_mask = np.ones(n_total, dtype=bool)
rest_mask[sub_idx] = False
rest_idx = np.where(rest_mask)[0]
dists, knn_idx = tree.query(X[rest_idx], k=5)
weights = 1.0 / (dists + 1e-10)
weights /= weights.sum(axis=1, keepdims=True)
Y2 = np.empty((n_total, 2), dtype=np.float32)
Y3 = np.empty((n_total, 3), dtype=np.float32)
Y2[sub_idx] = Y2_sub
Y3[sub_idx] = Y3_sub
Y2[rest_idx] = np.einsum("nk,nkd->nd", weights, Y2_sub[knn_idx]).astype(np.float32)
Y3[rest_idx] = np.einsum("nk,nkd->nd", weights, Y3_sub[knn_idx]).astype(np.float32)
np.save("{out2_path}", Y2)
np.save("{out3_path}", Y3)
print(f"  PaCMAP done ({{time.time()-t0:.1f}}s, {{n_sub:,}}-cell fit)")
'''
                    else:
                        _do_native_pm = native_2d
                        script = f'''
import numpy as np, pacmap, time
X = np.load("{input_path}")
nn = min(10, max(2, X.shape[0] // 50))
t0 = time.time()
do_native_2d = {_do_native_pm}
if do_native_2d:
    print(f"  PaCMAP 2D ({{X.shape[0]:,}} cells, {{nn}} neighbors)...")
    Y2 = pacmap.PaCMAP(n_components=2, n_neighbors=nn, verbose=False).fit_transform(X).astype(np.float32)
    np.save("{out2_path}", Y2)
print(f"  PaCMAP 3D{{'' if do_native_2d else ' (2D will be derived)'}}")
Y3 = pacmap.PaCMAP(n_components=3, n_neighbors=nn, verbose=False).fit_transform(X).astype(np.float32)
np.save("{out3_path}", Y3)
if not do_native_2d:
    np.save("{out2_path}", Y3[:, :2].copy())
print(f"  PaCMAP done ({{time.time()-t0:.1f}}s)")
'''
                    sub_label = f" (subsample {PACMAP_SUBSAMPLE:,})" if PACMAP_SUBSAMPLE > 0 and n_total > PACMAP_SUBSAMPLE else ""
                    print(f"  Computing PaCMAP in background subprocess ({N_CELLS:,} cells{sub_label})...")
                    proc = subprocess.run(
                        [sys.executable, "-c", script],
                        capture_output=True, text=True
                    )
                    if proc.stdout:
                        for line in proc.stdout.strip().split("\n"):
                            print(line)
                    if proc.returncode != 0:
                        stderr = (proc.stderr or "").strip()[-300:]
                        print(f"  PaCMAP subprocess failed (exit {proc.returncode}): {stderr}")
                        PACMAP_COMPUTING = False
                        return

                    PACMAP_2D = np.load(out2_path)
                    PACMAP_3D = np.load(out3_path)
                    src = "computed_quick" if _pm_suffix == "_quick" else "computed_full"
                    EMBEDDING_SOURCES["pacmap_2d"] = src
                    EMBEDDING_SOURCES["pacmap_3d"] = src

                # Save to cache with method-specific keys
                cached[f"pacmap_2d{_pm_suffix}"] = PACMAP_2D
                cached[f"pacmap_3d{_pm_suffix}"] = PACMAP_3D
                try:
                    np.savez(_cache_write_path(".fleek_cache.npz"), **cached)
                    print(f"  PaCMAP cached")
                except Exception:
                    pass
            except Exception as e:
                print(f"  PaCMAP background compute failed: {e}")
            finally:
                PACMAP_COMPUTING = False
        threading.Thread(target=_bg_pacmap_subprocess, daemon=True).start()
    else:
        PACMAP_2D = None
        PACMAP_3D = None
        PACMAP_COMPUTING = False

    # Build fast lookup structures
    GENE_INDEX = {g: i for i, g in enumerate(ADATA.var_names)}

    # Detect HVG names (instant if pre-annotated, fast subsample if not)
    # Used for UI chips only — gene lookups use full CSC once built
    HVG_NAMES = []
    HVG_VAR = {}
    ALL_GENE_VAR = None
    ALL_GENE_VAR_METRIC = "dispersions_norm"
    GENE_CUTOFF_FLAG = None
    _has_hvg_col = "highly_variable" in ADATA.var.columns
    if _has_hvg_col:
        hvg_mask = ADATA.var["highly_variable"].values.astype(bool)
        HVG_NAMES = [str(ADATA.var_names[i]) for i in np.where(hvg_mask)[0]]
        # Try to get dispersion/variance from adata.var (for ALL genes, not just HVGs)
        _var_col = None
        for vcol in ["dispersions_norm", "dispersions", "variances_norm", "variances"]:
            if vcol in ADATA.var.columns:
                _var_col = vcol
                ALL_GENE_VAR_METRIC = vcol
                vals = ADATA.var[vcol].values
                ALL_GENE_VAR = np.asarray(vals, dtype=np.float32)
                for i in np.where(hvg_mask)[0]:
                    HVG_VAR[str(ADATA.var_names[i])] = float(vals[i])
                print(f"  Found {len(HVG_NAMES):,} pre-annotated HVGs (variability from '{vcol}', all {len(vals):,} genes)")
                HVG_NAMES.sort(key=lambda g: HVG_VAR.get(g, 0), reverse=True)
                break
        else:
            print(f"  Found {len(HVG_NAMES):,} pre-annotated HVGs (no variability column)")

    if not BACKED and sp.issparse(ADATA.X):
        CSC_BUILDING = True
        _csc_cr, _ = _cache_read_path(".csc_cache.npz")
        csc_cache_path = _csc_cr if _csc_cr else _cache_write_path(".csc_cache.npz")
        def _build_gene_index():
            global HVG_NAMES, HVG_VAR, ALL_GENE_VAR, ALL_GENE_VAR_METRIC, GENE_CUTOFF_FLAG, X_CSC, CSC_BUILDING, CSC_CACHED, CSC_TIME
            _csc_t0 = time.time()

            # ── Gene Variability: Locally Standardized Trend Residual ──
            #
            # Goal: Rank genes by biological variability, corrected for both
            # the mean-variance trend AND the expression-dependent noise level.
            #
            # PROBLEM 1 (mean-variance dependence): In scRNA-seq, variance
            # scales with mean expression. Housekeeping genes (ACTB, Malat1,
            # ribosomal proteins) dominate raw variance rankings simply because
            # they're highly expressed, not because they're biologically
            # informative.
            #
            # PROBLEM 2 (heteroscedastic residuals): Even after correcting for
            # the mean-variance trend, the *scatter* around the trend is not
            # uniform — it's much wider at low expression than at moderate
            # expression. A gene expressed in 2/50,000 cells can produce a raw
            # residual of +2.0 from pure Poisson noise, while a residual of
            # +2.0 at moderate expression would be genuinely remarkable. Simple
            # trend residuals therefore remain biased toward barely-detected genes.
            #
            # Common approaches and their limitations:
            #   - Raw variance: biased toward highly expressed genes
            #   - Coefficient of variation (CV = σ/μ): biased toward barely-
            #     detected genes (one stochastic count → enormous CV)
            #   - Binned z-score (scanpy's dispersions_norm): discretizes the
            #     smooth mean-variance curve into 20 bins and z-scores within
            #     each. Requires arbitrary cutoffs (min_mean, max_mean, min_frac)
            #     to exclude noisy genes — ad hoc and can silently discard signal
            #   - Simple polynomial residual: corrects the trend but not the
            #     heteroscedastic scatter — low-expression noise still inflates
            #     rankings for barely-detected genes
            #
            # OUR APPROACH: Two-stage polynomial fitting with local standardization.
            #
            # Stage 1 — Trend correction:
            #   1. Compute per-gene mean (μ) and variance (σ²) across cells
            #   2. Compute Fano factor: F = σ²/μ for each gene with μ > 0
            #   3. Transform to log-log space: log₁₀(μ) vs log₁₀(F)
            #   4. Fit a degree-3 polynomial: log₁₀(F̂) = p₁(log₁₀(μ))
            #   5. Raw residual r = log₁₀(F) − log₁₀(F̂)
            #
            # Stage 2 — Local standardization:
            #   6. Fit a degree-3 polynomial to |r| vs log₁₀(μ) to estimate
            #      the local scatter envelope: σ̂(μ) = p₂(log₁₀(μ))
            #   7. Final score z = r / max(σ̂, floor)
            #      where floor = 0.1 prevents division by near-zero scatter
            #      at expression levels where genes sit tightly on the trend
            #
            # WHY THIS WORKS:
            #   - No arbitrary cutoffs: every expressed gene gets a z-score
            #   - Housekeeping genes: high variance fully explained by trend → z ≈ 0
            #   - Barely-detected genes: large raw residual, but divided by the
            #     equally large local scatter → z ≈ 0 (noise, not signal)
            #   - True markers: high residual at moderate expression where scatter
            #     is tight → z >> 0 (signal stands out from noise)
            #   - Adapts to any normalization (raw counts, log1p, SCTransform)
            #   - Two np.polyfit calls on ~15k genes: <2ms total
            #
            if not _has_hvg_col:
                print("  Computing gene variability (locally standardized trend residual)...")
                t0 = time.time()
                import scipy.sparse as _sp
                # Route to the ideal matrix (raw counts > log-normalized > denoised).
                route = MATRIX_ROUTING.get("variability", {}) if MATRIX_ROUTING else {}
                src_path = route.get("path") or "X"
                X = _get_matrix_by_path(src_path) if src_path else ADATA.X
                if X is None:
                    X = ADATA.X
                    src_path = "X"
                print(f"  Variability: reading from {route.get('matrix_label') or src_path}")
                n_sub = min(50000, X.shape[0])
                if n_sub < X.shape[0]:
                    rng = np.random.default_rng(42)
                    idx_sub = rng.choice(X.shape[0], n_sub, replace=False)
                    X_sub = X[idx_sub]
                else:
                    X_sub = X

                # Step 1: Compute per-gene mean, variance, and detection rate
                if _sp.issparse(X_sub):
                    mean = np.asarray(X_sub.mean(axis=0)).ravel()
                    sq_mean = np.asarray(X_sub.multiply(X_sub).mean(axis=0)).ravel()
                    var = sq_mean - mean ** 2
                    n_expressing = np.asarray((X_sub > 0).sum(axis=0)).ravel()
                else:
                    X_dense = np.asarray(X_sub)
                    mean = np.mean(X_dense, axis=0)
                    var = np.var(X_dense, axis=0)
                    n_expressing = np.sum(X_dense > 0, axis=0)
                del X_sub

                mean = mean.astype(np.float64)
                var = var.astype(np.float64)

                # Step 2: Fano factor (F = variance / mean) for expressed genes
                fano = np.zeros_like(mean)
                pos = mean > 0
                fano[pos] = var[pos] / mean[pos]

                # Step 3: Log-log transform for trend fitting
                fit_mask = pos & (fano > 0)
                log_mean = np.full_like(mean, np.nan)
                log_fano = np.full_like(mean, np.nan)
                log_mean[fit_mask] = np.log10(mean[fit_mask])
                log_fano[fit_mask] = np.log10(fano[fit_mask])

                fit_idx = np.where(fit_mask)[0]
                scores = np.zeros(len(mean), dtype=np.float64)

                if len(fit_idx) > 10:
                    # Step 4: Fit degree-3 polynomial to the mean-Fano trend
                    coeffs_trend = np.polyfit(log_mean[fit_idx], log_fano[fit_idx], deg=3)
                    predicted = np.polyval(coeffs_trend, log_mean[fit_idx])

                    # Step 5: Raw residuals
                    raw_residuals = log_fano[fit_idx] - predicted

                    # Step 6: Fit degree-3 polynomial to |residual| vs log₁₀(μ)
                    # This estimates the local scatter envelope — how much
                    # residual noise is *expected* at each expression level.
                    # At low expression, scatter is wide (Poisson noise);
                    # at moderate expression, scatter is tight (signal stands out)
                    abs_resid = np.abs(raw_residuals)
                    coeffs_scatter = np.polyfit(log_mean[fit_idx], abs_resid, deg=3)
                    local_sigma = np.polyval(coeffs_scatter, log_mean[fit_idx])

                    # Step 7: Locally standardized z-score
                    # z = raw_residual / local_σ
                    # Floor of 0.1 prevents division by near-zero at expression
                    # levels where genes happen to sit tightly on the trend
                    local_sigma = np.maximum(local_sigma, 0.1)
                    scores[fit_idx] = raw_residuals / local_sigma

                ALL_GENE_VAR = scores.astype(np.float32)
                ALL_GENE_VAR_METRIC = "trend_residual"

                # Compute which genes WOULD have been excluded by traditional
                # scanpy-style cutoffs (min_mean, max_mean, min_frac).
                # Flagged red in the UI so the researcher can evaluate whether
                # the locally standardized scores handle them correctly.
                # Stored as a bitmask so the UI can show specific reasons:
                #   bit 0 (1) = not expressed (mean == 0)
                #   bit 1 (2) = low mean (0 < mean < 0.0125)
                #   bit 2 (4) = high mean (mean > adaptive max)
                #   bit 3 (8) = low fraction (< 0.5% of cells expressing)
                frac_expr = n_expressing / n_sub
                mean_95 = np.percentile(mean[pos], 95) if pos.sum() > 0 else 3.0
                _max_mean = max(3.0, mean_95)
                flags = np.zeros(len(mean), dtype=np.int8)
                flags[~pos] |= 1                            # not expressed
                flags[pos & (mean < 0.0125)] |= 2           # low mean
                flags[pos & (mean > _max_mean)] |= 4        # high mean
                flags[frac_expr < 0.005] |= 8               # low fraction
                GENE_CUTOFF_FLAG = flags
                n_flagged = int((flags > 0).sum())

                n_hvg = min(2000, int(pos.sum()))
                hvg_idx = np.argsort(ALL_GENE_VAR)[-n_hvg:][::-1]
                HVG_NAMES = [str(ADATA.var_names[i]) for i in hvg_idx]
                HVG_VAR = {str(ADATA.var_names[i]): float(ALL_GENE_VAR[i]) for i in hvg_idx}
                n_scored = int(fit_mask.sum())
                print(f"  Scored {n_scored:,}/{len(mean):,} genes by locally standardized trend residual, {n_flagged:,} flagged by cutoffs ({time.time()-t0:.1f}s, {n_sub:,}-cell subsample)")

            # Try loading CSC from cache
            import scipy.sparse as _sp
            loaded_from_cache = False
            if csc_cache_path.exists():
                print("  Loading cached CSC index...")
                try:
                    t1 = time.time()
                    X_CSC = _sp.load_npz(str(csc_cache_path))
                    if X_CSC.shape != ADATA.X.shape:
                        print(f"  CSC cache shape mismatch ({X_CSC.shape} vs {ADATA.X.shape}), rebuilding...")
                        X_CSC = None
                    else:
                        loaded_from_cache = True
                        CSC_CACHED = True
                        print(f"  CSC index loaded from cache ({time.time()-t1:.1f}s)")
                except Exception as e:
                    print(f"  CSC cache load failed ({e}), rebuilding...")
                    X_CSC = None

            if not loaded_from_cache:
                # Build full CSC — single pass, no intermediate copies
                print("  Building CSC column index...")
                t1 = time.time()
                X_CSC = ADATA.X.tocsc()
                elapsed = time.time() - t1
                print(f"  CSC index ready ({elapsed:.0f}s) — all gene lookups now fast")
                # Save to cache (use write path, may differ from read path)
                try:
                    csc_wp = _cache_write_path(".csc_cache.npz")
                    _sp.save_npz(str(csc_wp), X_CSC, compressed=False)
                    cache_mb = os.path.getsize(str(csc_wp)) / 1e6
                    print(f"  CSC cache saved ({cache_mb:.0f} MB)")
                except Exception as e:
                    print(f"  CSC cache save failed ({e}) — will rebuild next time")

            CSC_TIME = round(time.time() - _csc_t0, 1)
            CSC_BUILDING = False
        threading.Thread(target=_build_gene_index, daemon=True).start()
    else:
        X_CSC = None
        CSC_BUILDING = False

    # Record load settings for the client
    load_elapsed = time.time() - t0
    pacmap_suffix = "_quick" if PACMAP_SETTINGS.get("fast", True) else "_full"
    LOAD_SETTINGS.clear()
    LOAD_SETTINGS["backed"] = BACKED
    LOAD_SETTINGS["quickmap"] = use_fast
    LOAD_SETTINGS["quickmap_n"] = fast_umap_subsample if use_fast else 0
    LOAD_SETTINGS["quickdeg"] = DEG_SETTINGS.get("fast", True)
    LOAD_SETTINGS["quickdeg_n"] = DEG_SETTINGS.get("max_cells", 50000) if DEG_SETTINGS.get("fast", True) else 0
    LOAD_SETTINGS["quickpacmap"] = PACMAP_SETTINGS.get("fast", True)
    LOAD_SETTINGS["native_2d"] = native_2d
    LOAD_SETTINGS["umap_cached"] = umap_from_cache
    LOAD_SETTINGS["umap_method"] = "quick" if use_fast else "full"
    LOAD_SETTINGS["pca_cached"] = ("pca_full" in cached or "pca_2d" in cached)
    LOAD_SETTINGS["leiden_cached"] = ("leiden_ids" in cached)
    LOAD_SETTINGS["pacmap_cached"] = (f"pacmap_2d{pacmap_suffix}" in cached or "pacmap_2d" in cached)
    LOAD_SETTINGS["pacmap_method"] = "quick" if PACMAP_SETTINGS.get("fast", True) else "full"
    LOAD_SETTINGS["csc_cached"] = CSC_CACHED
    LOAD_SETTINGS["load_time"] = round(load_elapsed, 1)
    LOAD_SETTINGS["cache_time"] = round(_cache_time, 2)
    LOAD_SETTINGS["cache_mode"] = CACHE_MODE
    LOAD_SETTINGS["cache_dir"] = str(CACHE_DIR) if CACHE_DIR else ""
    # Check annotation caches
    for key, suffix in [("annot_markers_cached", ".annot_markers.json"), ("annot_llm_cached", ".annot_llm.json"), ("annot_lineage_cached", ".annot_lineage.json")]:
        ap, _ = _cache_read_path(suffix)
        ap_exists = ap is not None
        LOAD_SETTINGS[key] = ap_exists

    # Load marker database for cell type annotation
    load_marker_db()
    load_go_db()
    detect_counts()

    PROGRESS["status"] = "done"
    PROGRESS["pct"] = 100
    PROGRESS["message"] = f"Ready: {N_CELLS:,} cells, {len(CLUSTER_NAMES)} types, {len(GENE_NAMES_LIST):,} genes"
    print(f"  Ready: {N_CELLS:,} cells, {len(CLUSTER_NAMES)} types, {len(GENE_NAMES_LIST):,} genes")


def _stratified_subsample(cluster_ids, n_target, rng):
    """Stratified subsample preserving cluster proportions, with minimum per cluster."""
    unique = np.unique(cluster_ids)
    n_total = len(cluster_ids)
    indices = []
    min_per_cluster = max(5, n_target // (len(unique) * 4))  # floor per cluster

    for cid in unique:
        mask = np.where(cluster_ids == cid)[0]
        # Proportional share
        share = max(min_per_cluster, int(round(len(mask) / n_total * n_target)))
        share = min(share, len(mask))  # can't take more than exist
        chosen = rng.choice(mask, size=share, replace=False)
        indices.append(chosen)

    result = np.concatenate(indices)
    # If we overshot, trim randomly
    if len(result) > n_target * 1.2:
        result = rng.choice(result, size=n_target, replace=False)
    result.sort()
    return result


def _ensure_neighbors(adata):
    """Build neighbor graph if missing."""
    import scanpy as sc

    has_conn = "connectivities" in adata.obsp
    has_uns = "neighbors" in adata.uns

    if has_conn and has_uns:
        return

    # If we have connectivities but no uns metadata, synthesize it
    if has_conn and not has_uns:
        adata.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances",
                                  "params": {"n_neighbors": 15, "method": "umap"}}
        print("    Synthesized missing neighbors metadata from existing connectivities")
        return

    if "scvi" in adata.obsm:
        emb = adata.obsm["scvi"]
        nan_mask = np.isnan(emb).any(axis=1)
        if nan_mask.sum() > 0:
            print(f"    Removing {nan_mask.sum()} NaN embedding cells")
            # Can't easily remove in place, just zero them
            emb[nan_mask] = 0
        print("    Building neighbors from scVI embedding...")
        sc.pp.neighbors(adata, use_rep="scvi", n_neighbors=15)
    elif "X_pca" in adata.obsm:
        print("    Building neighbors from existing PCA...")
        sc.pp.neighbors(adata, n_pcs=50, n_neighbors=15)
    else:
        print("    Running PCA + neighbors from scratch...")
        tmp = adata.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        # Use cell_ranger flavor — doesn't require skmisc
        sc.pp.highly_variable_genes(tmp, n_top_genes=2000)
        tmp = tmp[:, tmp.var.highly_variable].copy()
        sc.pp.scale(tmp, max_value=10)
        sc.tl.pca(tmp, n_comps=50)
        sc.pp.neighbors(tmp, n_pcs=50, n_neighbors=15)
        adata.obsp["distances"] = tmp.obsp["distances"]
        adata.obsp["connectivities"] = tmp.obsp["connectivities"]
        adata.uns["neighbors"] = tmp.uns["neighbors"]


_CELL_LIB_SIZES = None  # lazily computed per-cell totals for normalize+log1p

def _get_lib_sizes(src_mat):
    """Return per-cell library sizes for the given matrix, caching the result so
    repeated gene lookups don't rescan the matrix. Keyed by `id(src_mat)`."""
    global _CELL_LIB_SIZES
    if _CELL_LIB_SIZES is not None and _CELL_LIB_SIZES[0] == id(src_mat):
        return _CELL_LIB_SIZES[1]
    lib = np.asarray(src_mat.sum(axis=1)).ravel().astype(np.float64)
    lib[lib == 0] = 1.0
    _CELL_LIB_SIZES = (id(src_mat), lib)
    return lib


def get_gene_expression(gene_name):
    """Get normalized 0-1 expression for a single gene, routed to the best
    available matrix per MATRIX_ROUTING['gene_display']. Applies normalize+log1p
    on the fly if the chosen matrix is raw/denoised counts."""
    import scipy.sparse as sp

    idx = GENE_INDEX.get(gene_name)
    if idx is None:
        return None

    # Which matrix + transform are we using for display?
    route = MATRIX_ROUTING.get("gene_display", {}) if MATRIX_ROUTING else {}
    src_path = route.get("path") or "X"
    transform = route.get("transform", "none")

    # Fast path: the CSC index was built from adata.X. If the route points to a
    # different matrix OR asks for a transform, we need to read raw values from
    # the routed matrix directly.
    use_csc = X_CSC is not None and src_path == "X" and transform == "none"

    if use_csc:
        col = X_CSC[:, idx].toarray().ravel().astype(np.float32)
    else:
        src_mat = _get_matrix_by_path(src_path) if src_path else None
        if src_mat is None:
            src_mat = ADATA.X
        # Read the column
        if BACKED and src_path == "X":
            # Can't use CSC when backed; fall back to per-batch reads of raw X.
            n = ADATA.n_obs
            col = np.empty(n, dtype=np.float32)
            batch = 50000
            for bi in range(0, n, batch):
                be = min(bi + batch, n)
                chunk = src_mat[bi:be, idx]
                if sp.issparse(chunk):
                    chunk = chunk.toarray().ravel()
                else:
                    chunk = np.asarray(chunk).ravel()
                col[bi:be] = chunk
        else:
            col_src = src_mat[:, idx]
            if sp.issparse(col_src):
                col = col_src.toarray().ravel()
            else:
                col = np.asarray(col_src).ravel()
            col = col.astype(np.float32)

        if transform == "normalize+log1p":
            lib = _get_lib_sizes(src_mat)
            col = np.asarray(col, dtype=np.float64)
            col = col / lib * 1e4
            col = np.log1p(col).astype(np.float32)

    # Normalize to 0-1 for the client's colormap.
    pos = col[col > 0]
    if len(pos) > 10:
        pmax = np.percentile(pos, 98)
    else:
        pmax = col.max()
    if pmax > 0:
        col = np.clip(col / pmax, 0, 1)
    return col


def build_init_payload():
    """Build the initial binary payload with coords + clusters."""
    parts = []

    avail = []
    if UMAP_2D is not None:
        avail.append(2)
    if UMAP_3D is not None:
        avail.append(3)

    # Build embedding manifest: {method: [dims]}
    emb_manifest = {}
    emb_manifest["umap"] = list(avail)
    pca_dims = []
    if PCA_2D is not None:
        pca_dims.append(2)
    if PCA_3D is not None:
        pca_dims.append(3)
    if pca_dims:
        emb_manifest["pca"] = pca_dims

    pacmap_dims = []
    if PACMAP_2D is not None:
        pacmap_dims.append(2)
    if PACMAP_3D is not None:
        pacmap_dims.append(3)
    if pacmap_dims:
        emb_manifest["pacmap"] = pacmap_dims

    header = {
        "version": 4,
        "n_cells": N_CELLS,
        "n_dims": max(avail) if avail else 2,
        "available_dims": avail,
        "available_embeddings": emb_manifest,
        "n_genes": 0,  # genes loaded on demand
        "n_clusters": len(CLUSTER_NAMES),
        "gene_names": GENE_NAMES_LIST[:500],  # send first 500 for search
        "total_genes": len(GENE_NAMES_LIST),
        "cluster_names": CLUSTER_NAMES,
        "cluster_colors": CLUSTER_COLORS,
        "metadata_columns": {},
        "obs_categories": {col: {"categories": v["categories"], "counts": v["counts"]}
                           for col, v in OBS_COLS.items()},
        "load_settings": dict(LOAD_SETTINGS),
        "loaded_path": LOADED_PATH,
        "hvg_genes": HVG_NAMES[:200],  # top HVGs for quick-access UI
        "hvg_var": {g: round(HVG_VAR[g], 4) for g in HVG_NAMES[:200] if g in HVG_VAR},
        "organism": _DETECTED_ORGANISM[0] if _DETECTED_ORGANISM else "Human",
        "organism_reason": _DETECTED_ORGANISM[1] if _DETECTED_ORGANISM else "",
        "go_available": GO_DB is not None,
        "go_organism": GO_DB.get("label", "") if GO_DB else "",
        "counts_available": COUNTS_MATRIX is not None and (COUNTS_INFO is None or COUNTS_INFO.get("pb_ok", True)),
        "counts_label": COUNTS_LABEL or "",
        "counts_info": COUNTS_INFO,
        "pydeseq2_available": PYDESEQ2_AVAILABLE,
        "matrix_registry": MATRIX_REGISTRY,
        "matrix_routing": MATRIX_ROUTING,
        "embedding_sources": EMBEDDING_SOURCES,
    }

    # Include cached annotations if available
    # Validation: cluster_names must be a superset of current names (handles subsets)
    def _annot_valid(cached_obj):
        cn = cached_obj.get("cluster_names")
        nc = cached_obj.get("n_cells")
        if cn == CLUSTER_NAMES and nc == N_CELLS:
            return True  # exact match
        if cn and set(CLUSTER_NAMES).issubset(set(cn)):
            return True  # subset of parent — annotations still valid
        return False

    # Parent stem fallback: subsets store their parent's stem in uns so we
    # can look up the parent's annotation caches directly.
    _parent_stem = None
    if ADATA is not None:
        _parent_stem = ADATA.uns.get("fleek_parent_stem")

    def _find_annot_cache(suffix):
        """Find annotation cache: try normal path, then parent stem fallback."""
        ap, _ = _cache_read_path(suffix)
        if ap:
            return ap
        if _parent_stem and CACHE_DIR:
            pp = CACHE_DIR / f"{_parent_stem}{suffix}"
            if pp.exists():
                return pp
        return None

    for annot_key, suffix in [("cached_markers", ".annot_markers.json"), ("cached_llm", ".annot_llm.json")]:
        if LOADED_PATH:
            ap = _find_annot_cache(suffix)
            if ap:
                try:
                    with open(ap) as f:
                        cached = json.load(f)
                    if _annot_valid(cached):
                        # Remap cluster_ids by name so a subset's shifted cluster
                        # indexing doesn't cross-wire annotations to wrong clusters.
                        cached = _remap_cached_annotations(cached)
                        header[annot_key] = cached.get("results", [])
                        _sfx = " (remapped for subset)" if cached.get("remapped_from_parent") else ""
                        print(f"  Init: loaded {annot_key} from {ap.name}{_sfx}")
                    else:
                        print(f"  Init: {ap.name} FAILED validation (cached_names={cached.get('cluster_names')}, current={CLUSTER_NAMES})")
                except Exception as e:
                    print(f"  Init: error reading {ap}: {e}")

    # Include cached lineage tree if available
    if LOADED_PATH:
        lp = _find_annot_cache(".annot_lineage.json")
        if lp:
            try:
                with open(lp) as f:
                    lcached = json.load(f)
                if _annot_valid(lcached):
                    header["cached_lineage"] = lcached.get("tree")
                    print(f"  Init: loaded cached_lineage from {lp.name}")
                else:
                    print(f"  Init: {lp.name} lineage FAILED validation (cached_names={lcached.get('cluster_names')}, current={CLUSTER_NAMES})")
            except Exception as e:
                print(f"  Init: error reading lineage {lp}: {e}")

    header_bytes = json.dumps(_json_safe(header), separators=(",", ":")).encode("utf-8")

    buf = io.BytesIO()
    buf.write(b"SCRN")
    # Pad header to 4-byte alignment
    while len(header_bytes) % 4 != 0:
        header_bytes += b" "
    buf.write(struct.pack("<I", len(header_bytes)))
    buf.write(header_bytes)

    # Write embeddings in manifest order: umap first, then pca, etc.
    _emb_arrays = {
        "umap": {2: UMAP_2D, 3: UMAP_3D},
        "pca": {2: PCA_2D, 3: PCA_3D},
        "pacmap": {2: PACMAP_2D, 3: PACMAP_3D},
    }
    for method, dims in emb_manifest.items():
        for nd in dims:
            arr = _emb_arrays.get(method, {}).get(nd)
            if arr is not None:
                for d in range(nd):
                    buf.write(arr[:, d].astype(np.float32).tobytes())

    buf.write(CLUSTER_IDS.tobytes())

    # Write obs column codes (Int16, one array per column, in obs_categories order)
    for col in OBS_COLS:
        buf.write(OBS_COLS[col]["codes"].tobytes())

    return buf.getvalue()


def run_deg(group_a, group_b, test="wilcoxon"):
    """Run DEG between two cell groups with Cohen's d.
    
    Speed optimization: filters to genes expressed in >=3 cells across both
    groups before testing.  On a 60k-gene dataset this typically reduces the
    test set to 10–20k genes → 3–6× faster.
    """
    import scanpy as sc
    import scipy.sparse as sp

    MIN_CELLS_EXPR = 3  # gene must be nonzero in >=3 cells to be tested

    t0 = time.time()
    indices_a = list(group_a)
    indices_b = list(group_b)
    indices = indices_a + indices_b
    labels = ["A"] * len(indices_a) + ["B"] * len(indices_b)

    # Pick the right source matrix via the routing table so we don't silently run
    # on denoised / scaled data when raw counts are available elsewhere.
    route = MATRIX_ROUTING.get("deg_sc", {}) if MATRIX_ROUTING else {}
    src_path = route.get("path") or "X"
    src_mat = _get_matrix_by_path(src_path) if src_path else None
    if src_mat is None:
        src_mat = ADATA.X
        src_path = "X"

    sub = ADATA[indices]
    if BACKED:
        try:
            sub = sub.to_memory()
        except AttributeError:
            sub = sub.copy()
    else:
        sub = sub.copy()
    sub.obs["deg_group"] = labels

    # If the chosen source isn't already the subset's X, materialize it now.
    if src_path != "X":
        try:
            src_sub = src_mat[indices, :]
            if not sp.issparse(src_sub) and hasattr(src_sub, "toarray"):
                src_sub = src_sub.toarray()
            sub.X = src_sub
            print(f"  DEG: using {route.get('matrix_label') or src_path} as source matrix")
        except Exception as _e:
            print(f"  DEG: failed to switch to {src_path} ({_e}), falling back to adata.X")

    transform = route.get("transform", "none")
    if transform == "normalize+log1p":
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)

    # ── Filter to expressed genes (big speedup) ──
    n_all_orig = sub.n_vars
    if sp.issparse(sub.X):
        nnz = np.asarray((sub.X != 0).sum(axis=0)).ravel()
    else:
        nnz = np.asarray((np.asarray(sub.X) != 0).sum(axis=0)).ravel()
    keep_mask = nnz >= MIN_CELLS_EXPR
    n_kept = int(keep_mask.sum())
    if n_kept < n_all_orig and n_kept > 0:
        sub = sub[:, keep_mask].copy()
        print(f"  DEG: filtered {n_all_orig:,} → {n_kept:,} expressed genes")

    t_prep = time.time() - t0

    method = "wilcoxon" if test == "wilcoxon" else "t-test"
    n_all = sub.n_vars
    na, nb = len(indices_a), len(indices_b)
    print(f"  DEG: testing {n_all:,} genes, {na} vs {nb} cells, "
          f"{method}... (prep {t_prep:.1f}s)")
    sc.tl.rank_genes_groups(sub, groupby="deg_group", groups=["A"],
                            reference="B", method=method, n_genes=n_all)

    result = sub.uns["rank_genes_groups"]
    names = [str(x) for x in result["names"]["A"]]
    logfc = [float(x) for x in result["logfoldchanges"]["A"]]
    pvals = [float(x) for x in result["pvals"]["A"]]
    padj = [float(x) for x in result["pvals_adj"]["A"]]

    # Cohen's d — computed on sparse directly, no densification
    print(f"  DEG: computing Cohen's d...")
    X_a = sub.X[:na]
    X_b = sub.X[na:]
    if sp.issparse(X_a):
        mean_a = np.asarray(X_a.mean(axis=0)).ravel()
        mean_b = np.asarray(X_b.mean(axis=0)).ravel()
        sq_a = X_a.multiply(X_a)
        sq_b = X_b.multiply(X_b)
        var_a = np.asarray(sq_a.mean(axis=0)).ravel() - mean_a**2
        var_b = np.asarray(sq_b.mean(axis=0)).ravel() - mean_b**2
        if na > 1: var_a *= na / (na - 1)
        if nb > 1: var_b *= nb / (nb - 1)
        del sq_a, sq_b
    else:
        X_a = np.asarray(X_a, dtype=np.float64)
        X_b = np.asarray(X_b, dtype=np.float64)
        mean_a = X_a.mean(axis=0)
        mean_b = X_b.mean(axis=0)
        var_a = X_a.var(axis=0, ddof=1) if na > 1 else np.zeros(n_all)
        var_b = X_b.var(axis=0, ddof=1) if nb > 1 else np.zeros(n_all)
    pooled_sd = np.sqrt(((na - 1) * var_a + (nb - 1) * var_b) / max(na + nb - 2, 1))
    cohens_d_all = np.where(pooled_sd > 0, (mean_a - mean_b) / pooled_sd, 0.0)

    gene_to_varidx = {str(g): i for i, g in enumerate(sub.var_names)}

    clean = []
    for i in range(len(names)):
        fc = logfc[i]
        p = padj[i]
        if not (np.isfinite(fc) and np.isfinite(p) and p > 0):
            continue
        vidx = gene_to_varidx.get(names[i])
        cd = float(cohens_d_all[vidx]) if vidx is not None else 0.0
        if not np.isfinite(cd):
            cd = 0.0
        clean.append({"gene": names[i], "log2fc": round(fc, 4),
                      "cohens_d": round(cd, 4), "pval": pvals[i], "padj": p})

    clean.sort(key=lambda x: x["padj"])
    elapsed = time.time() - t0
    n_sig = sum(1 for g in clean if abs(g["log2fc"]) >= 0.5 and g["padj"] < 0.05)
    n_filt_note = f" (filtered from {n_all_orig:,})" if n_kept < n_all_orig else ""
    print(f"  DEG: done in {elapsed:.1f}s — {len(clean):,} genes{n_filt_note}, {n_sig:,} significant")

    return {
        "genes": clean, "test": method, "n_a": na, "n_b": nb,
        "n_total_genes": n_all, "n_total_genes_unfiltered": n_all_orig,
        "elapsed": round(elapsed, 1),
        "note": "P-values from single-sample comparison — pseudoreplication caveat applies."
    }


def export_h5ad_subset(cell_indices, group_name):
    """Export a subset of cells as h5ad, return file path.

    Embeddings from the full dataset are sliced and stored in obsm so the
    subset loads instantly without recomputing UMAP / PaCMAP / PCA.
    """
    import tempfile
    sub = ADATA[cell_indices]
    if BACKED:
        try:
            sub = sub.to_memory()
        except AttributeError:
            sub = sub.copy()
    else:
        sub = sub.copy()

    # ── Embed sliced embeddings so the subset loads without recomputing ──
    # N_CELLS is the size of the currently loaded dataset; all global embedding
    # arrays must match it, so we use it as the sanity check.
    idx = np.asarray(cell_indices)
    _nc = N_CELLS  # capture once

    # UMAP: store 3D if available (obsm fallback slices [:, :2] / [:, :3])
    if UMAP_3D is not None and UMAP_3D.shape[0] == _nc:
        sub.obsm["X_umap"] = UMAP_3D[idx]
    elif UMAP_2D is not None and UMAP_2D.shape[0] == _nc:
        sub.obsm["X_umap"] = UMAP_2D[idx]

    # PaCMAP: same strategy
    if PACMAP_3D is not None and PACMAP_3D.shape[0] == _nc:
        sub.obsm["X_pacmap"] = PACMAP_3D[idx]
    elif PACMAP_2D is not None and PACMAP_2D.shape[0] == _nc:
        sub.obsm["X_pacmap"] = PACMAP_2D[idx]

    # PCA: adata.obsm["X_pca"] is already preserved by copy() if it was set,
    # but if it only lived in the cache (not obsm), add it now.
    if "X_pca" not in sub.obsm:
        if PCA_3D is not None and PCA_3D.shape[0] == _nc:
            sub.obsm["X_pca"] = PCA_3D[idx]
        elif PCA_2D is not None and PCA_2D.shape[0] == _nc:
            sub.obsm["X_pca"] = PCA_2D[idx]

    # ── Bake colors + annotations into the subset ─────────────────────────
    # Store name→color pairs so colors survive subsetting (positional lists
    # break when categories are removed by subsetting). Use parallel arrays
    # instead of a plain dict because HDF5 can't have "/" in group-member
    # names and anndata silently drops any dict key containing a slash on
    # write — which would discard cluster names like "Paneth-like/DCS".
    if CLUSTER_NAMES and CLUSTER_COLORS:
        n_pairs = min(len(CLUSTER_NAMES), len(CLUSTER_COLORS))
        sub.uns["fleek_color_map"] = {
            "names": [str(CLUSTER_NAMES[i]) for i in range(n_pairs)],
            "colors": [str(CLUSTER_COLORS[i]) for i in range(n_pairs)],
        }

    # Store parent stem so subset can find parent's annotation caches
    _ps = Path(LOADED_PATH).stem
    sub.uns["fleek_parent_stem"] = _ps
    print(f"  Export: stored fleek_parent_stem = {_ps!r}")

    n_emb = sum(k in sub.obsm for k in ("X_umap", "X_pacmap", "X_pca"))
    print(f"  Export: {len(cell_indices)} cells, {n_emb} embeddings in obsm")

    safe_name = "".join(c if c.isalnum() or c in "._-" else "_" for c in group_name)
    fname = f"subset_{safe_name}_{len(cell_indices)}cells.h5ad"
    fpath = os.path.join(tempfile.gettempdir(), fname)
    sub.write_h5ad(fpath)
    print(f"  Export: {len(cell_indices)} cells -> {fpath}")

    # ── Pre-write cache to server cache dir ──────────────────────────────────
    # The exported h5ad lands wherever the user saves it (unknown to the server),
    # so we can't write a .fleek_cache.npz next to it at download time.
    # Instead, we pre-write to ~/.fleek_cache/ keyed by the subset filename.
    # When the subset is loaded from any location, _cache_read_path falls back
    # to the server cache dir and finds it instantly — no recomputation needed.
    try:
        _init_cache_dir()
        if CACHE_DIR is not None:
            stem = fname[:-5]  # strip ".h5ad"
            server_cache = CACHE_DIR / f"{stem}.fleek_cache.npz"
            subset_cache = {}
            # UMAP — use legacy keys so they're found regardless of quick/full suffix
            if UMAP_3D is not None and UMAP_3D.shape[0] == _nc:
                subset_cache["umap_2d"] = UMAP_3D[idx, :2]
                subset_cache["umap_3d"] = UMAP_3D[idx]
            elif UMAP_2D is not None and UMAP_2D.shape[0] == _nc:
                subset_cache["umap_2d"] = UMAP_2D[idx]
            # PaCMAP
            if PACMAP_3D is not None and PACMAP_3D.shape[0] == _nc:
                subset_cache["pacmap_2d"] = PACMAP_3D[idx, :2]
                subset_cache["pacmap_3d"] = PACMAP_3D[idx]
            elif PACMAP_2D is not None and PACMAP_2D.shape[0] == _nc:
                subset_cache["pacmap_2d"] = PACMAP_2D[idx]
            # PCA — prefer the full-dimensional version from obsm if available
            pca_src = sub.obsm.get("X_pca")
            if pca_src is not None:
                pca_src = pca_src.astype(np.float32)
                subset_cache["pca_full"] = pca_src
                subset_cache["pca_2d"] = pca_src[:, :2]
                subset_cache["pca_3d"] = pca_src[:, :3] if pca_src.shape[1] >= 3 else pca_src
            if subset_cache:
                np.savez(server_cache, **subset_cache)
                print(f"  Export: pre-wrote cache ({', '.join(subset_cache)}) -> {server_cache}")
            # Copy annotation caches keyed by PARENT stem so subsets can
            # find them via fleek_parent_stem regardless of file renaming
            import shutil
            _parent = Path(LOADED_PATH).stem
            for _annot_suffix in [".annot_markers.json", ".annot_llm.json", ".annot_lineage.json"]:
                _annot_src, _ = _cache_read_path(_annot_suffix)
                if _annot_src:
                    try:
                        _annot_dst = CACHE_DIR / f"{_parent}{_annot_suffix}"
                        if not _annot_dst.exists() or str(_annot_src) != str(_annot_dst):
                            shutil.copy2(str(_annot_src), str(_annot_dst))
                        print(f"  Export: ensured {_annot_suffix} at {_annot_dst.name}")
                    except Exception as _ae:
                        print(f"  Export: failed to copy {_annot_suffix}: {_ae}")
                else:
                    print(f"  Export: no source found for {_annot_suffix}")
    except Exception as e:
        print(f"  Export: cache pre-write failed (non-fatal): {e}")

    return fpath, fname



# ─── HTTP Server ───

class FleekHandler(SimpleHTTPRequestHandler):
    """Serves the visualizer HTML and API endpoints."""

    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path
        params = parse_qs(parsed.query)

        if path == "/" or path == "/index.html":
            self._serve_html()
        elif path == "/api/init":
            self._serve_init()
        elif path == "/api/gene":
            gene = params.get("name", [""])[0]
            self._serve_gene(gene)
        elif path == "/api/search":
            q = params.get("q", [""])[0]
            self._serve_search(q)
        elif path == "/api/go-search":
            q = params.get("q", [""])[0]
            self._serve_go_search(q)
        elif path == "/api/go-genes":
            go_id = params.get("id", [""])[0]
            self._serve_go_genes(go_id)
        elif path == "/api/go-gene-terms":
            gene = params.get("gene", [""])[0]
            self._serve_go_gene_terms(gene)
        elif path == "/api/progress":
            self._serve_progress()
        elif path == "/api/embedding":
            self._serve_embedding(params.get("method", [""])[0])
        elif path == "/api/browse":
            self._serve_browse(params.get("path", [""])[0])
        elif path == "/api/key-status":
            self._serve_key_status()
        elif path == "/api/complete":
            self._serve_complete(params.get("partial", [""])[0])
        elif path == "/api/gene-var":
            self._serve_gene_var()
        elif path == "/api/cluster-genes":
            self._serve_cluster_genes(params)
        else:
            self.send_error(404)

    def do_POST(self):
        parsed = urlparse(self.path)
        path = parsed.path
        if path == "/api/upload":
            self._serve_upload()
            return
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length)
        req = json.loads(body)
        if path == "/api/deg":
            self._serve_deg(req)
        elif path == "/api/pseudobulk-deg":
            self._serve_pseudobulk_deg(req)
        elif path == "/api/export":
            self._serve_export(req)
        elif path == "/api/load":
            self._serve_load(req)
        elif path == "/api/annotate":
            self._serve_annotate(req)
        elif path == "/api/annotate-llm":
            self._serve_annotate_llm(req)
        elif path == "/api/compute-embedding":
            self._serve_compute_embedding(req)
        elif path == "/api/abort":
            self._serve_abort(req)
        elif path == "/api/unload":
            self._serve_unload(req)
        elif path == "/api/heartbeat":
            self._serve_heartbeat(req)
        elif path == "/api/lineage":
            self._serve_lineage(req)
        elif path == "/api/clear-cache":
            self._serve_clear_cache(req)
        elif path == "/api/cache-settings":
            self._serve_cache_settings(req)
        elif path == "/api/set-key":
            self._serve_set_key(req)
        elif path == "/api/clear-key":
            self._serve_clear_key()
        elif path == "/api/test-key":
            self._serve_test_key(req)
        elif path == "/api/sessions":
            self._serve_sessions()
        elif path == "/api/session-save":
            self._serve_session_save(req)
        elif path == "/api/session-load":
            self._serve_session_load(req)
        elif path == "/api/session-delete":
            self._serve_session_delete(req)
        else:
            self.send_error(404)

    def _serve_html(self):
        html_path = Path(__file__).parent / "fleek.html"
        if not html_path.exists():
            html_path = _BUNDLE_DIR / "rna_fleek" / "fleek.html"
        if not html_path.exists():
            html_path = _BUNDLE_DIR / "fleek.html"
        if not html_path.exists():
            self.send_error(500, "fleek.html not found next to server script")
            return
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate, max-age=0")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        self.end_headers()
        self.wfile.write(html_path.read_bytes())

    def _serve_init(self):
        if ADATA is None:
            resp = json.dumps({"no_data": True}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(resp)))
            self.end_headers()
            self.wfile.write(resp)
            return
        data = build_init_payload()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(data)

    def _serve_gene(self, gene_name):
        if ADATA is None:
            self.send_error(503, "No dataset loaded");return
        expr = get_gene_expression(gene_name)
        if expr is None:
            self.send_error(404, f"Gene '{gene_name}' not found")
            return
        data = expr.tobytes()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_embedding(self, method):
        """Return embedding coordinates as binary float32."""
        _emb_arrays = {
            "umap": {2: UMAP_2D, 3: UMAP_3D},
            "pca": {2: PCA_2D, 3: PCA_3D},
            "pacmap": {2: PACMAP_2D, 3: PACMAP_3D},
        }
        emb = _emb_arrays.get(method)
        if not emb:
            err = json.dumps({"error": f"Unknown method: {method}"}).encode()
            self.send_response(404)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)
            return
        dims = []
        if emb.get(2) is not None: dims.append(2)
        if emb.get(3) is not None: dims.append(3)
        if not dims:
            err = json.dumps({"error": f"{method} not computed yet"}).encode()
            self.send_response(404)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)
            return
        # Binary format: 4-byte header_len, JSON header, then float32 arrays
        header = {"method": method, "dims": dims, "n_cells": N_CELLS}
        header_bytes = json.dumps(header, separators=(",", ":")).encode()
        while len(header_bytes) % 4 != 0:
            header_bytes += b" "
        buf = io.BytesIO()
        buf.write(struct.pack("<I", len(header_bytes)))
        buf.write(header_bytes)
        for nd in dims:
            arr = emb[nd]
            for d in range(nd):
                buf.write(arr[:, d].astype(np.float32).tobytes())
        data = buf.getvalue()
        self.send_response(200)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_browse(self, req_path):
        """List directories and h5ad files for the file browser."""
        try:
            # Default to directory of loaded file, or cwd
            if not req_path:
                if LOADED_PATH:
                    req_path = str(Path(LOADED_PATH).parent)
                else:
                    req_path = os.getcwd()
            elif req_path == "~":
                req_path = str(Path.home())

            p = Path(req_path).expanduser().resolve()
            if not p.is_dir():
                p = p.parent

            entries = []
            # Parent directory link
            parent = str(p.parent)
            if parent != str(p):  # not at root
                entries.append({"name": "..", "type": "dir", "path": parent})

            try:
                items = sorted(p.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower()))
            except PermissionError:
                items = []

            for item in items:
                if item.name.startswith("."):
                    continue  # skip hidden files
                if item.is_dir():
                    entries.append({"name": item.name, "type": "dir", "path": str(item)})
                elif item.suffix == ".h5ad":
                    size_mb = round(item.stat().st_size / 1e6, 1)  # decimal MB (matches macOS)
                    cache_path = item.with_suffix(".fleek_cache.npz")
                    csc_path = item.with_suffix(".csc_cache.npz")
                    annot_m = item.with_suffix(".annot_markers.json")
                    annot_l = item.with_suffix(".annot_llm.json")
                    caches = []
                    if cache_path.exists(): caches.append("fleek")
                    if csc_path.exists(): caches.append("csc")
                    if annot_m.exists(): caches.append("markers")
                    if annot_l.exists(): caches.append("claude")
                    annot_lin = item.with_suffix(".annot_lineage.json")
                    if annot_lin.exists(): caches.append("lineage")
                    is_loaded = (LOADED_PATH and str(item.resolve()) == str(Path(LOADED_PATH).resolve()))
                    entries.append({
                        "name": item.name, "type": "h5ad", "path": str(item),
                        "size_mb": size_mb, "caches": caches, "loaded": is_loaded
                    })

            result = {"path": str(p), "entries": entries}
            data = json.dumps(result).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(500)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_complete(self, partial):
        """Return path completions for the browse path input."""
        try:
            if not partial:
                data = json.dumps({"completions": []}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)
                return

            p = Path(partial).expanduser()
            # If partial ends with /, complete children of that directory
            if partial.endswith("/") or partial.endswith(os.sep):
                parent = p
                prefix = ""
            else:
                parent = p.parent
                prefix = p.name  # case-sensitive on Unix

            completions = []
            if parent.is_dir():
                try:
                    for item in sorted(parent.iterdir(), key=lambda x: x.name):
                        if item.name.startswith("."):
                            continue
                        if prefix and not item.name.startswith(prefix):
                            continue
                        if item.is_dir():
                            completions.append(str(item) + "/")
                        elif item.suffix == ".h5ad":
                            completions.append(str(item))
                        if len(completions) >= 15:
                            break
                except PermissionError:
                    pass

            data = json.dumps({"completions": completions}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            err = json.dumps({"completions": []}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_progress(self):
        p = dict(PROGRESS)
        p["csc_building"] = CSC_BUILDING
        p["csc_cached"] = CSC_CACHED
        p["csc_time"] = CSC_TIME
        p["hvg_ready"] = len(HVG_NAMES) > 0
        p["pacmap_computing"] = PACMAP_COMPUTING
        p["pacmap_ready"] = PACMAP_2D is not None
        p["annot"] = dict(ANNOT_STATUS)
        if len(HVG_NAMES) > 0:
            p["hvg_genes"] = HVG_NAMES[:200]
            if HVG_VAR:
                p["hvg_var"] = {g: round(HVG_VAR[g], 4) for g in HVG_NAMES[:200] if g in HVG_VAR}
        data = json.dumps(p).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(data)

    def _serve_gene_var(self):
        """Return all genes sorted by variability (descending).
        Each gene entry: [name, score, cutoff_flag_bitmask]
        cutoff_flag bitmask: 1=not expressed, 2=low mean, 4=high mean, 8=low fraction
        0 = passes all traditional cutoffs.
        """
        try:
            if ALL_GENE_VAR is None or GENE_NAMES_LIST is None:
                result = {"genes": [], "metric": "dispersions_norm"}
            else:
                order = np.argsort(ALL_GENE_VAR)[::-1]
                has_flags = GENE_CUTOFF_FLAG is not None
                genes = []
                for i in order:
                    v = float(ALL_GENE_VAR[i])
                    flag = int(GENE_CUTOFF_FLAG[i]) if has_flags else 0
                    genes.append([str(GENE_NAMES_LIST[i]), round(v, 4), flag])
                result = {"genes": genes, "metric": ALL_GENE_VAR_METRIC, "total": len(GENE_NAMES_LIST)}
            data = json.dumps(_json_safe(result), separators=(",", ":")).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            err = json.dumps({"genes": [], "error": str(e)}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_cluster_genes(self, params):
        """Return full gene rankings for a specific cluster (or all clusters) from the shared one-vs-rest DEG.
        Each gene has: name, log2fc, padj, cohens_d.  Triggers DEG computation if not cached."""
        try:
            if ADATA is None or CLUSTER_NAMES is None:
                raise ValueError("No dataset loaded")
            # Ensure shared DEG is computed (may trigger computation)
            if _DEG_CACHE.get("full_rankings") is None:
                print(f"  Cluster genes: computing shared DEG for full rankings...")
                _get_shared_deg(top_n=50, test="wilcoxon")
            rankings = _DEG_CACHE.get("full_rankings", {})
            from collections import Counter
            counts = Counter(CLUSTER_IDS.tolist())
            id_param = params.get("id", ["0"])[0]
            if id_param == "all":
                # Return all clusters in one response
                all_results = {}
                for cluster_id in range(len(CLUSTER_NAMES)):
                    cname = CLUSTER_NAMES[cluster_id]
                    cdata = rankings.get(cname)
                    n_cells = counts.get(cluster_id, 0)
                    if cdata is None:
                        all_results[str(cluster_id)] = {"genes": [], "log2fc": [], "padj": [], "cohens_d": [],
                                                         "cluster": cname, "cluster_id": cluster_id, "n_cells": n_cells,
                                                         "n_genes": 0, "small": True}
                    else:
                        all_results[str(cluster_id)] = {
                            "genes": cdata["genes"], "log2fc": cdata["log2fc"],
                            "padj": cdata["padj"], "cohens_d": cdata["cohens_d"],
                            "cluster": cname, "cluster_id": cluster_id,
                            "n_cells": n_cells, "n_genes": len(cdata["genes"]),
                        }
                result = {"all": True, "clusters": all_results}
            else:
                cluster_id = int(id_param)
                if cluster_id < 0 or cluster_id >= len(CLUSTER_NAMES):
                    raise ValueError(f"Invalid cluster id {cluster_id}")
                cname = CLUSTER_NAMES[cluster_id]
                cdata = rankings.get(cname)
                n_cells = counts.get(cluster_id, 0)
                if cdata is None:
                    result = {"genes": [], "log2fc": [], "padj": [], "cohens_d": [],
                              "cluster": cname, "cluster_id": cluster_id, "n_cells": n_cells,
                              "n_genes": 0, "small": True}
                else:
                    result = {
                        "genes": cdata["genes"], "log2fc": cdata["log2fc"],
                        "padj": cdata["padj"], "cohens_d": cdata["cohens_d"],
                        "cluster": cname, "cluster_id": cluster_id,
                        "n_cells": n_cells, "n_genes": len(cdata["genes"]),
                    }
            data = json.dumps(_json_safe(result), separators=(",", ":")).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_search(self, query):
        if ADATA is None:
            data = json.dumps([]).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        q = query.strip().upper()
        if not q:
            matches = []
        else:
            # Exact matches first, then prefix, then substring
            exact = [g for g in GENE_NAMES_LIST if g.upper() == q]
            prefix = [g for g in GENE_NAMES_LIST if g.upper().startswith(q) and g.upper() != q]
            substr = [g for g in GENE_NAMES_LIST if q in g.upper() and not g.upper().startswith(q)]
            matches = (exact + prefix + substr)[:50]
        data = json.dumps(matches).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_go_search(self, query):
        """Search GO term names, return matches with dataset gene counts."""
        if GO_DB is None:
            data = json.dumps({"error": "No GO database loaded", "results": []}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        q = query.strip().lower()
        if not q or len(q) < 2:
            data = json.dumps({"results": []}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        terms = GO_DB.get("terms", {})
        resolver = GO_SYNONYM_MAP or {}
        results = []
        for go_id, t in terms.items():
            name = t.get("name", "")
            if q in name.lower():
                genes = t.get("genes", [])
                # Count how many of this term's genes resolve to dataset genes
                n_dataset = sum(1 for g in genes if g in resolver)
                if n_dataset > 0:
                    results.append({
                        "id": go_id,
                        "name": name,
                        "ns": t.get("ns", ""),
                        "n_total": len(genes),
                        "n_dataset": n_dataset
                    })
        # Sort by dataset gene count descending, then by name
        results.sort(key=lambda x: (-x["n_dataset"], x["name"]))
        results = results[:30]
        data = json.dumps({"results": results}, separators=(",", ":")).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_go_genes(self, go_id):
        """Return genes for a GO term, intersected with the current dataset."""
        if GO_DB is None:
            data = json.dumps({"error": "No GO database loaded"}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        terms = GO_DB.get("terms", {})
        t = terms.get(go_id)
        if not t:
            data = json.dumps({"error": f"Unknown GO term: {go_id}"}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        resolver = GO_SYNONYM_MAP or {}
        all_genes = t.get("genes", [])
        # Resolve GO gene names to dataset gene names
        matched = []
        for g in all_genes:
            if g in resolver:
                matched.append(resolver[g])
        # Deduplicate while preserving order
        seen = set()
        unique = []
        for g in matched:
            if g not in seen:
                seen.add(g)
                unique.append(g)
        result = {
            "id": go_id,
            "name": t.get("name", ""),
            "ns": t.get("ns", ""),
            "genes": sorted(unique),
            "n_total": len(all_genes),
            "n_dataset": len(unique)
        }
        data = json.dumps(result, separators=(",", ":")).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_go_gene_terms(self, gene):
        """Return GO terms for a specific gene."""
        if GO_DB is None:
            data = json.dumps([]).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        gene_to_terms = GO_DB.get("gene_to_terms", {})
        synonyms = GO_DB.get("synonyms", {})
        terms = GO_DB.get("terms", {})
        # Try exact, then uppercase, then synonym
        tids = gene_to_terms.get(gene)
        if not tids:
            tids = gene_to_terms.get(gene.upper())
        if not tids:
            canonical = synonyms.get(gene) or synonyms.get(gene.upper())
            if canonical:
                tids = gene_to_terms.get(canonical)
        if not tids:
            tids = []
        result = []
        for tid in tids:
            t = terms.get(tid)
            if t:
                result.append({"id": tid, "name": t.get("name", ""), "ns": t.get("ns", "")})
        # Sort by namespace then name
        result.sort(key=lambda x: (x["ns"], x["name"]))
        data = json.dumps(result, separators=(",", ":")).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_pseudobulk_deg(self, req):
        """Run pseudo-bulk DEG between two selection groups."""
        try:
            group_a_cells = req.get("group_a_cells")
            group_b_cells = req.get("group_b_cells")
            replicate_col = req.get("replicate_col")
            group_a_name = req.get("group_a_name", "Group A")
            group_b_name = req.get("group_b_name", "Group B")
            method_pref = req.get("method_pref", "auto")

            if not group_a_cells or not group_b_cells or not replicate_col:
                raise ValueError("Missing required fields: group_a_cells, group_b_cells, replicate_col")

            t0 = time.time()
            result = run_pseudobulk_deg(group_a_cells, group_b_cells, replicate_col,
                                         group_a_name=group_a_name, group_b_name=group_b_name,
                                         method_pref=method_pref)
            result["elapsed"] = round(time.time() - t0, 1)

            data = json.dumps(_json_safe(result), separators=(",", ":")).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except BrokenPipeError:
            pass
        except AbortError:
            global ABORT_REQUESTED
            ABORT_REQUESTED = False  # clear so the next run isn't pre-aborted
            print("  Pseudo-bulk aborted by user")
            try:
                err = json.dumps({"error": "aborted", "aborted": True}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(err)))
                self.end_headers()
                self.wfile.write(err)
            except BrokenPipeError:
                pass
        except Exception as e:
            import traceback
            traceback.print_exc()
            try:
                err = json.dumps({"error": str(e)}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(err)))
                self.end_headers()
                self.wfile.write(err)
            except BrokenPipeError:
                pass

    def _serve_deg(self, req):
        """Run DEG between two cell groups."""
        if ADATA is None:
            err = json.dumps({"error": "No dataset loaded"}).encode("utf-8")
            self.send_response(503)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)
            return
        try:
            result = run_deg(
                group_a=req["group_a"],
                group_b=req["group_b"],
                test=req.get("test", "wilcoxon"),
            )
            data = json.dumps(_json_safe(result)).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            import traceback
            traceback.print_exc()
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(500)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_export(self, req):
        """Export cell subset as h5ad for download."""
        try:
            indices = req["indices"]
            name = req.get("name", "subset")
            fpath, fname = export_h5ad_subset(indices, name)
            fsize = os.path.getsize(fpath)
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Disposition", f'attachment; filename="{fname}"')
            self.send_header("Content-Length", str(fsize))
            self.end_headers()
            # Stream in chunks — macOS sendall fails on buffers > ~2GB
            CHUNK = 4 * 1024 * 1024  # 4MB
            with open(fpath, "rb") as f:
                while True:
                    chunk = f.read(CHUNK)
                    if not chunk:
                        break
                    self.wfile.write(chunk)
            os.remove(fpath)
        except Exception as e:
            import traceback
            traceback.print_exc()
            # Clean up temp file if it exists
            try:
                if 'fpath' in locals() and os.path.exists(fpath):
                    os.remove(fpath)
            except OSError:
                pass
            err = json.dumps({"error": str(e)}).encode("utf-8")
            try:
                self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(err)))
                self.end_headers()
                self.wfile.write(err)
            except Exception:
                pass  # headers may already be sent

    def _serve_load(self, req):
        """Load an h5ad file from a local path (no upload needed)."""
        try:
            fpath = req.get("path", "")
            print(f"  Load request: path={fpath!r}")
            if not fpath or not os.path.exists(fpath):
                raise FileNotFoundError(f"File not found: {fpath}")
            if not fpath.endswith(".h5ad"):
                raise ValueError("File must be .h5ad")
            fast = req.get("fast_umap", True)
            fast_sub = int(req.get("fast_umap_n", 50000))
            backed = req.get("backed", "off")
            _native_2d = req.get("native_2d", False)
            DEG_SETTINGS["fast"] = req.get("fast_deg", True)
            DEG_SETTINGS["max_cells"] = int(req.get("fast_deg_n", 50000))
            PACMAP_SETTINGS["fast"] = req.get("fast_pacmap", True)

            _progress("Loading from local path...", 1)
            def _bg():
                try:
                    with PROCESSING_LOCK:
                        load_and_prepare(fpath, max_cells=0, n_dims_list=[2, 3],
                                         fast_umap=fast, fast_umap_subsample=fast_sub,
                                         backed=backed, native_2d=_native_2d)
                except AbortError:
                    print("  Load aborted by user")
                    _reset_all()
                    PROGRESS["status"] = "aborted"
                    PROGRESS["message"] = "Load cancelled"
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    PROGRESS["status"] = "error"
                    PROGRESS["message"] = str(e)
            threading.Thread(target=_bg, daemon=True).start()

            result = {"ok": True, "path": fpath}
            data = json.dumps(result).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            import traceback
            traceback.print_exc()
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(500)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_annotate(self, req):
        """Run cell type annotation against marker database."""
        try:
            top_n = req.get("top_n", 50)
            test = req.get("test", "wilcoxon")
            result = annotate_clusters(top_n=top_n, test=test)
            data = json.dumps(_json_safe(result)).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except BrokenPipeError:
            pass
        except Exception as e:
            import traceback
            traceback.print_exc()
            try:
                err = json.dumps({"error": str(e)}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(err)))
                self.end_headers()
                self.wfile.write(err)
            except BrokenPipeError:
                pass

    def _serve_annotate_llm(self, req):
        """Run cell type annotation using Claude LLM."""
        try:
            if not ANTHROPIC_API_KEY:
                raise ValueError("No API key configured. Set ANTHROPIC_API_KEY env var or add it to ~/.fleek.env.")
            top_n = req.get("top_n", 30)
            tissue = req.get("tissue", "")
            result = annotate_clusters_llm(top_n=top_n, tissue_hint=tissue)
            data = json.dumps(_json_safe(result)).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except BrokenPipeError:
            pass  # Client disconnected (e.g. abort button) — silently ignore
        except Exception as e:
            import traceback
            traceback.print_exc()
            try:
                err = json.dumps({"error": str(e)}).encode("utf-8")
                self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(err)))
                self.end_headers()
                self.wfile.write(err)
            except BrokenPipeError:
                pass

    def _serve_clear_cache(self, req):
        """Delete a specific cache file."""
        try:
            if not LOADED_PATH:
                raise ValueError("No dataset loaded")
            cache_type = req.get("type", "")
            suffix_map = {
                "fleek": ".fleek_cache.npz",
                "csc": ".csc_cache.npz",
                "markers": ".annot_markers.json",
                "claude": ".annot_llm.json",
                "lineage": ".annot_lineage.json",
                "cluster_genes": ".fleek_cluster_genes.json",
            }
            if cache_type not in suffix_map:
                raise ValueError(f"Unknown cache type: {cache_type}")
            # Check both dataset and server cache locations, delete both
            dp = _dataset_cache_path(suffix_map[cache_type])
            sp = _server_cache_path(suffix_map[cache_type])
            deleted = False
            for cp in [dp, sp]:
                if cp and cp.exists():
                    try:
                        cp.unlink()
                        deleted = True
                        print(f"  Deleted cache: {cp}")
                    except (OSError, PermissionError):
                        print(f"  Cannot delete read-only cache: {cp}")
            # Update LOAD_SETTINGS
            settings_map = {
                "fleek": ["umap_cached", "pca_cached", "pacmap_cached", "leiden_cached"],
                "csc": ["csc_cached"],
                "markers": ["annot_markers_cached"],
                "claude": ["annot_llm_cached"],
                "lineage": ["annot_lineage_cached"],
            }
            for key in settings_map.get(cache_type, []):
                LOAD_SETTINGS[key] = False
            # Also wipe the in-memory shared-DEG cache when its disk copy is
            # cleared, so the next /api/cluster-genes forces a fresh compute.
            if cache_type == "cluster_genes":
                _DEG_CACHE["key"] = None
                _DEG_CACHE["cluster_genes"] = None
                _DEG_CACHE["small_clusters"] = None
                _DEG_CACHE["full_rankings"] = None
            result = {"ok": True, "deleted": deleted, "type": cache_type}
            data = json.dumps(result).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_cache_settings(self, req):
        """Get or set cache location mode."""
        try:
            global CACHE_MODE, CACHE_DIR
            new_mode = req.get("mode")
            new_dir = req.get("dir")
            if new_mode and new_mode in ("dataset", "server"):
                CACHE_MODE = new_mode
                LOAD_SETTINGS["cache_mode"] = CACHE_MODE
                print(f"  Cache mode set to: {CACHE_MODE}")
            if new_dir:
                CACHE_DIR = Path(new_dir)
                _init_cache_dir()
                LOAD_SETTINGS["cache_dir"] = str(CACHE_DIR)
                print(f"  Cache dir set to: {CACHE_DIR}")
            result = {
                "ok": True,
                "cache_mode": CACHE_MODE,
                "cache_dir": str(CACHE_DIR) if CACHE_DIR else "",
            }
            data = json.dumps(result).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_key_status(self):
        """Return whether an API key is currently loaded (never the key itself)."""
        result = json.dumps({"has_key": bool(ANTHROPIC_API_KEY)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_set_key(self, req):
        """Save a new Anthropic API key to ~/.fleek.env and reload it in memory."""
        global ANTHROPIC_API_KEY
        try:
            key = req.get("key", "").strip()
            if not key:
                raise ValueError("Key is empty.")
            if not key.startswith("sk-ant-"):
                raise ValueError("Key doesn't look like an Anthropic key (expected sk-ant-...).")
            env_path = Path.home() / ".fleek.env"
            # Preserve any other lines in the file
            lines = []
            if env_path.exists():
                for line in env_path.read_text().splitlines():
                    if not line.strip().startswith("ANTHROPIC_API_KEY="):
                        lines.append(line)
            lines.append(f"ANTHROPIC_API_KEY={key}")
            env_path.write_text("\n".join(lines) + "\n")
            env_path.chmod(0o600)
            ANTHROPIC_API_KEY = key
            print("  Claude API key updated via settings UI.")
            result = json.dumps({"ok": True}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_clear_key(self):
        """Remove the Anthropic API key from memory and ~/.fleek.env."""
        global ANTHROPIC_API_KEY
        try:
            ANTHROPIC_API_KEY = ""
            env_path = Path.home() / ".fleek.env"
            if env_path.exists():
                lines = [l for l in env_path.read_text().splitlines()
                         if not l.strip().startswith("ANTHROPIC_API_KEY=")]
                if lines:
                    env_path.write_text("\n".join(lines) + "\n")
                else:
                    env_path.unlink()
            print("  Claude API key cleared.")
            result = json.dumps({"ok": True}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_test_key(self, req):
        """Verify an Anthropic API key by making a minimal /v1/messages call.
        If `key` is provided in the body, test that specific value without
        persisting. Otherwise test the key currently loaded in memory.
        Returns the raw API response body on failure so users can see the
        actual reason (401, 429, 403 workspace-scoped, DNS, proxy, etc.)."""
        import urllib.request, urllib.error
        try:
            candidate = (req.get("key") or "").strip() if isinstance(req, dict) else ""
            key = candidate or ANTHROPIC_API_KEY
            if not key:
                raise ValueError("No key to test. Paste a key into the field, or save one first.")
            if not key.startswith("sk-ant-"):
                raise ValueError("Key doesn't look like an Anthropic key (expected sk-ant-...).")
            body = json.dumps({
                "model": "claude-haiku-4-5-20251001",
                "max_tokens": 8,
                "messages": [{"role": "user", "content": "ping"}],
            }).encode("utf-8")
            rq = urllib.request.Request(
                "https://api.anthropic.com/v1/messages",
                data=body,
                headers={
                    "Content-Type": "application/json",
                    "x-api-key": key,
                    "anthropic-version": "2023-06-01",
                },
                method="POST",
            )
            try:
                with urllib.request.urlopen(rq, timeout=15) as resp:
                    data = json.loads(resp.read().decode("utf-8", errors="replace"))
                result = json.dumps({"ok": True, "model": data.get("model", "")}).encode("utf-8")
            except urllib.error.HTTPError as e:
                raw = e.read().decode("utf-8", errors="replace")
                # Surface the verbatim API error so the user can see the real cause.
                result = json.dumps({"ok": False, "error": f"HTTP {e.code}: {raw}"}).encode("utf-8")
            except urllib.error.URLError as e:
                result = json.dumps({"ok": False, "error": f"Network error: {e.reason}"}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    # ── Session management ──

    def _session_suffix(self, name):
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)
        return f".fleek_session_{safe}.json"

    def _serve_sessions(self):
        """List all saved sessions for the current dataset."""
        try:
            if not LOADED_PATH:
                raise ValueError("No dataset loaded")
            _init_cache_dir()
            stem = Path(LOADED_PATH).stem
            seen = {}
            # Check dataset dir
            dp = Path(LOADED_PATH).parent
            if dp.exists():
                for f in dp.glob(f"{stem}.fleek_session_*.json"):
                    name = f.stem.replace(f"{stem}.fleek_session_", "")
                    seen[name] = {"name": name, "path": str(f),
                                  "mtime": f.stat().st_mtime, "auto": name == "_auto"}
            # Check server cache dir
            if CACHE_DIR and CACHE_DIR.exists():
                for f in CACHE_DIR.glob(f"{stem}.fleek_session_*.json"):
                    name = f.stem.replace(f"{stem}.fleek_session_", "")
                    if name not in seen:
                        seen[name] = {"name": name, "path": str(f),
                                      "mtime": f.stat().st_mtime, "auto": name == "_auto"}
            sessions = sorted(seen.values(), key=lambda s: (not s["auto"], s["name"]))
            result = json.dumps({"sessions": sessions}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"sessions": [], "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_session_save(self, req):
        """Save a named session."""
        try:
            if not LOADED_PATH:
                raise ValueError("No dataset loaded")
            name = req.get("name", "").strip()
            if not name:
                raise ValueError("Session name is empty")
            state = req.get("state", {})
            state["_meta"] = {"n_cells": N_CELLS, "n_clusters": len(CLUSTER_NAMES or []),
                              "dataset": Path(LOADED_PATH).name,
                              "saved": time.time()}
            suffix = self._session_suffix(name)
            wp = _cache_write_path(suffix)
            with open(wp, "w") as f:
                json.dump(state, f, separators=(",", ":"))
            print(f"  Session '{name}' saved to {wp}")
            result = json.dumps({"ok": True, "name": name}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_session_load(self, req):
        """Load a named session."""
        try:
            if not LOADED_PATH:
                raise ValueError("No dataset loaded")
            name = req.get("name", "").strip()
            if not name:
                raise ValueError("Session name is empty")
            suffix = self._session_suffix(name)
            path, _ = _cache_read_path(suffix)
            if not path:
                raise ValueError(f"Session '{name}' not found")
            with open(path) as f:
                state = json.load(f)
            meta = state.get("_meta", {})
            if meta.get("n_cells") and meta["n_cells"] != N_CELLS:
                raise ValueError(f"Session has {meta['n_cells']} cells but dataset has {N_CELLS}")
            result = json.dumps({"ok": True, "state": state}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_session_delete(self, req):
        """Delete a named session from both cache locations."""
        try:
            if not LOADED_PATH:
                raise ValueError("No dataset loaded")
            name = req.get("name", "").strip()
            if not name:
                raise ValueError("Session name is empty")
            suffix = self._session_suffix(name)
            dp = _dataset_cache_path(suffix)
            sp = _server_cache_path(suffix)
            deleted = False
            for cp in [dp, sp]:
                if cp and cp.exists():
                    cp.unlink()
                    deleted = True
                    print(f"  Session '{name}' deleted: {cp}")
            result = json.dumps({"ok": True, "deleted": deleted}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_lineage(self, req):
        """Construct a developmental lineage tree for the current cell types using Claude."""
        try:
            if not ANTHROPIC_API_KEY:
                raise ValueError("No API key configured. Set ANTHROPIC_API_KEY env var or add it to ~/.fleek.env.")
            if ADATA is None or CLUSTER_NAMES is None:
                raise ValueError("No dataset loaded")

            # Check cache (with parent stem fallback for subsets)
            _parent_stem = ADATA.uns.get("fleek_parent_stem") if ADATA is not None else None
            cache_path, _cw = _cache_read_path(".annot_lineage.json")
            if not cache_path and _parent_stem and CACHE_DIR:
                _pp = CACHE_DIR / f"{_parent_stem}.annot_lineage.json"
                if _pp.exists():
                    cache_path = _pp
            if cache_path:
                try:
                    with open(cache_path) as f:
                        cached = json.load(f)
                    _cn = cached.get("cluster_names")
                    _nc = cached.get("n_cells")
                    _valid = (_cn == CLUSTER_NAMES and _nc == N_CELLS) or (_cn and set(CLUSTER_NAMES).issubset(set(_cn)))
                    if _valid:
                        print(f"  Lineage loaded from cache: {cache_path.name if hasattr(cache_path, 'name') else cache_path}")
                        result = {"ok": True, "tree": cached["tree"], "cached": True}
                        data = json.dumps(result).encode("utf-8")
                        self.send_response(200)
                        self.send_header("Content-Type", "application/json")
                        self.send_header("Content-Length", str(len(data)))
                        self.end_headers()
                        self.wfile.write(data)
                        return
                except Exception:
                    pass

            # Build prompt with cluster names and cell counts
            cluster_info = []
            from collections import Counter
            counts = Counter(CLUSTER_IDS.tolist())

            # Load Claude annotations for inline context
            # Priority: cache file > parent stem fallback > client-supplied predictions
            llm_preds = {}  # cluster_id -> predicted cell type
            llm_cache_path, _ = _cache_read_path(".annot_llm.json")
            if not llm_cache_path and _parent_stem and CACHE_DIR:
                _pp = CACHE_DIR / f"{_parent_stem}.annot_llm.json"
                if _pp.exists():
                    llm_cache_path = _pp
            if llm_cache_path:
                try:
                    with open(llm_cache_path) as f:
                        llm_data = json.load(f)
                    for r in llm_data.get("results", []):
                        cid = r.get("cluster_id")
                        preds = r.get("predictions", [])
                        if preds and isinstance(preds[0], dict):
                            llm_preds[cid] = preds[0].get("cell_type", "")
                except Exception:
                    pass
            # Fallback: use client-supplied predictions (covers case where cache wasn't saved)
            if not llm_preds:
                client_preds = req.get("llm_predictions", {})
                for k, v in client_preds.items():
                    if v:
                        llm_preds[int(k)] = v

            for i, name in enumerate(CLUSTER_NAMES):
                n_cells = counts.get(i, 0)
                if n_cells == 0:
                    continue  # Skip empty clusters — no cells, no lineage
                # When Claude annotations exist, use inferred name as primary
                # (raw cluster labels like "R" or "RV" can mislead lineage construction)
                if i in llm_preds and llm_preds[i]:
                    entry = f'  Cluster {i}: "{llm_preds[i]}" ({n_cells:,} cells) [original label: "{name}"]'
                else:
                    entry = f'  Cluster {i}: "{name}" ({n_cells:,} cells)'
                cluster_info.append(entry)

            n_in_lineage = len(cluster_info)
            n_skipped = len(CLUSTER_NAMES) - n_in_lineage
            if n_skipped > 0:
                print(f"  Lineage: using {n_in_lineage} non-empty clusters (skipping {n_skipped} empty)")

            llm_context = ""

            prompt = f"""You are a developmental biology expert. Given these cell types from a single-cell RNA-seq dataset, construct a developmental lineage tree.

Cell types in this dataset:
{chr(10).join(cluster_info)}
{llm_context}

Instructions:
1. Build a tree showing the developmental hierarchy of these cell types
2. The root should be the last common developmental ancestor of all cell types present
3. Include intermediate progenitor stages that connect the cell types, even if they are not present in the dataset — mark these as "inferred"
4. Each node that corresponds to an actual cluster in the dataset should reference it by cluster_id
5. A single dataset cluster may map to one node. If two clusters are the same developmental identity (e.g. just numbered differently), still give each its own leaf
6. Keep inferred intermediates minimal — only add nodes needed to connect the lineage, not an exhaustive textbook tree
7. Use standard developmental biology nomenclature for inferred nodes

Return ONLY a JSON object with this structure (no markdown, no explanation):
{{
  "root": {{
    "name": "Ancestor Name",
    "present": false,
    "cluster_ids": [],
    "children": [
      {{
        "name": "Progenitor",
        "present": false,
        "cluster_ids": [],
        "children": [
          {{
            "name": "Cell Type A",
            "present": true,
            "cluster_ids": [0],
            "children": []
          }}
        ]
      }}
    ]
  }}
}}

"present" means the cell type exists as a cluster in the dataset (true) or is an inferred developmental intermediate (false).
"cluster_ids" is an array of cluster indices that map to this node (empty for inferred nodes).
Every non-empty cluster_id listed above must appear exactly once in the tree. Empty clusters (0 cells) have been excluded."""

            print(f"  Requesting lineage tree from Claude ({n_in_lineage} cell types, {len(llm_preds)} with inferred names)...")
            import urllib.request, urllib.error
            request_body = json.dumps({
                "model": "claude-sonnet-4-20250514",
                "max_tokens": 16384,
                "messages": [{"role": "user", "content": prompt}]
            }).encode("utf-8")

            api_req = urllib.request.Request(
                "https://api.anthropic.com/v1/messages",
                data=request_body,
                headers={
                    "Content-Type": "application/json",
                    "x-api-key": ANTHROPIC_API_KEY,
                    "anthropic-version": "2023-06-01"
                },
                method="POST"
            )

            try:
                with urllib.request.urlopen(api_req, timeout=300) as resp:
                    resp_data = json.loads(resp.read().decode("utf-8"))
            except urllib.error.HTTPError as e:
                body = e.read().decode("utf-8", errors="replace")
                raise ValueError(f"Claude API error {e.code}: {body[:200]}")

            usage = resp_data.get("usage", {})
            stop = resp_data.get("stop_reason", "")
            print(f"  Lineage response: {usage.get('input_tokens', 0)} in, {usage.get('output_tokens', 0)} out tokens (stop: {stop})")

            if stop == "max_tokens":
                raise ValueError(f"Lineage response truncated at {usage.get('output_tokens', 0)} tokens — tree too large for {n_in_lineage} cell types. Try loading a smaller subset.")

            # Parse response
            text = ""
            for block in resp_data.get("content", []):
                if block.get("type") == "text":
                    text += block["text"]

            text = text.strip()
            if text.startswith("```"):
                text = text.split("\n", 1)[1] if "\n" in text else text[3:]
            if text.endswith("```"):
                text = text[:-3].strip()

            tree = json.loads(text)
            # Accept either {root: ...} or direct node
            if "root" in tree:
                tree = tree["root"]

            # Validate: every cluster_id should appear
            def _collect_ids(node):
                ids = set(node.get("cluster_ids", []))
                for ch in node.get("children", []):
                    ids |= _collect_ids(ch)
                return ids
            found_ids = _collect_ids(tree)
            expected_ids = {i for i, name in enumerate(CLUSTER_NAMES) if counts.get(i, 0) > 0}
            missing = expected_ids - found_ids
            if missing:
                print(f"  Warning: lineage tree missing cluster_ids: {missing}")

            # Cache
            wp = _cache_write_path(".annot_lineage.json")
            if wp:
                try:
                    cache_obj = {
                        "cluster_names": CLUSTER_NAMES,
                        "n_cells": N_CELLS,
                        "tree": tree
                    }
                    with open(wp, "w") as f:
                        json.dump(cache_obj, f)
                    print(f"  Lineage cached to {wp.name}")
                except Exception as e:
                    print(f"  Lineage cache write failed: {e}")

            result = {"ok": True, "tree": tree, "cached": False}
            data = json.dumps(result).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        except BrokenPipeError:
            pass
        except Exception as e:
            import traceback
            traceback.print_exc()
            try:
                err = json.dumps({"error": str(e)}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(err)))
                self.end_headers()
                self.wfile.write(err)
            except BrokenPipeError:
                pass

    def _serve_compute_embedding(self, req):
        """Compute an embedding on demand (e.g. PaCMAP) and return coordinates."""
        global PACMAP_2D, PACMAP_3D
        method = req.get("method", "")
        try:
            if ADATA is None:
                raise ValueError("No dataset loaded")

            if method == "pacmap":
                # Check if already computed
                if PACMAP_2D is not None and PACMAP_3D is not None:
                    result = {"ok": True, "cached": True}
                    data = json.dumps(result).encode("utf-8")
                    self.send_response(200)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Content-Length", str(len(data)))
                    self.end_headers()
                    self.wfile.write(data)
                    return

                pca_input = ADATA.obsm.get("X_pca")
                if pca_input is None:
                    raise ValueError("No PCA coordinates available")
                pca_input = pca_input[:, :30].astype(np.float32)
                n_total = pca_input.shape[0]
                PACMAP_SUB = PACMAP_SETTINGS.get("max_cells", 50000) if PACMAP_SETTINGS.get("fast", True) else 0

                # Run PaCMAP in a subprocess to avoid numba/threading segfaults
                import tempfile, subprocess
                with tempfile.TemporaryDirectory() as tmpdir:
                    input_path = os.path.join(tmpdir, "pca_in.npy")
                    out2_path = os.path.join(tmpdir, "pacmap_2d.npy")
                    out3_path = os.path.join(tmpdir, "pacmap_3d.npy")
                    np.save(input_path, pca_input)

                    if PACMAP_SUB > 0 and n_total > PACMAP_SUB:
                        sub_idx = _stratified_subsample(CLUSTER_IDS, PACMAP_SUB,
                                                        np.random.default_rng(42))
                        sub_path = os.path.join(tmpdir, "sub_idx.npy")
                        np.save(sub_path, np.array(sub_idx, dtype=np.int32))
                        script = f"""
import numpy as np, pacmap, time
from scipy.spatial import cKDTree
X = np.load("{input_path}")
sub_idx = np.load("{sub_path}")
X_sub = X[sub_idx]
nn = min(10, max(2, len(sub_idx) // 50))
t0 = time.time()
print(f"  PaCMAP 2D ({{len(sub_idx):,}} / {{X.shape[0]:,}} cells)...")
Y2_sub = pacmap.PaCMAP(n_components=2, n_neighbors=nn, verbose=False).fit_transform(X_sub).astype(np.float32)
print(f"  PaCMAP 3D...")
Y3_sub = pacmap.PaCMAP(n_components=3, n_neighbors=nn, verbose=False).fit_transform(X_sub).astype(np.float32)
print(f"  Projecting remaining cells (k=5 weighted)...")
tree = cKDTree(X_sub)
rest_mask = np.ones(X.shape[0], dtype=bool)
rest_mask[sub_idx] = False
rest_idx = np.where(rest_mask)[0]
dists, knn_idx = tree.query(X[rest_idx], k=5)
weights = 1.0 / (dists + 1e-10)
weights /= weights.sum(axis=1, keepdims=True)
Y2 = np.empty((X.shape[0], 2), dtype=np.float32)
Y3 = np.empty((X.shape[0], 3), dtype=np.float32)
Y2[sub_idx] = Y2_sub; Y3[sub_idx] = Y3_sub
Y2[rest_idx] = np.einsum("nk,nkd->nd", weights, Y2_sub[knn_idx]).astype(np.float32)
Y3[rest_idx] = np.einsum("nk,nkd->nd", weights, Y3_sub[knn_idx]).astype(np.float32)
np.save("{out2_path}", Y2); np.save("{out3_path}", Y3)
print(f"  PaCMAP done ({{time.time()-t0:.1f}}s)")
"""
                    else:
                        script = f"""
import numpy as np, pacmap, time
X = np.load("{input_path}")
nn = min(10, max(2, X.shape[0] // 50))
t0 = time.time()
print(f"  PaCMAP 2D ({{X.shape[0]:,}} cells, {{nn}} neighbors)...")
Y2 = pacmap.PaCMAP(n_components=2, n_neighbors=nn, verbose=False).fit_transform(X).astype(np.float32)
np.save("{out2_path}", Y2)
print(f"  PaCMAP 3D...")
Y3 = pacmap.PaCMAP(n_components=3, n_neighbors=nn, verbose=False).fit_transform(X).astype(np.float32)
np.save("{out3_path}", Y3)
print(f"  PaCMAP done ({{time.time()-t0:.1f}}s)")
"""
                    t0 = time.time()
                    print(f"  Computing PaCMAP in subprocess ({N_CELLS:,} cells)...")
                    proc = subprocess.run(
                        [sys.executable, "-c", script],
                        capture_output=True, text=True
                    )
                    # Print subprocess output
                    if proc.stdout:
                        for line in proc.stdout.strip().split("\n"):
                            print(line)
                    if proc.returncode != 0:
                        stderr = proc.stderr.strip()[-500:] if proc.stderr else "unknown error"
                        raise RuntimeError(f"PaCMAP subprocess failed (exit {proc.returncode}): {stderr}")

                    PACMAP_2D = np.load(out2_path)
                    PACMAP_3D = np.load(out3_path)
                    elapsed = round(time.time() - t0, 1)
                    print(f"  PaCMAP total: {elapsed}s")

                # Save to cache with method-specific keys
                if LOADED_PATH:
                    _pm_sfx = "_quick" if PACMAP_SETTINGS.get("fast", True) else "_full"
                    cache_path = _cache_write_path(".fleek_cache.npz")
                    if cache_path.exists():
                        try:
                            cached = dict(np.load(cache_read, allow_pickle=False))
                            cached[f"pacmap_2d{_pm_sfx}"] = PACMAP_2D
                            cached[f"pacmap_3d{_pm_sfx}"] = PACMAP_3D
                            np.savez(_cache_write_path(".fleek_cache.npz"), **cached)
                            print(f"  PaCMAP cached ({_pm_sfx.strip('_')})")
                        except Exception:
                            pass

                result = {"ok": True, "elapsed": elapsed}
            else:
                raise ValueError(f"Unknown embedding method: {method}")

            data = json.dumps(result).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            import traceback
            traceback.print_exc()
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(500)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def _serve_abort(self, req):
        """Abort in-progress loading."""
        global ABORT_REQUESTED
        ABORT_REQUESTED = True
        print("  Abort requested by client")
        result = {"ok": True}
        data = json.dumps(result).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_unload(self, req):
        """Unload current dataset and free all memory."""
        global ABORT_REQUESTED
        ABORT_REQUESTED = True  # Stop any background processing
        print("  Unload requested by client")
        _reset_all()
        PROGRESS["status"] = "idle"
        PROGRESS["message"] = ""
        PROGRESS["pct"] = 0
        import gc; gc.collect()
        result = {"ok": True}
        data = json.dumps(result).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_heartbeat(self, req):
        """Client heartbeat — resets the auto-unload timer."""
        global AUTO_UNLOAD, AUTO_UNLOAD_TIMEOUT, _LAST_HEARTBEAT, _HEARTBEAT_TIMER
        _LAST_HEARTBEAT = time.time()
        print(f"  ♥ heartbeat (auto_unload={req.get('auto_unload') if req else '?'}, timeout={req.get('timeout') if req else '?'})")
        # Cancel any pending unload timer
        if _HEARTBEAT_TIMER:
            _HEARTBEAT_TIMER.cancel()
            _HEARTBEAT_TIMER = None
        # Parse settings updates from client
        if req:
            if "auto_unload" in req:
                AUTO_UNLOAD = bool(req["auto_unload"])
            if "timeout" in req:
                AUTO_UNLOAD_TIMEOUT = max(0.5, float(req["timeout"]))
        # Schedule next unload check
        if AUTO_UNLOAD and ADATA is not None:
            def _check_unload():
                global _HEARTBEAT_TIMER
                _HEARTBEAT_TIMER = None
                elapsed = time.time() - _LAST_HEARTBEAT
                if elapsed >= AUTO_UNLOAD_TIMEOUT and AUTO_UNLOAD and ADATA is not None:
                    print(f"  Auto-unloading dataset (no heartbeat for {elapsed:.1f}s)")
                    _reset_all()
                    PROGRESS["status"] = "idle"
                    PROGRESS["message"] = "Dataset auto-unloaded"
                    import gc; gc.collect()
            _HEARTBEAT_TIMER = threading.Timer(AUTO_UNLOAD_TIMEOUT + 0.5, _check_unload)
            _HEARTBEAT_TIMER.daemon = True
            _HEARTBEAT_TIMER.start()
        data = json.dumps({"ok": True, "auto_unload": AUTO_UNLOAD, "timeout": AUTO_UNLOAD_TIMEOUT}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_upload(self):
        """Receive an h5ad file upload, save to disk, then process in background."""
        import tempfile
        try:
            length = int(self.headers.get("Content-Length", 0))
            fname = self.headers.get("X-Filename", "upload.h5ad")
            if not fname.endswith(".h5ad"):
                raise ValueError("File must be .h5ad")
            _progress("Receiving upload...", 1)
            fpath = os.path.join(tempfile.gettempdir(), fname)
            received = 0
            with open(fpath, "wb") as f:
                while received < length:
                    chunk = self.rfile.read(min(65536, length - received))
                    if not chunk:
                        break
                    f.write(chunk)
                    received += len(chunk)
            _progress("Upload saved, starting processing...", 3)
            fast = self.headers.get("X-Fast-UMAP", "0") == "1"
            fast_sub = int(self.headers.get("X-Fast-UMAP-N", "50000"))
            backed = self.headers.get("X-Backed", "off")  # "auto", "on", "off"
            _native_2d = self.headers.get("X-Native-2D", "0") == "1"
            # QuickDEG settings
            DEG_SETTINGS["fast"] = self.headers.get("X-Fast-DEG", "1") == "1"
            DEG_SETTINGS["max_cells"] = int(self.headers.get("X-Fast-DEG-N", "50000"))
            PACMAP_SETTINGS["fast"] = self.headers.get("X-Fast-PaCMAP", "1") == "1"

            # Start processing in background thread
            def _bg():
                try:
                    with PROCESSING_LOCK:
                        load_and_prepare(fpath, max_cells=0, n_dims_list=[2, 3],
                                         fast_umap=fast, fast_umap_subsample=fast_sub,
                                         backed=backed, native_2d=_native_2d)
                except AbortError:
                    print("  Load aborted by user")
                    _reset_all()
                    PROGRESS["status"] = "aborted"
                    PROGRESS["message"] = "Load cancelled"
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    PROGRESS["status"] = "error"
                    PROGRESS["message"] = str(e)

            threading.Thread(target=_bg, daemon=True).start()

            # Respond immediately — client will poll /api/progress
            result = {"processing": True}
            data = json.dumps(result).encode("utf-8")
            self.send_response(202)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        except Exception as e:
            import traceback
            traceback.print_exc()
            PROGRESS["status"] = "error"
            PROGRESS["message"] = str(e)
            err = json.dumps({"error": str(e)}).encode("utf-8")
            self.send_response(500)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(err)))
            self.end_headers()
            self.wfile.write(err)

    def log_message(self, format, *args):
        # Quieter logging
        if "/api/" in str(args[0]):
            return
        super().log_message(format, *args)


def main():
    parser = argparse.ArgumentParser(description="RNA-FLEEK: serve h5ad directly to browser")
    parser.add_argument("h5ad", type=str, nargs="?", default=None,
                        help="Path to h5ad file (optional — can upload via browser)")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--max-cells", type=int, default=0,
                        help="Subsample to this many cells (0 = all)")
    parser.add_argument("--dims", type=str, default="2,3",
                        help="UMAP dims to compute, comma-separated (default: 2,3)")
    parser.add_argument("--host", type=str, default="127.0.0.1")
    parser.add_argument("--fast-umap", action="store_true",
                        help="Use subsample-and-project for large datasets (>50k cells)")
    parser.add_argument("--fast-umap-n", type=int, default=50000,
                        help="Number of cells for UMAP subsample (default: 50000)")
    parser.add_argument("--backed", type=str, default="off", choices=["auto", "on", "off"],
                        help="Backed mode: auto (detect), on (force disk), off (force RAM, default)")
    parser.add_argument("--cache-mode", type=str, default="dataset", choices=["dataset", "server"],
                        help="Where to store caches: 'dataset' (alongside h5ad, default) or 'server' (~/.fleek_cache)")
    parser.add_argument("--cache-dir", type=str, default=None,
                        help="Custom server cache directory (default: ~/.fleek_cache)")
    parser.add_argument("--no-browser", action="store_true",
                        help="Don't auto-open browser (useful for remote/SSH sessions)")
    args = parser.parse_args()

    global CACHE_MODE, CACHE_DIR
    CACHE_MODE = args.cache_mode
    if args.cache_dir:
        CACHE_DIR = Path(args.cache_dir)
    _init_cache_dir()

    dims = [int(d) for d in args.dims.split(",")]

    if args.h5ad:
        load_and_prepare(args.h5ad, max_cells=args.max_cells, n_dims_list=dims,
                         fast_umap=args.fast_umap, fast_umap_subsample=args.fast_umap_n,
                         backed=args.backed)

    server = type('ThreadingHTTPServer', (ThreadingMixIn, HTTPServer), {'daemon_threads': True})((args.host, args.port), FleekHandler)
    url = f"http://{args.host}:{args.port}"
    print(f"\n{'='*50}")
    print(f"  RNA-FLEEK running at: {url}")
    if ADATA is not None:
        print(f"  {N_CELLS:,} cells · {len(CLUSTER_NAMES)} types · {len(GENE_NAMES_LIST):,} genes")
    else:
        print(f"  No dataset loaded — upload via browser")
    key_status = "loaded" if ANTHROPIC_API_KEY else "not set (Claude annotation disabled)"
    print(f"  Claude API key: {key_status}")
    print(f"{'='*50}\n")

    # Try to open browser (skip if --no-browser, or if running over SSH)
    is_ssh = os.environ.get("SSH_CLIENT") or os.environ.get("SSH_CONNECTION") or os.environ.get("SSH_TTY")
    if not args.no_browser and not is_ssh:
        try:
            import webbrowser
            webbrowser.open(url)
        except:
            pass
    elif is_ssh:
        print(f"  (SSH session detected — open {url} in your local browser)")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()