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

# Kept in sync with rna_fleek/__init__.py, pyproject.toml, fleek.html FLEEK_VERSION.
FLEEK_VERSION = "0.4.16"

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
X_CSC = None  # CSC copy of adata.X for fast column access
# Secondary CSC index used when the gene_display routing picks a matrix other
# than adata.X (e.g. log_normalized stored in a layer, or raw.X when X holds
# something else). Built lazily on load. Without it, gene lookups on that
# route would fall into CSR column-slicing — O(nnz) per call, catastrophic
# on >1M-cell datasets.
DISPLAY_CSC = None
DISPLAY_CSC_PATH = None
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
# Pathway / gene-set databases for over-representation analysis. Built on
# dataset load from the existing GO_DB plus any GMT files in
# data/gmt/{organism}/. Each entry: {label, source, n_pathways, terms:
# {id: {name, description, url, genes:[...]}}}. Memory: typical ~50 MB
# for GO + Reactome, mostly the gene-name strings.
PATHWAY_DBS = {}            # {db_id: {label, source, terms: {...}}}
PATHWAY_ORGANISM = None     # organism name PATHWAY_DBS was built for
_PATHWAY_ORA_CACHE = {}     # {(genes_hash, db_id, n_bg): result_dict} — cleared on dataset change
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
AUTO_UNLOAD = False   # If True, auto-quit (unload dataset + exit process) when no client heartbeat. Off by default — users opt in via the Settings toggle.
AUTO_UNLOAD_TIMEOUT = 20.0  # seconds without heartbeat before auto-quit
_LAST_HEARTBEAT = 0.0  # time.time() of last heartbeat
_HEARTBEAT_TIMER = None  # threading.Timer for delayed unload

# ── Cache location management ──
# CACHE_MODE: "dataset" = store caches alongside .h5ad file (default)
#             "server"  = store caches in a server-local directory
# Falls back to server dir when dataset directory is not writable.
# Reads check both locations (dataset dir first, then server dir).
CACHE_MODE = "dataset"
CACHE_DIR = None  # Path object, set at startup (defaults to ~/.fleek_cache)

# ── Remote-server mode ────────────────────────────────────────────────
# When REMOTE_MODE is True the server is being run as a multi-user lab
# tool behind nginx + supervisor. USER_DATA_DIR is the per-user (or
# test-user) writable root; SHARED_DATA_DIR is the canonical, RO root
# the admin pre-populates with shared atlases + pre-baked caches.
# USER_QUOTA_MB caps how much each user may write into their dir
# (settings + sessions + uploads + exports + caches combined).
# All four are set at startup from CLI flags; downstream code reads
# them via the _user_path() helper rather than touching the globals
# directly. In the default single-user mode none of this is set and
# the existing /tmp paths apply.
REMOTE_MODE = False
USER_DATA_DIR = None      # Path; required when REMOTE_MODE is True
SHARED_DATA_DIR = None    # Path; required when REMOTE_MODE is True
USER_QUOTA_MB = 100

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

def _user_path(subdir):
    """Return <USER_DATA_DIR>/<subdir>/, ensuring the dir exists.

    The canonical place for per-user writable artifacts in remote-mode:
    sessions, uploads, exports, generated caches. In single-user mode
    USER_DATA_DIR is None and this returns None — callers fall back to
    their pre-existing default paths (typically tempfile.gettempdir()).
    """
    if USER_DATA_DIR is None:
        return None
    p = USER_DATA_DIR / subdir
    p.mkdir(parents=True, exist_ok=True)
    return p

def _user_dir_used_bytes():
    """Recursively sum the byte size of USER_DATA_DIR. Used by the
    per-user write-quota check before accepting an upload / export /
    session-save. Cheap (fits-in-RAM 5-figure file counts) for the
    quotas we'll set; revisit if quotas grow into the GB range and
    walking the tree becomes a measurable cost."""
    if USER_DATA_DIR is None or not USER_DATA_DIR.exists():
        return 0
    total = 0
    for root, _dirs, files in os.walk(USER_DATA_DIR):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total

def _is_under(path, root):
    """True iff `path` resolves to something inside `root`. Both args
    are resolved before comparison, so symlink shenanigans (../, link
    farms, etc.) don't escape. Used by every endpoint that takes a
    user-controllable path in remote mode — browse, load, clear-cache.
    Returns False on any resolution error so the default in remote mode
    is "deny"."""
    try:
        p = Path(path).expanduser().resolve()
        r = Path(root).resolve()
        return r in p.parents or p == r
    except Exception:
        return False

def _path_in_allowed_roots(path):
    """True iff `path` is under USER_DATA_DIR or SHARED_DATA_DIR. In
    non-remote mode every path is allowed (preserves single-user
    behaviour). In remote mode this is the gatekeeper for /api/browse,
    /api/load, /api/clear-cache, etc."""
    if not REMOTE_MODE:
        return True
    return _is_under(path, USER_DATA_DIR) or _is_under(path, SHARED_DATA_DIR)

def _path_in_shared(path):
    """True iff `path` resolves under SHARED_DATA_DIR. Used to refuse
    delete operations on shared caches — the admin owns those."""
    if SHARED_DATA_DIR is None:
        return False
    return _is_under(path, SHARED_DATA_DIR)

def _quota_check_or_raise(extra_bytes=0):
    """In remote mode, refuse the operation if the new write would push
    USER_DATA_DIR over USER_QUOTA_MB. extra_bytes is the size of the
    pending write. In non-remote mode this is a no-op so the existing
    single-user flow is unaffected."""
    if not REMOTE_MODE:
        return
    cap = USER_QUOTA_MB * 1024 * 1024
    if _user_dir_used_bytes() + extra_bytes > cap:
        raise RuntimeError(
            f"Quota exceeded: this user's directory is at the {USER_QUOTA_MB} MB cap. "
            f"Delete some sessions / uploads / exports and try again."
        )

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
    In 'server' mode, always uses server dir.
    Returns None if no writable path is resolvable (e.g. LOADED_PATH unset)."""
    _init_cache_dir()
    if CACHE_MODE == "dataset":
        dp = _dataset_cache_path(suffix)
        if dp and _is_writable(dp):
            return dp
        print(f"  Cache: dataset dir not writable, using server cache dir")
    return _server_cache_path(suffix)

def _safe_savez(suffix, **kwargs):
    """Resolve a write path and np.savez into it — but never raise. Returning
    None for suffix is always possible in edge cases (unset LOADED_PATH, wiped
    CACHE_DIR). Previously these paths flowed straight into np.savez, producing
    "expected str, bytes or os.PathLike object, not NoneType" and failing the
    whole load. Now cache-save failures just log and continue."""
    try:
        wp = _cache_write_path(suffix)
        if wp is None:
            print(f"  Cache write skipped ({suffix}): no writable path (LOADED_PATH or CACHE_DIR unset)")
            return None
        np.savez(wp, **kwargs)
        return wp
    except Exception as _e:
        print(f"  Cache write failed ({suffix}): {_e} — continuing")
        return None

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


def _fetch_ncbi_gene_summary(entrez_id, timeout=8):
    """Pull the authoritative gene summary straight from NCBI E-utilities
    (the same text rendered on the public https://www.ncbi.nlm.nih.gov/gene/
    page). mygene.info's `summary` field is NCBI RefSeq-curated and often
    empty for model organisms — the Alliance of Genome Resources / MGI / RGD
    / ZFIN-contributed summaries that NCBI shows live in the `esummary`
    response instead. Low call volume (cached per-gene per-organism), so we
    run unauthenticated; E-utilities rate-limits at 3 req/s without an API
    key which we'll never approach here. Returns a stripped string or "". """
    if not entrez_id:
        return ""
    try:
        import urllib.request, urllib.parse
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
               "?db=gene&retmode=json&id=" + urllib.parse.quote(str(entrez_id)))
        req = urllib.request.Request(url, headers={"User-Agent": "RNA-FLEEK"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = json.loads(r.read().decode("utf-8"))
        rec = (data.get("result") or {}).get(str(entrez_id)) or {}
        s = (rec.get("summary") or "").strip()
        return s
    except Exception as e:
        print(f"  NCBI esummary fallback failed for {entrez_id}: {e}")
        return ""


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

    # Build pathway databases (GO subdivided by namespace + any GMT files
    # found under data/gmt/{organism}/). Wrapping happens after GO load so
    # GO sub-databases reuse the already-parsed term table.
    _build_pathway_dbs(organism)


def _build_pathway_dbs(organism):
    """Populate the pathway-database registry.

    Three data sources are merged into one registry, keyed by short id:
      • GO BP / MF / CC — derived from GO_DB by splitting on namespace.
      • GMT files in data/gmt/{organism}/*.gmt — one database per file,
        id = file stem (e.g. data/gmt/human/hallmark.gmt -> "hallmark").

    Each registry entry: {label, source, n_pathways, terms: {id: {name,
    description, url, genes: [...]}}}. Stored as plain dicts so JSON
    serialization for /api/pathway-databases is trivial.
    """
    global PATHWAY_DBS, PATHWAY_ORGANISM, _PATHWAY_ORA_CACHE
    PATHWAY_DBS = {}
    _PATHWAY_ORA_CACHE = {}
    PATHWAY_ORGANISM = organism
    # ---- GO sub-databases ----
    if GO_DB and GO_DB.get("terms"):
        ns_buckets = {"BP": {}, "MF": {}, "CC": {}}
        ns_label = {"BP": "Biological Process", "MF": "Molecular Function", "CC": "Cellular Component"}
        for tid, term in GO_DB["terms"].items():
            ns = term.get("ns")
            if ns not in ns_buckets:
                continue
            ns_buckets[ns][tid] = {
                "name": term.get("name", tid),
                "description": "",
                "url": "https://www.ebi.ac.uk/QuickGO/term/" + tid,
                "genes": term.get("genes", []),
            }
        for ns_code, terms in ns_buckets.items():
            if not terms:
                continue
            db_id = "go_" + ns_code.lower()
            PATHWAY_DBS[db_id] = {
                "label": "GO " + ns_label[ns_code],
                "source": "QuickGO",
                "n_pathways": len(terms),
                "terms": terms,
            }
    # ---- GMT files in data/gmt/{organism}/ ----
    script_dir = Path(__file__).parent
    gmt_roots = [
        script_dir / "data" / "gmt" / organism,
        _BUNDLE_DIR / "rna_fleek" / "data" / "gmt" / organism,
        _BUNDLE_DIR / "data" / "gmt" / organism,
        Path("data") / "gmt" / organism,
    ]
    for gmt_root in gmt_roots:
        if not Path(gmt_root).is_dir():
            continue
        for gmt_path in sorted(Path(gmt_root).glob("*.gmt")):
            try:
                terms = {}
                with open(gmt_path) as f:
                    for line in f:
                        parts = line.rstrip("\n").rstrip("\r").split("\t")
                        if len(parts) < 3:
                            continue
                        pid = parts[0].strip()
                        if not pid:
                            continue
                        descr_raw = parts[1].strip()
                        url = ""
                        description = descr_raw
                        # Heuristic: if column 2 looks like a URL, treat it as such; else as description.
                        if descr_raw.startswith("http://") or descr_raw.startswith("https://"):
                            url = descr_raw
                            description = ""
                        genes = [g.strip() for g in parts[2:] if g.strip()]
                        if not genes:
                            continue
                        terms[pid] = {
                            "name": pid.replace("_", " ").title() if "_" in pid else pid,
                            "description": description,
                            "url": url,
                            "genes": genes,
                        }
                if not terms:
                    continue
                db_id = gmt_path.stem.lower()
                # Pretty label: Hallmark, Reactome, etc. (file stem with first letter cap)
                label = gmt_path.stem.replace("_", " ").title()
                PATHWAY_DBS[db_id] = {
                    "label": label,
                    "source": "GMT (" + gmt_path.name + ")",
                    "n_pathways": len(terms),
                    "terms": terms,
                }
                print(f"  Loaded pathway database: {db_id} ({len(terms)} sets) from {gmt_path}")
            except Exception as e:
                print(f"  Warning: failed to load GMT {gmt_path}: {e}")
    if PATHWAY_DBS:
        print(f"  Pathway databases ready: {len(PATHWAY_DBS)} total ({', '.join(sorted(PATHWAY_DBS.keys()))})")
    else:
        print(f"  No pathway databases available for {organism}")


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


def _is_mmap_backed(arr):
    """True if `arr` or any of its base chain is an np.memmap. Used to
    detect whether a sparse matrix's data array was loaded via mmap vs
    read into RAM — our fast path builds arrays with np.frombuffer over a
    memmap, so the array itself isn't a memmap but its .base is."""
    if arr is None:
        return False
    cur = arr
    seen = 0
    while cur is not None and seen < 8:
        if isinstance(cur, np.memmap):
            return True
        cur = getattr(cur, "base", None)
        seen += 1
    return False


def _load_csc_mmap(npz_path):
    """Open an uncompressed scipy-sparse .npz (CSC/CSR) and memory-map each
    of its three array members in place — no buffered read through Python's
    zipfile stream, no full pull into RAM. Returns a scipy.sparse matrix
    whose data/indices/indptr are np.memmap views of the file on disk;
    pages fault in only when a column is actually touched. Near-instant
    startup replaces the previous full-RAM load, which for a 1.6M-cell
    matrix was several GB of I/O routed through the zipfile reader.

    The scipy .npz format is a standard zip archive containing four .npy
    members: data, indices, indptr, shape, and format. When saved with
    `compressed=False` the members are stored with ZIP_STORED, which
    means the bytes inside each .npy are at a known offset in the outer
    file — perfect for np.memmap. Compressed caches (which FLEEK doesn't
    create but someone might have copied in) fall back to the slow
    sp.load_npz path via the caller.

    Returns None on any parse error or compressed-entry encounter — the
    caller then does scipy's regular load_npz.
    """
    import zipfile
    import struct
    from numpy.lib import format as _npfmt
    import scipy.sparse as _sp
    try:
        entries = {}
        data_offsets = {}
        with zipfile.ZipFile(str(npz_path), "r") as zf:
            for info in zf.infolist():
                if info.compress_type != zipfile.ZIP_STORED:
                    return None  # compressed — needs decompression, skip mmap
                entries[info.filename] = info
        # Resolve each entry's absolute file offset by reading its local
        # file header (30 bytes + variable-length filename + extra fields).
        # ZipInfo uses __slots__ so we can't attach attributes to it; keep a
        # parallel dict instead.
        with open(str(npz_path), "rb") as fp:
            for fname, info in entries.items():
                fp.seek(info.header_offset)
                lh = fp.read(30)
                if len(lh) < 30 or lh[:4] != b"PK\x03\x04":
                    return None
                fname_len = struct.unpack("<H", lh[26:28])[0]
                extra_len = struct.unpack("<H", lh[28:30])[0]
                data_offsets[fname] = info.header_offset + 30 + fname_len + extra_len
        import mmap as _mmap
        _page = _mmap.ALLOCATIONGRANULARITY
        def _mmap_npy(name):
            if name not in entries:
                return None
            with open(str(npz_path), "rb") as fp:
                fp.seek(data_offsets[name])
                version = _npfmt.read_magic(fp)
                if version == (1, 0):
                    shape, fortran, dtype = _npfmt.read_array_header_1_0(fp)
                elif version == (2, 0):
                    shape, fortran, dtype = _npfmt.read_array_header_2_0(fp)
                else:
                    return None
                array_offset = fp.tell()
            if fortran:
                return None
            # np.memmap(offset=N) requires N to be a multiple of the OS
            # allocation granularity (4K on Linux/mac, 64K on Windows). .npy
            # array payloads inside a STORED .npz are only 64-byte aligned
            # relative to their zip member, so in practice the absolute
            # offset is rarely page-aligned. Work around this by mmapping
            # a page-aligned raw-byte span that covers (padding + payload),
            # then slicing out the payload as a typed numpy view — still
            # zero-copy, backed by the same mmap pages.
            nbytes = int(np.dtype(dtype).itemsize) * int(np.prod(shape))
            if nbytes == 0:
                return np.empty(shape, dtype=dtype)
            aligned = (array_offset // _page) * _page
            pad = array_offset - aligned
            try:
                raw = np.memmap(str(npz_path), dtype=np.uint8, mode="r",
                                offset=aligned, shape=(pad + nbytes,))
            except Exception:
                return None
            return np.frombuffer(raw, dtype=dtype, count=int(np.prod(shape)),
                                 offset=pad).reshape(shape)
        data = _mmap_npy("data.npy")
        indices = _mmap_npy("indices.npy")
        indptr = _mmap_npy("indptr.npy")
        if data is None or indices is None or indptr is None:
            return None
        # shape + format are tiny; read normally via np.load.
        with np.load(str(npz_path)) as loaded:
            shape_arr = tuple(int(x) for x in loaded["shape"].tolist())
            fmt_raw = loaded["format"].item()
        fmt = fmt_raw.decode() if isinstance(fmt_raw, (bytes, bytearray)) else str(fmt_raw)
        if fmt == "csc":
            return _sp.csc_matrix((data, indices, indptr), shape=shape_arr)
        elif fmt == "csr":
            return _sp.csr_matrix((data, indices, indptr), shape=shape_arr)
        return None
    except Exception:
        return None


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
    """Load Anthropic API key from the most user-specific source available.

    Order of preference (first match wins):
      1. <data-dir>/api_key — only in remote mode, with USER_DATA_DIR set.
         This is where /api/set-key writes for an authenticated per-user
         setting; takes priority over the env var so a key the user
         updated via the UI overrides whatever the supervisor passed at
         spawn time (which came from the supervisor DB at signup).
      2. ANTHROPIC_API_KEY env var — set by the supervisor at spawn for
         remote mode, or by the user's shell for single-user mode.
      3. ~/.fleek.env — single-user fallback only; ignored in remote
         mode since the home-dir file is shared across all per-user
         fleek subprocesses on the same VM.

    Called once at module import (with no globals set, so just env +
    ~/.fleek.env are checked) and again inside main() after CLI args
    are parsed and USER_DATA_DIR is known — the latter is what picks
    up the data-dir/api_key file in remote mode.
    """
    # 1. Per-user file (remote mode, after main() has set USER_DATA_DIR)
    if REMOTE_MODE and USER_DATA_DIR is not None:
        kf = USER_DATA_DIR / "api_key"
        if kf.exists():
            try:
                key = kf.read_text().strip()
                if key:
                    return key
            except Exception:
                pass
    # 2. Env var
    key = os.environ.get("ANTHROPIC_API_KEY", "")
    if key:
        return key
    # 3. ~/.fleek.env (single-user mode only)
    if not REMOTE_MODE:
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
        n_cells_deg = adata_deg.n_obs

        # Pre-compute per-cluster sum and sum-of-squares ONCE via one sparse
        # matmul with a cluster indicator matrix. Previously the loop did
        # X[~mask].power(2).mean() per cluster — O(n_clusters × nnz(X)) of
        # wasted work that cost several minutes on large datasets.
        valid_cnames = [c for c in CLUSTER_NAMES if c not in small_clusters]
        cname_to_cidx = {c: i for i, c in enumerate(valid_cnames)}
        sum_per_cluster = None
        sumsq_per_cluster = None
        n_per_cluster = None
        total_sum = None
        total_sumsq = None
        gene_to_varidx = {str(g): j for j, g in enumerate(adata_deg.var_names)}
        if valid_cnames:
            row_ids = np.array([cname_to_cidx.get(c, -1) for c in cluster_labels], dtype=np.int64)
            keep = row_ids >= 0
            rows = row_ids[keep]
            cols = np.where(keep)[0]
            C = sp.csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)),
                              shape=(len(valid_cnames), n_cells_deg))
            CX = C @ X
            sum_per_cluster = CX.toarray() if sp.issparse(CX) else np.asarray(CX)
            X_sq = X.power(2) if sp.issparse(X) else (np.asarray(X) ** 2)
            CX2 = C @ X_sq
            sumsq_per_cluster = CX2.toarray() if sp.issparse(CX2) else np.asarray(CX2)
            n_per_cluster = np.asarray(C.sum(axis=1)).ravel().astype(np.int64)
            total_sum = sum_per_cluster.sum(axis=0)
            total_sumsq = sumsq_per_cluster.sum(axis=0)

        for cname in CLUSTER_NAMES:
            if cname in small_clusters:
                full_rankings[cname] = None
                continue
            try:
                all_names = [str(g) for g in rgg["names"][cname]]
                all_lfc = rgg["logfoldchanges"][cname].astype(np.float64)
                all_padj = rgg["pvals_adj"][cname].astype(np.float64)
                cidx = cname_to_cidx.get(cname, -1)
                if cidx >= 0 and sum_per_cluster is not None:
                    n1 = int(n_per_cluster[cidx])
                    n2 = n_cells_deg - n1
                else:
                    n1 = n2 = 0
                if n1 >= 2 and n2 >= 2:
                    s1 = sum_per_cluster[cidx]
                    s2 = total_sum - s1
                    q1 = sumsq_per_cluster[cidx]
                    q2 = total_sumsq - q1
                    m1 = s1 / n1
                    m2 = s2 / n2
                    # sample variance: (sumsq - n*mean^2) / (n-1)
                    v1 = np.maximum(q1 - n1 * m1 * m1, 0.0) / (n1 - 1)
                    v2 = np.maximum(q2 - n2 * m2 * m2, 0.0) / (n2 - 1)
                    pooled = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
                    d_all = (m1 - m2) / np.maximum(pooled, 1e-9)
                    # Reorder Cohen's d to match rgg ranking order
                    idx_arr = np.fromiter((gene_to_varidx.get(g, -1) for g in all_names),
                                          dtype=np.int64, count=len(all_names))
                    safe_idx = np.where(idx_arr >= 0, idx_arr, 0)
                    all_d = np.where(idx_arr >= 0, d_all[safe_idx], 0.0)
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


# ── Claude model catalogue ──────────────────────────────────────────
# The set of models the FLEEK Settings UI offers users for the Claude-
# backed annotation + lineage features. Cost figures are approximate
# per-1M-token Anthropic list prices (input/output) as of late 2025;
# the dropdown shows a "$ per typical run" estimate derived from these
# plus the typical token-count of an annotation pass (~8k input + ~3k
# output for 50 clusters × 30 marker genes). Users see this BEFORE
# clicking Annotate so a $5 run isn't surprising.
#
# IDs map to Anthropic's `model` field; the API will accept the short
# alias (e.g. "claude-opus-4-7") and return the dated model in the
# response's `model` field, which fleek then displays in the tooltip.
CLAUDE_MODELS = {
    "haiku": {
        "id": "claude-haiku-4-5-20251001",
        "label": "Haiku 4.5",
        "cost_hint": "fast, cheap (~$0.01/run)",
    },
    "sonnet": {
        "id": "claude-sonnet-4-6",
        "label": "Sonnet 4.6",
        "cost_hint": "balanced (~$0.10/run, default)",
    },
    "opus": {
        "id": "claude-opus-4-7",
        "label": "Opus 4.7",
        "cost_hint": "most capable, expensive (~$0.50–1/run)",
    },
}
CLAUDE_DEFAULT_MODEL_KEY = "sonnet"


def _resolve_claude_model(req_dict_or_value):
    """Pick a Claude model ID from a request payload or a bare value.

    Accepts either the short key ('haiku' / 'sonnet' / 'opus') or a
    full model ID ('claude-sonnet-4-6'). Unknown values fall back to
    the default. The two-form input is a small ergonomic concession
    so the frontend can send a short key while CLI / scripts can
    pass a literal model ID."""
    if isinstance(req_dict_or_value, dict):
        v = (req_dict_or_value.get("model") or "").strip()
    else:
        v = (req_dict_or_value or "").strip() if isinstance(req_dict_or_value, str) else ""
    if v in CLAUDE_MODELS:
        return CLAUDE_MODELS[v]["id"]
    if v.startswith("claude-"):
        return v
    return CLAUDE_MODELS[CLAUDE_DEFAULT_MODEL_KEY]["id"]


def annotate_clusters_llm(top_n=30, test="wilcoxon", tissue_hint="", model=None):
    """Use Claude API to infer cell types from top marker genes per cluster.
    Batches clusters into groups to avoid token limits.

    `model` (optional): Claude model ID to use (e.g. 'claude-sonnet-4-6'
    or short key 'sonnet' / 'opus' / 'haiku'). Defaults to Sonnet 4.6.

    Returns same format as annotate_clusters for UI compatibility.
    """
    model_id = _resolve_claude_model(model)
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
            "model": model_id,
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
    global DISPLAY_CSC, DISPLAY_CSC_PATH
    DISPLAY_CSC = None; DISPLAY_CSC_PATH = None
    HVG_NAMES = []; HVG_VAR = {}; ALL_GENE_VAR = None; GENE_CUTOFF_FLAG = None; GENE_INDEX = {}
    CSC_BUILDING = False; CSC_CACHED = False; CSC_TIME = 0
    LOADED_PATH = ""; LOAD_SETTINGS.clear()
    COUNTS_MATRIX = None; COUNTS_LABEL = None; COUNTS_INFO = None
    # Pathway databases are organism-bound; clear them so a different
    # dataset doesn't see stale entries while reload is in flight.
    global PATHWAY_DBS, PATHWAY_ORGANISM, _PATHWAY_ORA_CACHE
    PATHWAY_DBS = {}; PATHWAY_ORGANISM = None; _PATHWAY_ORA_CACHE = {}
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

# ── Skipped-layer registry (populated by _read_h5ad_skip_dense_layers) ─
# Maps layer name → {"shape", "dtype", "nbytes_gb"} for any layer the
# read step omitted from the in-memory AnnData because loading it would
# have OOMed. Display gene queries for these layers can still be served
# from a pre-built display-CSC cache on disk; DEG / preprocessing on
# these layers is unsupported in this fallback path. The user-visible
# message (printed by the loader) names every skipped layer.
SKIPPED_LAYERS = {}


def _layer_display_csc_cache_path(layer_name):
    """Resolve the on-disk display-CSC cache file for a given layer name,
    or return None if no such cache exists. Mirrors the suffix-sanitisation
    in the existing display-CSC build path so cached files written by
    earlier runs (suffix derived from `layers['<name>']`) are recognised.
    Used to short-circuit sparse-layer reads at load time when an
    equivalent cache is already on disk."""
    _dpath = f"layers['{layer_name}']"
    _safe = (_dpath
             .replace("/", "_")
             .replace("[", "_")
             .replace("]", "_")
             .replace("'", "")
             .replace('"', ""))
    suffix = f".display_csc_{_safe}.npz"
    cp, _ = _cache_read_path(suffix)
    if cp is not None and cp.exists():
        return cp
    return None


def _read_h5ad_skip_dense_layers(path, dense_skip_gb=5.0, skip_X=False):
    """Read an h5ad while omitting any DENSE layer above `dense_skip_gb`.

    Why this exists: anndata's read_h5ad — even with backed='r' — eagerly
    materializes every entry in /layers into RAM. For datasets with a
    large dense layer (e.g. denoised 664k×20k float32 = 50 GB), this
    OOMs on Linux servers with no swap, even when .X itself is sparse
    and small. macOS hides the problem because of aggressive memory
    compression + dynamic swap; Linux is strict.

    What it does instead:
      1. Walks the h5ad with h5py
      2. For .X / .obs / .var / .obsm / .obsp / .varm / .varp / .uns /
         .raw — reads via anndata.experimental.read_elem (always small
         compared to a multi-GB dense layer)
      3. For /layers — reads each entry except DENSE ones above the
         threshold. Sparse layers (h5py.Group with CSR/CSC encoding)
         are loaded normally; their cost is bounded by nnz.
      4. Returns (adata, skipped_dict) where skipped_dict records
         {name: {shape, dtype, nbytes_gb}} for every omitted layer.

    Callers should record skipped_dict somewhere globally (e.g.
    SKIPPED_LAYERS) so downstream paths can decide whether to fall
    back to a disk-cached display CSC, ignore the layer, or surface
    a "not available" message in the UI.
    """
    import h5py
    # The read_elem helper has shifted homes between anndata versions —
    # try the public location first, then the internal ones.
    try:
        from anndata.experimental import read_elem
    except ImportError:
        try:
            from anndata._io.specs.registry import read_elem
        except ImportError:
            from anndata._io.specs import read_elem

    skipped = {}
    kwargs = {}

    # Helper to size a single h5ad component in GB — a sparse Group's
    # cost is data + indices + indptr; a dense Dataset is shape × dtype.
    # Used both for reporting in-progress reads and for deciding whether
    # to skip a dense layer.
    def _component_gb(elem):
        if isinstance(elem, h5py.Dataset):
            return (int(np.prod(elem.shape)) * elem.dtype.itemsize) / 1e9
        if isinstance(elem, h5py.Group):
            etype = elem.attrs.get("encoding-type", "")
            if etype in ("csr_matrix", "csc_matrix"):
                total = 0
                for sub in ("data", "indices", "indptr"):
                    if sub in elem:
                        total += elem[sub].size * elem[sub].dtype.itemsize
                return total / 1e9
        return 0.0

    with h5py.File(path, "r") as f:
        # Lightweight components — every reasonable h5ad has these well
        # below the cost of the heavy layer. Read them first so a
        # failure here is easy to diagnose vs a layer-read failure.
        # Each iteration also calls _check_abort so the user can cancel
        # mid-read if a particular component is taking forever (cloud
        # disk being slow, sparse arrays of GB-scale to deserialize),
        # and emits a _progress update so the loading UI shows what's
        # actually happening rather than a static "Loading h5ad file..."
        # for minutes at a time.
        for key in ("X", "obs", "var", "obsm", "obsp", "varm", "varp", "uns", "raw"):
            if key not in f:
                continue
            _check_abort()
            # Optionally skip /X — caller will set adata.X from the
            # csc_cache.npz mmap on disk after we return. This is the
            # main cold-load speedup: a 14 GB sparse-X read goes away
            # entirely when the CSC cache is already present.
            if key == "X" and skip_X:
                cost_gb = _component_gb(f[key])
                _progress(f"Skipping /X read (~{cost_gb:.1f} GB) — will mmap from CSC cache", 5)
                continue
            cost_gb = _component_gb(f[key])
            cost_str = f", ~{cost_gb:.1f} GB" if cost_gb >= 0.5 else ""
            _progress(f"Reading /{key}{cost_str}...", 5)
            try:
                kwargs[key] = read_elem(f[key])
            except Exception as e:
                print(f"  Warn: couldn't read /{key}: {e}")

        # Layers — selectively skip dense ones above the threshold.
        # Same abort-check + progress-update pattern so a slow sparse
        # layer (e.g. raw_counts at 17 GB) doesn't look like a hang.
        # Sparse-layer short-circuit: if a display-CSC cache file
        # already exists alongside the dataset for a given layer,
        # mmap it instead of re-reading the layer from h5ad. The
        # cache file IS the same data (just CSC instead of CSR), so
        # the AnnData ends up with the same .layers[name] content as
        # if we'd read it normally — but with no disk read for the
        # layer itself. For the user's adata_final.h5ad this turns
        # the 14 GB denoised_norm read into a near-instant mmap.
        if "layers" in f:
            layers_dict = {}
            for layer_name in f["layers"]:
                _check_abort()
                layer = f["layers"][layer_name]
                if isinstance(layer, h5py.Dataset):
                    nbytes_gb = (int(np.prod(layer.shape)) * layer.dtype.itemsize) / 1e9
                    if nbytes_gb > dense_skip_gb:
                        skipped[layer_name] = {
                            "shape": tuple(layer.shape),
                            "dtype": str(layer.dtype),
                            "nbytes_gb": nbytes_gb,
                        }
                        print(
                            f"  Skipping dense layer '{layer_name}' "
                            f"({tuple(layer.shape)} {layer.dtype}, {nbytes_gb:.1f} GB) — "
                            f"would OOM in-memory load. Display CSC cache "
                            f"(if present) still serves gene queries from this layer."
                        )
                        continue
                # Display-CSC short-circuit (sparse layers only): same
                # short-circuit pattern as the X-CSC one in v0.5.7. If
                # a cached CSC for this layer is on disk, mmap it and
                # use that as the layer value.
                if isinstance(layer, h5py.Group):
                    _disp_cache_path = _layer_display_csc_cache_path(layer_name)
                    if _disp_cache_path is not None:
                        cost_gb = _component_gb(layer)
                        try:
                            import scipy.sparse as _sp_local
                            _disp_mat = _load_csc_mmap(_disp_cache_path)
                            if _disp_mat is None:
                                _disp_mat = _sp_local.load_npz(str(_disp_cache_path))
                            _progress(
                                f"Skipping /layers/{layer_name} read (~{cost_gb:.1f} GB) — using cached display CSC",
                                6,
                            )
                            layers_dict[layer_name] = _disp_mat
                            continue
                        except Exception as _e:
                            # Cache load failed; fall through to a normal read.
                            print(f"  Display CSC mmap for /layers/{layer_name} failed ({_e}); reading layer from h5ad")
                # Sparse Group, or small dense Dataset — load normally.
                cost_gb = _component_gb(layer)
                cost_str = f", ~{cost_gb:.1f} GB" if cost_gb >= 0.5 else ""
                _progress(f"Reading /layers/{layer_name}{cost_str}...", 6)
                try:
                    layers_dict[layer_name] = read_elem(layer)
                except Exception as e:
                    print(f"  Warn: couldn't read /layers/{layer_name}: {e}")
            if layers_dict:
                kwargs["layers"] = layers_dict

    _check_abort()
    _progress("Constructing AnnData object...", 8)
    import anndata as _ad
    adata = _ad.AnnData(**kwargs)
    return adata, skipped


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
        # Inspect the h5ad's /X group via h5py BEFORE loading — h5ad
        # files use HDF5 compression, so a 50-GB dense float32 matrix
        # lands on disk as 5–10 GB. The naïve "file_size_gb > 0.4 *
        # ram" heuristic decides "plenty of room" for that case, then
        # anndata decompresses to 50 GB and OOMs. Computing the actual
        # in-memory cost from shape × dtype (or the sparse arrays' own
        # sizes) is the right comparison.
        in_memory_cost_gb = file_size_gb  # fallback if inspection fails
        try:
            import h5py
            with h5py.File(h5ad_path, "r") as _f:
                _X = _f.get("X")
                if isinstance(_X, h5py.Dataset):
                    # Dense X: shape × dtype itemsize is the full cost.
                    in_memory_cost_gb = (int(np.prod(_X.shape)) * _X.dtype.itemsize) / 1e9
                elif isinstance(_X, h5py.Group):
                    # Sparse X (CSR/CSC): cost is data + indices + indptr.
                    _dat = _X.get("data");  _idx = _X.get("indices");  _ptr = _X.get("indptr")
                    if _dat is not None and _idx is not None and _ptr is not None:
                        in_memory_cost_gb = (
                            _dat.size * _dat.dtype.itemsize
                            + _idx.size * _idx.dtype.itemsize
                            + _ptr.size * _ptr.dtype.itemsize
                        ) / 1e9
        except Exception as _e:
            print(f"  Auto-detect: couldn't inspect /X metadata ({_e}); using file size as proxy")
        # Use backed mode if the in-memory cost would exceed 50% of
        # available RAM. Threshold gives headroom for working memory
        # (CSC index, embeddings, supervisor + other users on shared VM).
        use_backed = in_memory_cost_gb > avail_ram_gb * 0.5
        print(f"  Auto-detect: file {file_size_gb:.1f}GB on disk, "
              f"~{in_memory_cost_gb:.1f}GB in RAM, "
              f"{avail_ram_gb:.0f}GB available → "
              f"{'BACKED mode' if use_backed else 'in-memory'}")
    # backed == "off" → use_backed stays False

    # Reset the skipped-layer registry; any prior dataset's entries
    # are stale once we start loading a new file.
    global SKIPPED_LAYERS
    SKIPPED_LAYERS = {}
    if use_backed:
        # NOTE: backed mode is incompatible with the skip-dense-layers
        # path because anndata's backed AnnData wraps the on-disk h5ad
        # with all layers eagerly loaded. If a user explicitly requests
        # --backed=on with a dataset that has a huge dense layer, they
        # still hit the OOM. Practical mitigation: if the auto-detect
        # picks backed mode but we'd hit a dense layer, fall back to
        # the skip-layers path. (Auto-detect almost never picks backed
        # for the in-memory cost we now compute, so this is rare.)
        _progress(f"Loading in backed mode ({file_size_gb:.1f}GB file, {avail_ram_gb:.0f}GB RAM avail)...", 5)
        adata = sc.read_h5ad(h5ad_path, backed="r")
        BACKED = True
    else:
        # X-CSC short-circuit: if we already have a CSC index cached on
        # disk, the bytes in /X are redundant — the cache *is* the same
        # data, just transposed for fast column access. Skip the /X
        # read entirely and load the cache instead, saving the slowest
        # single read in the pipeline (10+ GB of sparse data on cold
        # disks). After the helper returns adata with X=None we mmap
        # the cache and assign it; the rest of fleek treats adata.X
        # interchangeably whether it's CSR or CSC.
        _csc_pre, _ = _cache_read_path(".csc_cache.npz")
        skip_X = _csc_pre is not None and _csc_pre.exists()
        # Use the layer-skipping reader. For datasets with no dense
        # layers above the threshold, behaviour is identical to
        # sc.read_h5ad. For datasets with a huge dense layer (denoised
        # outputs, etc.), the offending layer is skipped and the
        # registry is populated so downstream code can route around
        # the gap (e.g. via the display-CSC cache for that layer).
        adata, SKIPPED_LAYERS = _read_h5ad_skip_dense_layers(
            h5ad_path, dense_skip_gb=5.0, skip_X=skip_X,
        )
        BACKED = False
        if skip_X:
            # Load the cached CSC into adata.X. Use mmap when possible
            # so the bytes don't actually fault into RAM until a gene
            # column is accessed; falls back to a regular load_npz on
            # any error. If the cache shape doesn't match the dataset
            # (e.g. you re-uploaded a different version of the h5ad
            # over an old cache), we re-read /X from the h5ad as the
            # safe fallback rather than serving wrong data.
            _progress("Loading X from CSC cache (mmap)...", 9)
            try:
                _x_cached = _load_csc_mmap(_csc_pre)
                if _x_cached is None:
                    _x_cached = sp.load_npz(str(_csc_pre))
                if _x_cached.shape == (adata.n_obs, adata.n_vars):
                    adata.X = _x_cached
                    _mmap_note = " (mmap)" if _is_mmap_backed(getattr(_x_cached, "data", None)) else ""
                    print(f"  X mmapped from CSC cache{_mmap_note} (skipped /X read)")
                else:
                    print(f"  X CSC cache shape mismatch ({_x_cached.shape} vs {(adata.n_obs, adata.n_vars)}); re-reading /X from h5ad")
                    _x_cached = None
                    raise ValueError("shape mismatch")
            except Exception as _e:
                # Fallback: read /X from h5ad after all
                _progress("Re-reading /X from h5ad (cache short-circuit failed)...", 9)
                try:
                    from anndata.experimental import read_elem as _read_elem_fb
                except ImportError:
                    try:
                        from anndata._io.specs.registry import read_elem as _read_elem_fb
                    except ImportError:
                        from anndata._io.specs import read_elem as _read_elem_fb
                import h5py as _h5py
                with _h5py.File(h5ad_path, "r") as _f:
                    if "X" in _f:
                        adata.X = _read_elem_fb(_f["X"])
        if SKIPPED_LAYERS:
            _progress(
                f"Loaded with {len(SKIPPED_LAYERS)} dense layer(s) skipped "
                f"to avoid OOM: {list(SKIPPED_LAYERS.keys())}",
                10,
            )
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
        # Cache key suffix marks how the UMAP was computed (full vs
        # subsampled-and-projected); both are valid coords for display.
        # Check the matching suffix first, fall back to the OTHER suffix
        # so a cache built in full mode is honoured when loading in fast
        # mode and vice versa. Without this fallback fleek would load
        # UMAP from adata.obsm['X_umap'] instead, mark need_save=True,
        # and re-write the entire ~115 MB cache file even though the
        # data was already on disk under a slightly different key.
        cache_key = f"umap_{nd}d{umap_suffix}"
        legacy_key = f"umap_{nd}d"
        other_suffix = "_full" if umap_suffix == "_quick" else "_quick"
        other_key = f"umap_{nd}d{other_suffix}"
        if cache_key in cached:
            _progress(f"Loaded {nd}D UMAP from cache ({'quick' if use_fast else 'full'})", 25)
            if nd == 2:
                UMAP_2D = cached[cache_key].astype(np.float32)
            else:
                UMAP_3D = cached[cache_key].astype(np.float32)
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "cache_quick" if use_fast else "cache_full"
        elif other_key in cached:
            _progress(f"Loaded {nd}D UMAP from cache ({other_suffix.lstrip('_')} variant)", 25)
            if nd == 2:
                UMAP_2D = cached[other_key].astype(np.float32)
            else:
                UMAP_3D = cached[other_key].astype(np.float32)
            EMBEDDING_SOURCES[f"umap_{nd}d"] = "cache_" + other_suffix.lstrip("_")
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
        _safe_savez(".fleek_cache.npz", **cached)

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
                    _safe_savez(".fleek_cache.npz", **cached)
                except Exception:
                    pass
                print(f"  PCA computed in background ({pca_arr.shape[1]} components)")
            except Exception as e:
                print(f"  Background PCA failed: {e}")
        threading.Thread(target=_bg_pca, daemon=True).start()

    # ── PaCMAP embedding (background subprocess to avoid numba/threading issues) ──
    # Same cross-suffix fallback pattern as UMAP — a cache from full mode
    # is valid coords in fast mode and vice versa, so don't recompute or
    # re-save when the data is already on disk under the other suffix.
    _pm_suffix = "_quick" if PACMAP_SETTINGS.get("fast", True) else "_full"
    _pm_other_suffix = "_full" if _pm_suffix == "_quick" else "_quick"
    _pm_key2 = f"pacmap_2d{_pm_suffix}"
    _pm_key3 = f"pacmap_3d{_pm_suffix}"
    _pm_other2 = f"pacmap_2d{_pm_other_suffix}"
    _pm_other3 = f"pacmap_3d{_pm_other_suffix}"
    if _pm_key2 in cached and _pm_key3 in cached:
        PACMAP_2D = cached[_pm_key2].astype(np.float32)
        PACMAP_3D = cached[_pm_key3].astype(np.float32)
        PACMAP_COMPUTING = False
        print(f"  PaCMAP loaded from cache ({'quick' if _pm_suffix == '_quick' else 'full'})")
        src = "cache_quick" if _pm_suffix == "_quick" else "cache_full"
        EMBEDDING_SOURCES["pacmap_2d"] = src
        EMBEDDING_SOURCES["pacmap_3d"] = src
    elif _pm_other2 in cached and _pm_other3 in cached:
        PACMAP_2D = cached[_pm_other2].astype(np.float32)
        PACMAP_3D = cached[_pm_other3].astype(np.float32)
        PACMAP_COMPUTING = False
        print(f"  PaCMAP loaded from cache ({_pm_other_suffix.lstrip('_')} variant)")
        EMBEDDING_SOURCES["pacmap_2d"] = "cache_" + _pm_other_suffix.lstrip("_")
        EMBEDDING_SOURCES["pacmap_3d"] = "cache_" + _pm_other_suffix.lstrip("_")
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
                    _safe_savez(".fleek_cache.npz", **cached)
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
            global HVG_NAMES, HVG_VAR, ALL_GENE_VAR, ALL_GENE_VAR_METRIC, GENE_CUTOFF_FLAG, X_CSC, CSC_BUILDING, CSC_CACHED, CSC_TIME, DISPLAY_CSC, DISPLAY_CSC_PATH
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
            # Try gene-variability disk cache first — the mean/var/detection
            # stats over a 50k-cell subsample of a sparse column matrix are
            # by far the expensive step (~30s at 1M cells), while the final
            # results are tiny (a float32 + int8 per gene). Cache fingerprint
            # guards against dataset changes: cell/gene counts, first/last
            # gene name, and the routed source matrix path.
            _var_from_cache = False
            _var_route_path = (MATRIX_ROUTING.get("variability", {}).get("path") or "X") if MATRIX_ROUTING else "X"
            if not _has_hvg_col:
                _vp_r, _ = _cache_read_path(".fleek_gene_var.npz")
                if _vp_r is not None and _vp_r.exists():
                    try:
                        _vc = np.load(str(_vp_r), allow_pickle=True)
                        _ok = (int(_vc["n_cells"]) == ADATA.n_obs and
                               int(_vc["n_genes"]) == ADATA.n_vars and
                               str(_vc["src_path"]) == _var_route_path and
                               str(_vc["first_gene"]) == str(ADATA.var_names[0]) and
                               str(_vc["last_gene"]) == str(ADATA.var_names[-1]))
                        if _ok:
                            ALL_GENE_VAR = np.asarray(_vc["all_gene_var"], dtype=np.float32)
                            GENE_CUTOFF_FLAG = np.asarray(_vc["gene_cutoff_flag"], dtype=np.int8)
                            ALL_GENE_VAR_METRIC = str(_vc["metric"])
                            _hvg_names = [str(x) for x in _vc["hvg_names"]]
                            HVG_NAMES = _hvg_names
                            HVG_VAR = {n: float(v) for n, v in zip(_hvg_names, _vc["hvg_values"])}
                            _var_from_cache = True
                            _n_flagged = int((GENE_CUTOFF_FLAG > 0).sum())
                            print(f"  Gene variability loaded from cache ({len(ALL_GENE_VAR):,} genes, {_n_flagged:,} flagged)")
                        else:
                            print("  Gene variability cache fingerprint mismatch, recomputing...")
                    except Exception as _ve:
                        print(f"  Gene variability cache load failed ({_ve}), recomputing...")

            if not _has_hvg_col and not _var_from_cache:
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

                # Save to disk — cheap (under 1 MB typically), huge payoff
                # on reload of multi-million-cell datasets.
                try:
                    _vp_w = _cache_write_path(".fleek_gene_var.npz")
                    _hvg_values = np.asarray([HVG_VAR[n] for n in HVG_NAMES], dtype=np.float32)
                    np.savez(str(_vp_w),
                             all_gene_var=ALL_GENE_VAR,
                             gene_cutoff_flag=GENE_CUTOFF_FLAG,
                             hvg_names=np.asarray(HVG_NAMES, dtype=object),
                             hvg_values=_hvg_values,
                             metric=ALL_GENE_VAR_METRIC,
                             n_cells=ADATA.n_obs,
                             n_genes=ADATA.n_vars,
                             src_path=_var_route_path,
                             first_gene=str(ADATA.var_names[0]),
                             last_gene=str(ADATA.var_names[-1]))
                    _vp_kb = os.path.getsize(str(_vp_w)) / 1024
                    print(f"  Gene variability cached ({_vp_kb:.0f} KB)")
                except Exception as _se:
                    print(f"  Gene variability cache save failed ({_se})")

            # Try loading CSC from cache
            import scipy.sparse as _sp
            loaded_from_cache = False
            if csc_cache_path.exists():
                print("  Loading cached CSC index...")
                try:
                    t1 = time.time()
                    # Prefer memory-mapped load: data/indices/indptr become
                    # np.memmap views and pages fault in only when a gene
                    # column is actually requested. Previously the full
                    # matrix (several GB for >1M-cell datasets) was read
                    # into RAM through Python's zipfile reader on every
                    # startup — that's the "long time to load" the user
                    # experienced. Falls back to scipy's load_npz on any
                    # parse error or compressed entry.
                    X_CSC = _load_csc_mmap(csc_cache_path)
                    if X_CSC is None:
                        X_CSC = _sp.load_npz(str(csc_cache_path))
                        _mmap_note = ""
                    else:
                        _mmap_note = " (mmap)"
                    if X_CSC.shape != ADATA.X.shape:
                        print(f"  CSC cache shape mismatch ({X_CSC.shape} vs {ADATA.X.shape}), rebuilding...")
                        X_CSC = None
                    else:
                        loaded_from_cache = True
                        CSC_CACHED = True
                        print(f"  CSC index loaded from cache{_mmap_note} ({time.time()-t1:.2f}s)")
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

            # If the gene-display route points somewhere other than adata.X
            # (common when log_normalized lives in a layer while X holds raw
            # counts), build a second CSC for that matrix so gene lookups on
            # the display route also stay O(nnz_column) instead of falling
            # into CSR column-slicing. Cached to disk separately.
            try:
                _dr = MATRIX_ROUTING.get("gene_display", {}) if MATRIX_ROUTING else {}
                _dpath = _dr.get("path") or "X"
                if _dpath != "X" and not BACKED:
                    _dmat = _get_matrix_by_path(_dpath)
                    if _dmat is not None and _sp.issparse(_dmat):
                        DISPLAY_CSC_PATH = _dpath
                        # Sanitize path into a filename-safe suffix.
                        _safe = _dpath.replace("/", "_").replace("[", "_").replace("]", "_").replace("'", "").replace('"', "")
                        _disp_suffix = f".display_csc_{_safe}.npz"
                        _disp_read, _ = _cache_read_path(_disp_suffix)
                        _loaded_disp = False
                        if _disp_read is not None:
                            try:
                                t2 = time.time()
                                DISPLAY_CSC = _load_csc_mmap(_disp_read)
                                if DISPLAY_CSC is None:
                                    DISPLAY_CSC = _sp.load_npz(str(_disp_read))
                                _disp_note = " (mmap)" if _is_mmap_backed(getattr(DISPLAY_CSC, "data", None)) else ""
                                if DISPLAY_CSC is not None and DISPLAY_CSC.shape == _dmat.shape:
                                    print(f"  Display CSC ({_dpath}) loaded from cache{_disp_note} ({time.time()-t2:.2f}s)")
                                    _loaded_disp = True
                                else:
                                    DISPLAY_CSC = None
                            except Exception:
                                DISPLAY_CSC = None
                        if not _loaded_disp:
                            t2 = time.time()
                            print(f"  Building display CSC index for route {_dpath!r}...")
                            DISPLAY_CSC = _dmat.tocsc()
                            print(f"  Display CSC ready ({time.time()-t2:.0f}s)")
                            try:
                                _disp_wp = _cache_write_path(_disp_suffix)
                                _sp.save_npz(str(_disp_wp), DISPLAY_CSC, compressed=False)
                                _dmb = os.path.getsize(str(_disp_wp)) / 1e6
                                print(f"  Display CSC cache saved ({_dmb:.0f} MB)")
                            except Exception as e:
                                print(f"  Display CSC cache save failed ({e}) — will rebuild next time")
            except Exception as e:
                print(f"  Display CSC build skipped ({e})")

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
    on the fly if the chosen matrix is raw/denoised counts.

    Fast-path philosophy: whenever the routed matrix has a CSC index available
    (either X_CSC, which covers adata.X, or DISPLAY_CSC, a lazily-built CSC
    for a routed non-X display matrix), we pull the column from the index in
    O(nnz_column) time and apply any transform AFTER extraction. Previous
    versions skipped the CSC path whenever a transform was requested, which
    forced the slow scipy-CSR column-slice (scans all nnz — ~100M+ entries
    for a 1.6M×20k matrix) and turned a single gene fetch into a 40+ second
    wait on large census datasets. Now the CSC is used even for the
    normalize+log1p display path.
    """
    import scipy.sparse as sp

    idx = GENE_INDEX.get(gene_name)
    if idx is None:
        return None

    route = MATRIX_ROUTING.get("gene_display", {}) if MATRIX_ROUTING else {}
    src_path = route.get("path") or "X"
    transform = route.get("transform", "none")

    col = None
    # Try a CSC-backed fast extraction first.
    if src_path == "X" and X_CSC is not None:
        col = X_CSC[:, idx].toarray().ravel().astype(np.float32)
    elif DISPLAY_CSC is not None and DISPLAY_CSC_PATH == src_path:
        col = DISPLAY_CSC[:, idx].toarray().ravel().astype(np.float32)

    # Fall back to matrix slicing if no CSC covers this route.
    if col is None:
        src_mat = _get_matrix_by_path(src_path) if src_path else None
        if src_mat is None:
            src_mat = ADATA.X
        if BACKED and src_path == "X":
            # Backed-mode batched reads of raw X (can't CSC a backed matrix).
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
            # WARNING: on large CSR matrices this is O(nnz) per call and is
            # the pathological slow path. DISPLAY_CSC should be built at
            # load time for any routed display matrix that isn't adata.X.
            print(f"  gene lookup slow path: no CSC for route {src_path!r}, column-slicing CSR")
            col_src = src_mat[:, idx]
            if sp.issparse(col_src):
                col = col_src.toarray().ravel()
            else:
                col = np.asarray(col_src).ravel()
            col = col.astype(np.float32)

    # Apply on-the-fly normalisation for routes that expose raw-ish counts.
    if transform == "normalize+log1p":
        src_mat = _get_matrix_by_path(src_path) if src_path else ADATA.X
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

    # Filename: <parent-stem>_subset_<group>_<N>cells_<YYYY-MM-DD_HH-MM-SS>.h5ad
    # The parent stem makes it obvious in a downloads folder which dataset
    # the subset came from; the timestamp guards against silent overwrites
    # when the user re-exports the same group with the same options.
    # Dashes + underscore are filename-safe on macOS/Linux/Windows/FAT.
    import datetime
    ts = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    safe_name = "".join(c if c.isalnum() or c in "._-" else "_" for c in group_name).strip("_") or "subset"
    parent_stem = Path(LOADED_PATH).stem if LOADED_PATH else "fleek"
    parent_safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in parent_stem).strip("_") or "fleek"
    fname = f"{parent_safe}_subset_{safe_name}_{len(cell_indices)}cells_{ts}.h5ad"
    # Remote mode: write subset under USER_DATA_DIR/exports so concurrent
    # users don't collide in /tmp and so the user's own quota (not the
    # system /tmp) is what bounds export volume. Single-user mode keeps
    # the historical /tmp path. The file is unlinked after streaming
    # in either case (see _serve_export), so neither is a long-term
    # storage commitment.
    exports_dir = _user_path("exports")
    if exports_dir is not None:
        fpath = str(exports_dir / fname)
    else:
        fpath = os.path.join(tempfile.gettempdir(), fname)
    # Quota guard: refuse the export pre-write if the user is already
    # at cap. We don't know the exact subset size until write_h5ad
    # finishes, but a hard pre-check at _quota_check_or_raise(0) at
    # least catches the "already over cap from previous artifacts"
    # case so we don't write a doomed file.
    _quota_check_or_raise(extra_bytes=0)
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
        elif path == "/api/server-info":
            self._serve_server_info()
        elif path == "/api/complete":
            self._serve_complete(params.get("partial", [""])[0])
        elif path == "/api/gene-var":
            self._serve_gene_var()
        elif path == "/api/cluster-genes":
            self._serve_cluster_genes(params)
        elif path == "/api/pathway-databases":
            self._serve_pathway_databases()
        elif path == "/api/gene-info":
            self._serve_gene_info(params.get("symbol", [""])[0])
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
        elif path == "/api/pathway-ora":
            self._serve_pathway_ora(req)
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
        """List directories and h5ad files for the file browser.

        In remote mode the browser is chrooted to USER_DATA_DIR ∪
        SHARED_DATA_DIR — any path the client sends is rejected unless
        it resolves under one of those two roots, and the special
        sentinel values "my" and "shared" map to those roots so the
        client doesn't have to know absolute paths. The "..": parent
        link is suppressed at the root of each chroot so the user
        can't walk above it (Path.parent stops at /, but our gate
        still rejects it so it'd appear and produce empty results)."""
        try:
            # Default + sentinel handling. In remote mode, "" / unset
            # defaults to the user's own data dir; the client can also
            # send "my" or "shared" to switch roots without typing a
            # full path.
            if REMOTE_MODE:
                if not req_path or req_path == "my":
                    req_path = str(USER_DATA_DIR)
                elif req_path == "shared":
                    req_path = str(SHARED_DATA_DIR)
                elif req_path == "~":
                    # ~ is meaningless on a server; redirect to user dir.
                    req_path = str(USER_DATA_DIR)
            else:
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

            # Remote-mode gate: anything outside the two roots is
            # rejected, including a path that .resolve()'d through a
            # symlink to an unrelated location. The client gets a
            # clear error instead of silent wrong data.
            if not _path_in_allowed_roots(str(p)):
                raise PermissionError(
                    f"Path not allowed in remote mode: {p}"
                )

            entries = []
            # Parent directory link — suppress when we'd walk above a
            # chroot boundary.
            parent = str(p.parent)
            at_chroot_root = REMOTE_MODE and (
                _is_under(p, USER_DATA_DIR) and Path(p) == USER_DATA_DIR
                or _is_under(p, SHARED_DATA_DIR) and Path(p) == SHARED_DATA_DIR
            )
            if parent != str(p) and not at_chroot_root:  # not at root, not at chroot
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
            # Remote mode: refuse completions outside the chroot. We
            # check the parent dir (the directory we'd iterate) since
            # a partial like "/data/users/alice/uplo" has parent
            # "/data/users/alice/" which IS allowed and should yield
            # "uploads/".
            if REMOTE_MODE:
                check = (p if (partial.endswith("/") or partial.endswith(os.sep)) else p.parent)
                if not _path_in_allowed_roots(str(check)):
                    data = json.dumps({"completions": []}).encode("utf-8")
                    self.send_response(200)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Content-Length", str(len(data)))
                    self.end_headers()
                    self.wfile.write(data)
                    return
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

    # Bump this when the gene-info payload schema changes so cached entries
    # from an older shape get re-fetched automatically. Entries without a
    # matching _v are treated as stale.
    _GENE_INFO_SCHEMA = 3

    def _serve_gene_info(self, symbol):
        """Fetch gene metadata from mygene.info for `symbol`, normalize into a
        compact payload, and cache per-organism to CACHE_DIR/gene_info_{species}.json
        so each symbol is only fetched once per organism across sessions.

        Fields returned: symbol, name, summary, summary_source, generifs,
        aliases, chromosome, map_location, type, entrez, ensembl, uniprot,
        hgnc, mim, species, found, _v. `generifs` is a list of up-to-4 short
        Gene Reference Into Function statements (NCBI-curated) used as a
        fallback when the Entrez summary field is empty — common for mouse
        genes and less-studied loci. `summary_source` is "entrez" or "generif"
        so the client can explain where the text came from.
        Errors degrade gracefully — a {"found": False, "error": ...} response
        still lets the client render the external-links grid (which needs only
        the symbol + organism)."""
        symbol = (symbol or "").strip()
        def _send(payload, status=200):
            data = json.dumps(payload).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
        if not symbol:
            _send({"error": "empty symbol", "found": False}); return
        # Map detected organism name -> mygene.info species key
        org_name = (_DETECTED_ORGANISM[0] if _DETECTED_ORGANISM else "Human").lower()
        species_map = {"human": "human", "mouse": "mouse", "rat": "rat",
                       "zebrafish": "zebrafish", "fruitfly": "fruitfly",
                       "nematode": "nematode"}
        species = species_map.get(org_name, "human")
        # Cross-dataset cache — keyed by organism, not dataset
        _init_cache_dir()
        cache_file = CACHE_DIR / f"gene_info_{species}.json" if CACHE_DIR else None
        cache = {}
        if cache_file and cache_file.exists():
            try:
                with open(cache_file) as f: cache = json.load(f)
            except Exception: cache = {}
        key = symbol.upper()
        cached = cache.get(key)
        if cached and cached.get("_v") == self._GENE_INFO_SCHEMA:
            _send(cached); return
        # Miss or stale-schema entry — fetch fresh
        try:
            import urllib.request, urllib.parse
            fields = ("symbol,name,summary,alias,other_names,genomic_pos,type_of_gene,"
                      "entrezgene,ensembl,uniprot,HGNC,MIM,map_location,generif")
            url = ("https://mygene.info/v3/query?q=" + urllib.parse.quote(symbol)
                   + "&species=" + species + "&fields=" + fields + "&size=1")
            req = urllib.request.Request(url, headers={"User-Agent": "RNA-FLEEK"})
            with urllib.request.urlopen(req, timeout=10) as r:
                data = json.loads(r.read().decode("utf-8"))
            hits = data.get("hits") or []
            if not hits:
                payload = {"symbol": symbol, "species": species, "found": False,
                           "_v": self._GENE_INFO_SCHEMA}
            else:
                h = hits[0]
                # The query endpoint returns a summary field in the hit, but when
                # genome annotations are rich it may only expose generif. If
                # summary is empty, fetch the full gene record too — /v3/gene/
                # returns everything including `entrezgene_summary` on occasion.
                pos = h.get("genomic_pos") or {}
                if isinstance(pos, list): pos = pos[0] if pos else {}
                chrom = str(pos.get("chr", "")) if isinstance(pos, dict) else ""
                ensembl = h.get("ensembl") or {}
                if isinstance(ensembl, list): ensembl = ensembl[0] if ensembl else {}
                ens_gene = ensembl.get("gene", "") if isinstance(ensembl, dict) else ""
                uniprot = h.get("uniprot") or {}
                up = ""
                if isinstance(uniprot, dict):
                    up = uniprot.get("Swiss-Prot") or uniprot.get("TrEMBL") or ""
                    if isinstance(up, list): up = up[0] if up else ""
                aliases = h.get("alias") or []
                if isinstance(aliases, str): aliases = [aliases]
                other_names = h.get("other_names") or []
                if isinstance(other_names, str): other_names = [other_names]
                # Extract generifs — each is {text, pubmed} or a string
                raw_generif = h.get("generif") or []
                if isinstance(raw_generif, dict): raw_generif = [raw_generif]
                generifs = []
                for g in raw_generif[:8]:
                    if isinstance(g, dict):
                        t = (g.get("text") or "").strip()
                        pm = g.get("pubmed") or g.get("pmid") or ""
                        if isinstance(pm, list): pm = pm[0] if pm else ""
                        if t:
                            generifs.append({"text": t, "pmid": str(pm)})
                    elif isinstance(g, str) and g.strip():
                        generifs.append({"text": g.strip(), "pmid": ""})
                # Summary priority:
                #   1. mygene.info `summary` (NCBI RefSeq-curated; strongest
                #      coverage for human).
                #   2. NCBI E-utilities `esummary` on db=gene — carries the
                #      Alliance-of-Genome-Resources / MGI / RGD / ZFIN text
                #      that NCBI Gene pages actually display for model
                #      organisms. mygene.info doesn't always sync this field,
                #      which is why Sort1 (mouse) and similar genes looked
                #      bare before. Quality match: what you see on the web
                #      page is what you get here.
                #   3. First 2-3 NCBI GeneRIFs assembled into a paragraph,
                #      kept as a last-resort fallback for genes that have
                #      neither of the above (rare but possible for lncRNAs
                #      and newly-annotated loci).
                entrez_summary = (h.get("summary") or "").strip()
                entrez_id = str(h.get("entrezgene") or "")
                summary_text = ""
                summary_source = ""
                if entrez_summary:
                    summary_text = entrez_summary
                    summary_source = "entrez"
                else:
                    ncbi_summary = _fetch_ncbi_gene_summary(entrez_id)
                    if ncbi_summary:
                        summary_text = ncbi_summary
                        summary_source = "ncbi"
                    elif generifs:
                        summary_text = " ".join(g["text"] for g in generifs[:3])
                        summary_source = "generif"
                payload = {
                    "symbol": h.get("symbol") or symbol,
                    "name": h.get("name") or "",
                    "summary": summary_text,
                    "summary_source": summary_source,
                    "generifs": generifs[:4],
                    "aliases": list(aliases),
                    "other_names": [n for n in other_names if isinstance(n, str)][:4],
                    "chromosome": chrom,
                    "map_location": str(h.get("map_location") or ""),
                    "type": h.get("type_of_gene") or "",
                    "entrez": entrez_id,
                    "ensembl": ens_gene,
                    "uniprot": up,
                    "hgnc": str(h.get("HGNC") or ""),
                    "mim": str(h.get("MIM") or ""),
                    "species": species,
                    "found": True,
                    "_v": self._GENE_INFO_SCHEMA,
                }
            if cache_file is not None:
                cache[key] = payload
                try:
                    with open(cache_file, "w") as f: json.dump(cache, f)
                except Exception as _e:
                    print(f"  gene-info cache write failed: {_e}")
            _send(payload)
        except Exception as e:
            _send({"symbol": symbol, "species": species, "found": False,
                   "error": str(e)})

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
            # Remote-mode chroot: the requested .h5ad must live under
            # the user's own dir or the shared dataset dir. Mirrors the
            # browse gate so a client can't load a file outside what
            # the browser would have shown them.
            if not _path_in_allowed_roots(fpath):
                raise PermissionError(
                    "This file is outside the allowed data roots."
                )
            fast = req.get("fast_umap", True)
            fast_sub = int(req.get("fast_umap_n", 50000))
            # Default 'auto' so a dataset that wouldn't otherwise fit in
            # RAM (dense .X bigger than ~40% of available memory) loads
            # in backed mode instead of OOMing the server. Small files
            # still go fully in-memory for speed. The frontend doesn't
            # currently expose a "backed mode" toggle in remote-mode UI;
            # this default is what makes "just click a big dataset and
            # it loads" work without manual intervention.
            backed = req.get("backed", "auto")
            _native_2d = req.get("native_2d", False)
            DEG_SETTINGS["fast"] = req.get("fast_deg", True)
            DEG_SETTINGS["max_cells"] = int(req.get("fast_deg_n", 50000))
            PACMAP_SETTINGS["fast"] = req.get("fast_pacmap", True)

            # Defensive reset: clear any leftover abort/error/aborted state
            # from a previous run BEFORE the background thread starts. Without
            # this, if a user aborted a prior load and then clicked Load
            # again quickly, the new load could either (a) inherit the
            # ABORT_REQUESTED=True flag and self-cancel before it begins, or
            # (b) the polling client could see the stale "aborted" or
            # "error" PROGRESS state and surface a confusing message.
            global ABORT_REQUESTED
            ABORT_REQUESTED = False
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
            # `model` may be a short key ('sonnet') or a full model ID;
            # _resolve_claude_model handles both. Defaults to Sonnet 4.6
            # when absent or unrecognized — see CLAUDE_MODELS for the
            # supported set + cost guidance the UI shows.
            model = req.get("model")
            result = annotate_clusters_llm(top_n=top_n, tissue_hint=tissue, model=model)
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
            # Remote mode: refuse to delete shared caches. Users can
            # still clear their own per-user caches (the dataset's
            # cache file in <user>/caches/ — _server_cache_path) but
            # the dataset-dir copy that lives in /data/shared/ is admin
            # property. The deletion loop below already handles the
            # case where one of the two paths errors with PermissionError
            # gracefully (RO filesystem perms would make the unlink fail
            # anyway); this gate fails fast with a clearer message.
            if REMOTE_MODE and _path_in_shared(LOADED_PATH):
                # Allow deletion only of the user-side server-cache copy
                # for the shared dataset; never touch the shared-dir copy.
                pass  # handled in the loop below by skipping shared paths
            cache_type = req.get("type", "")
            suffix_map = {
                "fleek": ".fleek_cache.npz",
                "csc": ".csc_cache.npz",
                "gene_var": ".fleek_gene_var.npz",
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
                    # In remote mode, refuse to even attempt a delete
                    # on shared caches — keeps the audit log clean and
                    # gives the user a clearer message than a silent
                    # PermissionError swallow.
                    if REMOTE_MODE and _path_in_shared(cp):
                        print(f"  Refusing to delete shared cache: {cp}")
                        continue
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
            # Remote mode: cache location is admin-controlled (set at
            # supervisor-spawn time). The endpoint stays accessible for
            # GET-style introspection but writes are refused so the
            # client can't redirect cache writes outside the user's
            # quota'd dir.
            if REMOTE_MODE and (new_mode or new_dir):
                raise PermissionError(
                    "Cache location is admin-controlled in remote mode."
                )
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

    def _serve_server_info(self):
        """Return server-mode metadata for the frontend so it can hide
        controls that don't apply in remote mode (cache-mode toggle,
        API-key save, browse-everywhere, etc.). Hit once at boot from
        the client; no per-request overhead. The 'restrictions' block
        is the canonical list — every entry maps to a hidden control
        in fleek.html. Adding a new restriction here without also
        teaching the frontend to honour it is a no-op (the endpoint
        still works); the lockdown actually lives on the server side."""
        info = {
            "remote_mode": REMOTE_MODE,
            "fleek_version": FLEEK_VERSION,
            # Catalogue of Claude models the user can pick from in the
            # Settings → Claude model dropdown. The frontend reads this
            # so adding a model server-side automatically appears in
            # the UI; cost_hint is shown inline as a per-option suffix
            # so the user can pick with eyes-open.
            "claude_models": [
                {"key": k, "id": m["id"], "label": m["label"], "cost_hint": m["cost_hint"]}
                for k, m in CLAUDE_MODELS.items()
            ],
            "claude_default_model": CLAUDE_DEFAULT_MODEL_KEY,
        }
        if REMOTE_MODE:
            info["user_quota_mb"] = USER_QUOTA_MB
            info["restrictions"] = {
                "cache_mode_toggle": True,    # hide Settings → cache-mode
                "cache_dir_edit":    True,    # hide Settings → cache-dir
                # api_key_edit: NOT in restrictions any more. In remote mode
                # /api/set-key writes a per-user <data-dir>/api_key file
                # (not the server-level ~/.fleek.env), so it's safe + useful
                # to expose the Settings → API key UI to users.
                "auto_quit_toggle":  True,    # admin-controlled
                "browse_everywhere": True,    # only show My data / Shared roots
                "delete_shared_caches": True, # × on shared cache badges hidden
            }
        result = json.dumps(info).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_set_key(self, req):
        """Save a new Anthropic API key.

        In single-user mode: writes to ~/.fleek.env (preserves other vars).
        In remote mode:      writes to <data-dir>/api_key — a per-user
                             file, mode 600. The supervisor reads this
                             file at next subprocess spawn (preferred over
                             its DB column), so the new key survives
                             idle-reap + respawn.

        Either way the in-memory ANTHROPIC_API_KEY global is updated so
        the change takes effect immediately for in-flight Claude calls.
        """
        global ANTHROPIC_API_KEY
        try:
            key = req.get("key", "").strip()
            if not key:
                raise ValueError("Key is empty.")
            if not key.startswith("sk-ant-"):
                raise ValueError("Key doesn't look like an Anthropic key (expected sk-ant-...).")
            if REMOTE_MODE:
                if USER_DATA_DIR is None:
                    raise RuntimeError("Remote mode without USER_DATA_DIR — refusing to write key")
                kf = USER_DATA_DIR / "api_key"
                kf.write_text(key)
                try:
                    kf.chmod(0o600)
                except OSError:
                    pass
                print(f"  Claude API key updated via settings UI (per-user → {kf}).")
            else:
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
                print("  Claude API key updated via settings UI (~/.fleek.env).")
            ANTHROPIC_API_KEY = key
            result = json.dumps({"ok": True}).encode("utf-8")
        except Exception as e:
            result = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        self.wfile.write(result)

    def _serve_clear_key(self):
        """Remove the Anthropic API key from memory and disk.
        Single-user mode: edits ~/.fleek.env; remote mode: deletes
        <data-dir>/api_key. ANTHROPIC_API_KEY is cleared in either case."""
        global ANTHROPIC_API_KEY
        try:
            ANTHROPIC_API_KEY = ""
            if REMOTE_MODE:
                if USER_DATA_DIR is not None:
                    kf = USER_DATA_DIR / "api_key"
                    if kf.exists():
                        try:
                            kf.unlink()
                        except OSError:
                            pass
                print("  Claude API key cleared (per-user file removed).")
            else:
                env_path = Path.home() / ".fleek.env"
                if env_path.exists():
                    lines = [l for l in env_path.read_text().splitlines()
                             if not l.strip().startswith("ANTHROPIC_API_KEY=")]
                    if lines:
                        env_path.write_text("\n".join(lines) + "\n")
                    else:
                        env_path.unlink()
                print("  Claude API key cleared (~/.fleek.env).")
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
        The model is taken from req["model"] (resolved via _resolve_claude_model)
        so the user can confirm their key works for the model they actually
        plan to use — testing the cheap one and getting a 200 doesn't
        guarantee the expensive one is enabled on the same key (Anthropic
        permissions are model-scoped). Falls back to Haiku when no model
        is sent."""
        import urllib.request, urllib.error
        try:
            candidate = (req.get("key") or "").strip() if isinstance(req, dict) else ""
            key = candidate or ANTHROPIC_API_KEY
            if not key:
                raise ValueError("No key to test. Paste a key into the field, or save one first.")
            if not key.startswith("sk-ant-"):
                raise ValueError("Key doesn't look like an Anthropic key (expected sk-ant-...).")
            # Pick the model the user wants to verify. Default Haiku
            # only when nothing was sent — preserves backwards-compat
            # for any pre-v0.5.11 client that doesn't include `model`.
            test_model_id = (
                _resolve_claude_model(req)
                if isinstance(req, dict) and req.get("model")
                else "claude-haiku-4-5-20251001"
            )
            body = json.dumps({
                "model": test_model_id,
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
        """Save a named session. Reject silent overwrites: if a session file
        with the same name already exists on disk and the request didn't pass
        `overwrite: true`, return {"ok":false, "code":"exists"} so the client
        can prompt and retry explicitly. The auto-save slot "_auto" is exempt
        — it's the point of an auto-save to overwrite itself on every write."""
        try:
            if not LOADED_PATH:
                raise ValueError("No dataset loaded")
            name = req.get("name", "").strip()
            if not name:
                raise ValueError("Session name is empty")
            overwrite = bool(req.get("overwrite", False))
            suffix = self._session_suffix(name)
            # Existence check across both dataset and server cache dirs
            existing, _ = _cache_read_path(suffix)
            if existing is not None and not overwrite and name != "_auto":
                result = json.dumps({"ok": False, "code": "exists",
                                     "error": f"Session '{name}' already exists"}).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(result)))
                self.end_headers()
                self.wfile.write(result)
                return
            state = req.get("state", {})
            state["_meta"] = {"n_cells": N_CELLS, "n_clusters": len(CLUSTER_NAMES or []),
                              "dataset": Path(LOADED_PATH).name,
                              "saved": time.time()}
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

    def _serve_pathway_databases(self):
        """List the pathway databases currently loaded for the dataset's
        organism. Cheap response — just metadata, no gene lists."""
        try:
            entries = []
            for db_id in sorted(PATHWAY_DBS.keys()):
                d = PATHWAY_DBS[db_id]
                entries.append({
                    "id": db_id,
                    "label": d.get("label", db_id),
                    "source": d.get("source", ""),
                    "n_pathways": d.get("n_pathways", 0),
                })
            payload = {
                "ok": True,
                "organism": PATHWAY_ORGANISM,
                "databases": entries,
            }
            data = json.dumps(payload).encode("utf-8")
        except Exception as e:
            data = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _serve_pathway_ora(self, req):
        """Over-representation analysis (hypergeometric).

        For each pathway in the chosen database we test the null
        hypothesis that the user's input gene list is a random sample
        from the background. The test is one-sided (we care about
        OVER-representation, not under-). p-values are BH-corrected
        across all pathways tested in this run.

        Inputs:
          genes:        ["ACTB", "MYH7", ...]    REQUIRED
          database:     "hallmark" | "go_bp" | …  REQUIRED
          background:   ["ACTB", ...] | null      OPTIONAL — null = use
                                                   the dataset's full gene list
          min_set_size: int (default 5)           skip tiny sets
          max_set_size: int (default 500)         skip huge sets
        """
        import hashlib
        import time as _t
        t0 = _t.time()
        try:
            genes_in = req.get("genes") or []
            db_id = (req.get("database") or "").strip()
            background_in = req.get("background")
            min_set = int(req.get("min_set_size", 5))
            max_set = int(req.get("max_set_size", 500))
            if not genes_in or not isinstance(genes_in, list):
                raise ValueError("'genes' must be a non-empty list")
            if db_id not in PATHWAY_DBS:
                avail = sorted(PATHWAY_DBS.keys())
                raise ValueError(f"Unknown database '{db_id}'. Available: {avail}")
            db = PATHWAY_DBS[db_id]
            terms = db.get("terms", {})
            # Background: explicit list if provided, else the dataset's
            # full gene list. Used to scope the hypergeometric universe.
            if background_in and isinstance(background_in, list) and len(background_in) > 0:
                bg_set = set(g for g in background_in if isinstance(g, str))
            elif GENE_NAMES_LIST:
                bg_set = set(GENE_NAMES_LIST)
            else:
                raise ValueError("No background available; load a dataset or pass 'background'.")
            if not bg_set:
                raise ValueError("Background is empty after deduplication.")
            # Case-insensitive matching: build an upper->canonical map
            # over the background. Pathway gene names from many sources
            # (MSigDB, Reactome) are uppercase regardless of organism;
            # mouse datasets typically use Title-case symbols. Without
            # this mapping, mouse + human-style GMT files would get zero
            # hits. The same trick is used by the existing GO loader.
            bg_upper = {}
            for g in bg_set:
                bg_upper[g.upper()] = g
            def _canon(g):
                return bg_upper.get(g.upper())
            input_set = set()
            for g in genes_in:
                if not isinstance(g, str):
                    continue
                c = _canon(g)
                if c is not None:
                    input_set.add(c)
            n = len(input_set)
            N = len(bg_set)
            warnings = []
            if n < 3:
                warnings.append(f"Input list resolves to only {n} genes after background intersection — results unreliable.")
            # Cache key: sorted input genes + db + background size.
            # Genes hash is stable across reorderings; background size is
            # a cheap proxy for "same dataset + same explicit background".
            key_genes = ",".join(sorted(input_set))
            key_bg = "" if background_in else "DATASET:" + str(N)
            cache_key = hashlib.sha1((db_id + "|" + key_genes + "|" + key_bg).encode("utf-8")).hexdigest()
            if cache_key in _PATHWAY_ORA_CACHE:
                cached = _PATHWAY_ORA_CACHE[cache_key]
                cached["cached"] = True
                cached["elapsed_s"] = round(_t.time() - t0, 4)
                data = json.dumps(cached).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)
                return
            # Run the test. SciPy's hypergeom is exact and fast — even at
            # 30k GO BP terms, this runs in ~30 ms.
            from scipy.stats import hypergeom
            results = []
            n_tested = 0
            for pid, term in terms.items():
                gene_names = term.get("genes") or []
                # Intersect pathway with background. Map case-insensitively.
                pathway_in_bg = set()
                for g in gene_names:
                    c = _canon(g)
                    if c is not None:
                        pathway_in_bg.add(c)
                K = len(pathway_in_bg)
                if K < min_set or K > max_set:
                    continue
                overlap = sorted(input_set & pathway_in_bg)
                k = len(overlap)
                if k == 0:
                    continue
                # Expected overlap under null = n*K/N. log2 fold-enrichment.
                expected = (n * K) / N if N > 0 else 0.0
                if expected > 0 and k > 0:
                    log2fe = float(np.log2(k / expected))
                else:
                    log2fe = 0.0
                # P(X >= k | hypergeom)  ==  hypergeom.sf(k - 1, N, K, n)
                p = float(hypergeom.sf(k - 1, N, K, n))
                results.append({
                    "id": pid,
                    "name": term.get("name", pid),
                    "description": term.get("description", ""),
                    "url": term.get("url", ""),
                    "n_pathway": K,
                    "n_overlap": k,
                    "log2fe": log2fe,
                    "p_value": p,
                    "p_adj": p,  # placeholder, BH applied next
                    "overlap_genes": overlap,
                })
                n_tested += 1
            # BH correction across all tested pathways. p_sorted_idx[i] = original index of i-th smallest p.
            if results:
                m = len(results)
                idx_sorted = sorted(range(m), key=lambda i: results[i]["p_value"])
                # Compute BH-adjusted p-values (Benjamini-Hochberg step-up)
                prev = 1.0
                for rank in range(m - 1, -1, -1):
                    i = idx_sorted[rank]
                    p = results[i]["p_value"]
                    adj = min(1.0, p * m / (rank + 1))
                    if adj < prev:
                        prev = adj
                    results[i]["p_adj"] = float(prev)
                # Sort by p_adj asc (stable, with p_value as tiebreak)
                results.sort(key=lambda r: (r["p_adj"], r["p_value"]))
            payload = {
                "ok": True,
                "database": db_id,
                "database_label": db.get("label", db_id),
                "n_input": n,
                "n_input_raw": len(genes_in),
                "n_background": N,
                "n_pathways_tested": n_tested,
                "n_results": len(results),
                "results": results,
                "warnings": warnings,
                "elapsed_s": round(_t.time() - t0, 4),
                "cached": False,
            }
            # Bound the cache so a runaway "paste a fresh list and run"
            # workflow doesn't grow it without limit. ~16 entries is fine.
            if len(_PATHWAY_ORA_CACHE) >= 16:
                # Drop oldest (insertion order). Python dict iteration is
                # insertion-ordered since 3.7.
                _PATHWAY_ORA_CACHE.pop(next(iter(_PATHWAY_ORA_CACHE)))
            _PATHWAY_ORA_CACHE[cache_key] = dict(payload)
            data = json.dumps(payload).encode("utf-8")
        except Exception as e:
            data = json.dumps({"ok": False, "error": str(e)}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

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
            # Same model-resolution pattern as the annotation endpoint —
            # short key ('sonnet') or full model ID, default to Sonnet
            # 4.6 if absent / unrecognized.
            lineage_model_id = _resolve_claude_model(req)

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
                "model": lineage_model_id,
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
                            _safe_savez(".fleek_cache.npz", **cached)
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
        """Client heartbeat — resets the auto-quit timer (unload + process exit)."""
        global AUTO_UNLOAD, AUTO_UNLOAD_TIMEOUT, _LAST_HEARTBEAT, _HEARTBEAT_TIMER
        _LAST_HEARTBEAT = time.time()
        # Cancel any pending timer
        if _HEARTBEAT_TIMER:
            _HEARTBEAT_TIMER.cancel()
            _HEARTBEAT_TIMER = None
        # Parse settings updates from client
        if req:
            if "auto_unload" in req:
                AUTO_UNLOAD = bool(req["auto_unload"])
            if "timeout" in req:
                AUTO_UNLOAD_TIMEOUT = max(0.5, float(req["timeout"]))
        # Schedule next quit check — runs even when no dataset is loaded, since
        # the server should shut down after a disconnect regardless. The timer
        # frees the dataset (if any) and then exits the process so the user
        # doesn't need to know a background server exists.
        if AUTO_UNLOAD:
            def _check_quit():
                global _HEARTBEAT_TIMER
                _HEARTBEAT_TIMER = None
                elapsed = time.time() - _LAST_HEARTBEAT
                if elapsed >= AUTO_UNLOAD_TIMEOUT and AUTO_UNLOAD:
                    if ADATA is not None:
                        print(f"  Auto-quit: unloading dataset (no heartbeat for {elapsed:.1f}s)")
                        try:
                            _reset_all()
                        except Exception as _e:
                            print(f"  Auto-quit: _reset_all failed ({_e})")
                        import gc; gc.collect()
                    print(f"  Auto-quit: stopping FLEEK server (no heartbeat for {elapsed:.1f}s)")
                    # os._exit avoids running atexit handlers / blocking on
                    # lingering threads. Safe here — we've already cleaned up.
                    import os as _os
                    _os._exit(0)
            _HEARTBEAT_TIMER = threading.Timer(AUTO_UNLOAD_TIMEOUT + 0.5, _check_quit)
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
            # Sanitize so a malicious X-Filename can't path-escape
            # (e.g. "../../etc/foo.h5ad"). basename() strips any
            # directory components; in remote mode this is essential
            # because the destination is a user-specific directory.
            fname = os.path.basename(fname)
            # Refuse early if accepting this upload would push the
            # user over quota. Length comes from Content-Length so
            # we know the size before reading.
            _quota_check_or_raise(extra_bytes=length)
            _progress("Receiving upload...", 1)
            # In remote mode the upload lands in the user's dedicated
            # uploads directory — no /tmp collisions between concurrent
            # users uploading files with the same name. Single-user
            # mode preserves the historical /tmp behaviour.
            uploads_dir = _user_path("uploads")
            if uploads_dir is not None:
                fpath = str(uploads_dir / fname)
            else:
                fpath = os.path.join(tempfile.gettempdir(), fname)
            received = 0
            with open(fpath, "wb") as f:
                while received < length:
                    chunk = self.rfile.read(min(65536, length - received))
                    if not chunk:
                        break
                    f.write(chunk)
                    received += len(chunk)
            # See the matching reset in _serve_load — same rationale: clear
            # any leftover abort/error state from a prior run so the new
            # load isn't pre-cancelled or surfacing stale UI text.
            global ABORT_REQUESTED
            ABORT_REQUESTED = False
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
    parser.add_argument("--backed", type=str, default="auto", choices=["auto", "on", "off"],
                        help="Backed mode: auto (detect; default — falls back to disk-backed read when the file is bigger than ~40%% of available RAM), on (force disk), off (force fully in-memory; will OOM on huge datasets)")
    parser.add_argument("--cache-mode", type=str, default="dataset", choices=["dataset", "server"],
                        help="Where to store caches: 'dataset' (alongside h5ad, default) or 'server' (~/.fleek_cache)")
    parser.add_argument("--cache-dir", type=str, default=None,
                        help="Custom server cache directory (default: ~/.fleek_cache)")
    parser.add_argument("--no-browser", action="store_true",
                        help="Don't auto-open browser (useful for remote/SSH sessions)")
    # ── Remote-server mode ──────────────────────────────────────────
    # --remote-mode is the single switch that flips FLEEK from "local
    # single-user tool" to "shared multi-user lab server". When set,
    # --data-dir and --shared-data are required and provide the per-
    # user writable root and the canonical RO dataset root respectively.
    # In Phase 1 these flags only affect filesystem routing (uploads,
    # exports, caches); endpoint lockdowns + frontend awareness arrive
    # in Phase 1B / 1C. The supervisor (Phase 2) is what actually sets
    # one --data-dir per user; manual invocation can also point this
    # at a single test directory.
    parser.add_argument("--remote-mode", action="store_true",
                        help="Run in multi-user-server mode (per-user data dir, shared dataset root, quotas, locked-down endpoints).")
    parser.add_argument("--data-dir", type=str, default=None,
                        help="Per-user writable root (required with --remote-mode). Holds caches, sessions, uploads, exports, settings.")
    parser.add_argument("--shared-data", type=str, default=None,
                        help="Read-only shared dataset root (required with --remote-mode). Admin pre-populates with shared atlases + caches.")
    parser.add_argument("--user-quota-mb", type=int, default=100,
                        help="Per-user write quota in MB (default 100). Applies in --remote-mode.")
    args = parser.parse_args()

    global CACHE_MODE, CACHE_DIR
    global REMOTE_MODE, USER_DATA_DIR, SHARED_DATA_DIR, USER_QUOTA_MB
    if args.remote_mode:
        # Validate required-with-remote-mode args. Failing fast here is
        # better than discovering at first request that --data-dir is
        # None and writes are silently going to /tmp.
        if not args.data_dir or not args.shared_data:
            parser.error("--remote-mode requires both --data-dir and --shared-data")
        REMOTE_MODE = True
        USER_DATA_DIR = Path(args.data_dir).expanduser().resolve()
        SHARED_DATA_DIR = Path(args.shared_data).expanduser().resolve()
        USER_QUOTA_MB = max(1, int(args.user_quota_mb))
        # The user dir must be writable; the shared dir must exist (RO
        # is fine — admin owns it). Create the user dir if missing so
        # a brand-new user signup doesn't have to pre-mkdir.
        USER_DATA_DIR.mkdir(parents=True, exist_ok=True)
        if not SHARED_DATA_DIR.exists():
            parser.error(f"--shared-data path does not exist: {SHARED_DATA_DIR}")
        # In remote mode, default --cache-dir to <user-dir>/caches/
        # so every cache the user generates lands in their own quota'd
        # space. The user can still override with an explicit --cache-dir
        # for testing, but the supervisor (Phase 2) won't.
        if not args.cache_dir:
            args.cache_dir = str(USER_DATA_DIR / "caches")
        # Force CACHE_MODE to "server" — the dataset-dir cache write
        # path in remote mode often hits the read-only shared dir and
        # falls back to CACHE_DIR anyway; explicitly using server mode
        # makes the intent obvious and avoids surprising fall-throughs.
        if args.cache_mode == "dataset":
            args.cache_mode = "server"
    CACHE_MODE = args.cache_mode
    if args.cache_dir:
        CACHE_DIR = Path(args.cache_dir)
    _init_cache_dir()

    # Re-load the Anthropic API key now that USER_DATA_DIR + REMOTE_MODE
    # globals are set. The first _load_api_key() call at module import
    # only saw the env + ~/.fleek.env paths; this second call also
    # picks up <data-dir>/api_key for remote-mode personal keys (set
    # by the user via /api/set-key from the FLEEK Settings UI).
    global ANTHROPIC_API_KEY
    ANTHROPIC_API_KEY = _load_api_key()

    dims = [int(d) for d in args.dims.split(",")]

    if args.h5ad:
        load_and_prepare(args.h5ad, max_cells=args.max_cells, n_dims_list=dims,
                         fast_umap=args.fast_umap, fast_umap_subsample=args.fast_umap_n,
                         backed=args.backed)

    server = type('ThreadingHTTPServer', (ThreadingMixIn, HTTPServer), {'daemon_threads': True})((args.host, args.port), FleekHandler)
    url = f"http://{args.host}:{args.port}"
    print(f"\n{'='*50}")
    print(f"  RNA-FLEEK v{FLEEK_VERSION} running at: {url}")
    if ADATA is not None:
        print(f"  {N_CELLS:,} cells · {len(CLUSTER_NAMES)} types · {len(GENE_NAMES_LIST):,} genes")
    else:
        print(f"  No dataset loaded — upload via browser")
    key_status = "loaded" if ANTHROPIC_API_KEY else "not set (Claude annotation disabled)"
    print(f"  Claude API key: {key_status}")
    if REMOTE_MODE:
        print(f"  Remote mode: ON")
        print(f"    user data: {USER_DATA_DIR}")
        print(f"    shared:    {SHARED_DATA_DIR}  (read-only)")
        print(f"    quota:     {USER_QUOTA_MB} MB / user")
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