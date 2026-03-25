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
import os
import struct
import sys
import time
import threading
from http.server import HTTPServer, SimpleHTTPRequestHandler
from socketserver import ThreadingMixIn
from pathlib import Path
from urllib.parse import urlparse, parse_qs

import numpy as np

# Globals set at startup
ADATA = None
UMAP_2D = None
UMAP_3D = None
CLUSTER_IDS = None
CLUSTER_NAMES = None
GENE_NAMES_LIST = None
N_CELLS = 0
HTML_DIR = None
PROGRESS = {"status": "idle", "message": "", "pct": 0}  # for client polling
PROCESSING_LOCK = threading.Lock()
BACKED = False  # True if adata.X is on disk (backed mode)
X_CSC = None  # CSC copy of expression matrix for fast column access
HVG_NAMES = []  # ordered list of HVG names for the client
GENE_INDEX = {}  # gene_name -> column index
CSC_BUILDING = False  # True while background CSC build is running
CSC_CACHED = False    # True if CSC was loaded from cache
CSC_TIME = 0          # seconds taken to build/load CSC index
LOAD_SETTINGS = {}  # Sent to client so it knows what options were active
LOADED_PATH = ""    # path to currently loaded h5ad, for annotation caching
MARKER_DB = None    # {cell_type: [genes]} — loaded from file or built-in

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
    # Look next to server script and in working directory
    script_dir = Path(__file__).parent
    paths_to_try.extend([
        script_dir / "cell_markers.json",
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


def annotate_clusters(top_n=50, test="wilcoxon"):
    """Run one-vs-rest DEG for each cluster, score against marker database.
    
    Returns list of {cluster_id, cluster_name, predictions: [{cell_type, score, pval, markers_found, markers_total}]}
    """
    import scanpy as sc
    import scipy.sparse as sp
    from scipy.stats import fisher_exact
    
    if ADATA is None or MARKER_DB is None:
        return {"error": "No data or marker database loaded"}
    
    if not MARKER_DB:
        load_marker_db()
    
    # Check cache
    cache_path = Path(LOADED_PATH).with_suffix(".annot_markers.json") if LOADED_PATH else None
    if cache_path and cache_path.exists():
        try:
            with open(cache_path) as f:
                cached = json.load(f)
            if cached.get("cluster_names") == CLUSTER_NAMES and cached.get("n_cells") == N_CELLS:
                print(f"  Marker annotation loaded from cache")
                cached["cached"] = True
                return cached
        except Exception as e:
            print(f"  Marker annotation cache load failed: {e}")
    
    t0 = time.time()
    
    # Create a temporary AnnData with cluster labels for rank_genes_groups
    # Work on a copy to avoid modifying the original
    adata = ADATA
    if "fleek_cluster" not in adata.obs.columns:
        adata.obs["fleek_cluster"] = [CLUSTER_NAMES[cid] for cid in CLUSTER_IDS]
    
    # Ensure we have a usable X (normalize if raw counts)
    # Check if data looks like raw counts
    if sp.issparse(adata.X):
        sample = adata.X[:min(1000, adata.n_obs)].toarray()
    else:
        sample = np.asarray(adata.X[:min(1000, adata.n_obs)])
    is_raw = np.allclose(sample, sample.astype(int)) and sample.max() > 20
    
    if is_raw:
        print("  Data appears to be raw counts — normalizing for DEG...")
        # Create a temporary normalized copy
        adata_work = adata.copy() if not BACKED else adata
        if not BACKED:
            sc.pp.normalize_total(adata_work, target_sum=1e4)
            sc.pp.log1p(adata_work)
        use_raw = False
    else:
        adata_work = adata
        use_raw = False
    
    print(f"  Running one-vs-rest {test} for {len(CLUSTER_NAMES)} clusters...")
    
    # Filter out clusters with too few cells for statistics
    min_cells = 5
    from collections import Counter
    cluster_counts = Counter(adata_work.obs["fleek_cluster"])
    valid_clusters = {c for c, n in cluster_counts.items() if n >= min_cells}
    small_clusters = {c for c, n in cluster_counts.items() if n < min_cells}
    if small_clusters:
        print(f"  Skipping {len(small_clusters)} clusters with <{min_cells} cells")
    mask = adata_work.obs["fleek_cluster"].isin(valid_clusters)
    adata_deg = adata_work[mask].copy()
    
    try:
        sc.tl.rank_genes_groups(adata_deg, groupby="fleek_cluster", method=test,
                                n_genes=top_n, use_raw=use_raw)
    except Exception as e:
        print(f"  rank_genes_groups failed: {e}, retrying with t-test...")
        try:
            sc.tl.rank_genes_groups(adata_deg, groupby="fleek_cluster", method="t-test",
                                    n_genes=top_n, use_raw=use_raw)
        except Exception as e2:
            return {"error": f"DEG failed: {str(e2)}"}
    
    # Build gene universe size (for Fisher's test)
    n_universe = len(GENE_NAMES_LIST)
    
    results = []
    for cid, cname in enumerate(CLUSTER_NAMES):
        n_cells_cluster = int(np.sum(CLUSTER_IDS == cid))
        # Get top upregulated genes for this cluster
        top_genes = set()
        if cname not in small_clusters:
            try:
                names = adata_deg.uns["rank_genes_groups"]["names"][cname][:top_n]
                scores = adata_deg.uns["rank_genes_groups"]["scores"][cname][:top_n]
                for i, (g, s) in enumerate(zip(names, scores)):
                    if s > 0:
                        top_genes.add(str(g).upper())
            except (KeyError, IndexError):
                pass
        
        if not top_genes:
            results.append({
                "cluster_id": cid,
                "cluster_name": cname,
                "n_cells": n_cells_cluster,
                "top_genes": [],
                "predictions": []
            })
            continue
        
        # Score against each cell type in the marker database
        predictions = []
        for cell_type, markers in MARKER_DB.items():
            marker_set = set(m.upper() for m in markers)
            if len(marker_set) < 2:
                continue  # skip cell types with too few markers
            
            overlap = top_genes & marker_set
            if not overlap:
                continue
            
            # Overlap coefficient: fraction of known markers found
            score = len(overlap) / len(marker_set)
            
            # Fisher's exact test for enrichment
            # 2x2 table: [found_in_both, in_markers_not_top] / [in_top_not_markers, neither]
            a = len(overlap)
            b = len(marker_set) - a
            c = len(top_genes) - a
            d = n_universe - a - b - c
            _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            
            predictions.append({
                "cell_type": cell_type,
                "score": round(score, 3),
                "pval": float(f"{pval:.2e}"),
                "markers_found": sorted(overlap),
                "markers_total": len(marker_set),
                "n_found": len(overlap)
            })
        
        # Sort by score descending, then by p-value ascending
        predictions.sort(key=lambda x: (-x["score"], x["pval"]))
        
        results.append({
            "cluster_id": cid,
            "cluster_name": cname,
            "n_cells": int(np.sum(CLUSTER_IDS == cid)),
            "top_genes": sorted(top_genes)[:20],  # send top 20 for display
            "predictions": predictions[:10]  # top 10 predictions
        })
    
    elapsed = round(time.time() - t0, 1)
    print(f"  Annotation complete ({elapsed}s) — {len(results)} clusters scored")
    
    # Clean up temporary objects
    if adata_deg is not adata_work and adata_deg is not adata:
        del adata_deg
    if adata_work is not adata:
        del adata_work
    
    result = {
        "results": results,
        "elapsed": elapsed,
        "n_clusters": len(CLUSTER_NAMES),
        "marker_db_size": len(MARKER_DB),
        "top_n": top_n,
        "cluster_names": CLUSTER_NAMES,
        "n_cells": N_CELLS
    }
    
    # Save cache
    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(result, f, separators=(",", ":"))
            print(f"  Marker annotation cached to {cache_path}")
        except Exception as e:
            print(f"  Marker annotation cache save failed: {e}")
    
    return result


# ── LLM-based cell type annotation ──

ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY", "")

def _get_cluster_top_genes(top_n=30, test="wilcoxon"):
    """Run one-vs-rest DEG and return {cluster_name: [top_genes]} dict."""
    import scanpy as sc
    import scipy.sparse as sp

    if ADATA is None:
        return {}

    adata = ADATA
    if "fleek_cluster" not in adata.obs.columns:
        adata.obs["fleek_cluster"] = [CLUSTER_NAMES[cid] for cid in CLUSTER_IDS]

    # Normalize if raw counts
    if sp.issparse(adata.X):
        sample = adata.X[:min(1000, adata.n_obs)].toarray()
    else:
        sample = np.asarray(adata.X[:min(1000, adata.n_obs)])
    is_raw = np.allclose(sample, sample.astype(int)) and sample.max() > 20

    if is_raw:
        adata_work = adata.copy() if not BACKED else adata
        if not BACKED:
            sc.pp.normalize_total(adata_work, target_sum=1e4)
            sc.pp.log1p(adata_work)
    else:
        adata_work = adata

    print(f"  Running one-vs-rest {test} for {len(CLUSTER_NAMES)} clusters...")
    
    # Filter out clusters with too few cells
    from collections import Counter
    cluster_counts = Counter(adata_work.obs["fleek_cluster"])
    valid_clusters = {c for c, n in cluster_counts.items() if n >= 5}
    small_clusters = {c for c, n in cluster_counts.items() if n < 5}
    if small_clusters:
        print(f"  Skipping {len(small_clusters)} clusters with <5 cells")
    mask = adata_work.obs["fleek_cluster"].isin(valid_clusters)
    adata_deg = adata_work[mask].copy()
    
    try:
        sc.tl.rank_genes_groups(adata_deg, groupby="fleek_cluster", method=test,
                                n_genes=top_n, use_raw=False)
    except Exception as e:
        print(f"  rank_genes_groups failed: {e}, retrying with t-test...")
        try:
            sc.tl.rank_genes_groups(adata_deg, groupby="fleek_cluster", method="t-test",
                                    n_genes=top_n, use_raw=False)
        except Exception as e2:
            return {}

    result = {}
    for cid, cname in enumerate(CLUSTER_NAMES):
        if cname in small_clusters:
            result[cname] = []
            continue
        try:
            names = adata_deg.uns["rank_genes_groups"]["names"][cname][:top_n]
            scores = adata_deg.uns["rank_genes_groups"]["scores"][cname][:top_n]
            genes = [str(g) for g, s in zip(names, scores) if s > 0]
            result[cname] = genes[:top_n]
        except (KeyError, IndexError):
            result[cname] = []

    if adata_deg is not adata_work and adata_deg is not adata:
        del adata_deg

    if adata_work is not adata:
        del adata_work

    return result


def annotate_clusters_llm(top_n=30, test="wilcoxon", tissue_hint=""):
    """Use Claude API to infer cell types from top marker genes per cluster.
    Batches clusters into groups to avoid token limits.

    Returns same format as annotate_clusters for UI compatibility.
    """
    import urllib.request
    import urllib.error

    if not ANTHROPIC_API_KEY:
        return {"error": "No API key. Set ANTHROPIC_API_KEY environment variable or pass --api-key."}

    if ADATA is None:
        return {"error": "No dataset loaded"}

    # Check cache
    cache_path = Path(LOADED_PATH).with_suffix(".annot_llm.json") if LOADED_PATH else None
    if cache_path and cache_path.exists():
        try:
            with open(cache_path) as f:
                cached = json.load(f)
            if cached.get("cluster_names") == CLUSTER_NAMES and cached.get("n_cells") == N_CELLS:
                print(f"  LLM annotation loaded from cache")
                cached["cached"] = True
                return cached
        except Exception as e:
            print(f"  LLM annotation cache load failed: {e}")

    t0 = time.time()

    # Step 1: Get top genes per cluster
    print("  LLM annotation: computing top genes per cluster...")
    cluster_genes = _get_cluster_top_genes(top_n=top_n, test=test)

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

    # Save cache
    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(result, f, separators=(",", ":"))
            print(f"  LLM annotation cached to {cache_path}")
        except Exception as e:
            print(f"  LLM annotation cache save failed: {e}")

    return result

def _progress(msg, pct=None):
    """Update global progress for client polling."""
    PROGRESS["message"] = msg
    if pct is not None:
        PROGRESS["pct"] = pct
    PROGRESS["status"] = "processing"
    print(f"  [{pct or '?'}%] {msg}")


def _get_available_ram_gb():
    """Estimate available RAM in GB."""
    try:
        import psutil
        return psutil.virtual_memory().available / 1e9
    except ImportError:
        # Assume 16GB if psutil not available
        return 16.0


def load_and_prepare(h5ad_path, max_cells=0, n_dims_list=[2, 3], fast_umap=False,
                     fast_umap_subsample=50000, backed="off"):
    """Load h5ad, subsample if needed, compute UMAP(s).
    
    backed: "auto" (use file size vs RAM), "on" (force backed), "off" (force in-memory)
    """
    global ADATA, UMAP_2D, UMAP_3D, CLUSTER_IDS, CLUSTER_NAMES, GENE_NAMES_LIST, N_CELLS, BACKED, X_CSC, HVG_NAMES, GENE_INDEX, CSC_BUILDING, CSC_CACHED, CSC_TIME, LOADED_PATH

    import scanpy as sc
    import scipy.sparse as sp
    import gc

    # Free old dataset before loading new one
    if ADATA is not None:
        _progress("Freeing previous dataset...", 2)
        ADATA = None
        UMAP_2D = None
        UMAP_3D = None
        CLUSTER_IDS = None
        CLUSTER_NAMES = None
        GENE_NAMES_LIST = None
        X_CSC = None
        HVG_NAMES = []
        GENE_INDEX = {}
        CSC_BUILDING = False
        CSC_CACHED = False
        CSC_TIME = 0
        N_CELLS = 0
        gc.collect()

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

    # Cluster IDs from cell_type or leiden
    if "cell_type" in adata.obs.columns:
        cats = adata.obs["cell_type"].astype("category")
        CLUSTER_IDS = cats.cat.codes.to_numpy().astype(np.int32)
        CLUSTER_NAMES = cats.cat.categories.tolist()
    elif "leiden" in adata.obs.columns:
        cats = adata.obs["leiden"].astype("category")
        CLUSTER_IDS = cats.cat.codes.to_numpy().astype(np.int32)
        CLUSTER_NAMES = [f"Cluster {c}" for c in cats.cat.categories]
    else:
        # Run leiden clustering
        print("  No cell_type or leiden found, computing clusters...")
        if "connectivities" not in adata.obsp:
            _ensure_neighbors(adata)
        sc.tl.leiden(adata, resolution=1.0)
        cats = adata.obs["leiden"].astype("category")
        CLUSTER_IDS = cats.cat.codes.to_numpy().astype(np.int32)
        CLUSTER_NAMES = [f"Cluster {c}" for c in cats.cat.categories]

    # Gene names (use var index, or feature_name if available)
    if "feature_name" in adata.var.columns:
        # Replace var index with gene names — ensure string type first
        adata.var.index = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()
        GENE_NAMES_LIST = adata.var_names.tolist()
    else:
        adata.var.index = adata.var.index.astype(str)
        adata.var_names_make_unique()
        GENE_NAMES_LIST = adata.var_names.tolist()

    # Compute UMAP(s) — with caching to avoid recomputation
    cache_path = Path(h5ad_path).with_suffix(".umap_cache.npz")
    cached = {}
    if cache_path.exists():
        _progress("Checking UMAP cache...", 20)
        try:
            cached = dict(np.load(cache_path, allow_pickle=False))
            for k in list(cached.keys()):
                if cached[k].shape[0] != N_CELLS:
                    _progress("Cache mismatch, will recompute", 20)
                    cached = {}
                    break
        except Exception as e:
            _progress(f"Cache load failed, recomputing", 20)
            cached = {}

    need_save = False
    neighbors_built = False

    # Decide if we should use subsample-and-project
    use_fast = fast_umap and N_CELLS > fast_umap_subsample

    # Figure out which dims need computing
    dims_needed = []
    for nd in n_dims_list:
        cache_key = f"umap_{nd}d"
        if cache_key in cached:
            _progress(f"Loaded {nd}D UMAP from cache", 25)
            if nd == 2:
                UMAP_2D = cached[cache_key].astype(np.float32)
            else:
                UMAP_3D = cached[cache_key].astype(np.float32)
        else:
            dims_needed.append(nd)

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
            cached[f"umap_{nd}d"] = full_emb
            need_save = True

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
            cached[f"umap_{nd}d"] = umap_arr
            need_save = True

    if need_save:
        _progress("Saving UMAP cache...", 95)
        np.savez(cache_path, **cached)

    ADATA = adata

    # Build fast lookup structures
    GENE_INDEX = {g: i for i, g in enumerate(ADATA.var_names)}

    # Detect HVG names (instant if pre-annotated, fast subsample if not)
    # Used for UI chips only — gene lookups use full CSC once built
    HVG_NAMES = []
    _has_hvg_col = "highly_variable" in ADATA.var.columns
    if _has_hvg_col:
        hvg_mask = ADATA.var["highly_variable"].values.astype(bool)
        HVG_NAMES = [str(ADATA.var_names[i]) for i in np.where(hvg_mask)[0]]
        print(f"  Found {len(HVG_NAMES):,} pre-annotated HVGs")

    if not BACKED and sp.issparse(ADATA.X):
        CSC_BUILDING = True
        csc_cache_path = Path(h5ad_path).with_suffix(".csc_cache.npz")
        def _build_gene_index():
            global HVG_NAMES, X_CSC, CSC_BUILDING, CSC_CACHED, CSC_TIME
            _csc_t0 = time.time()

            # Compute HVGs by variance if not pre-annotated (subsample for speed)
            if not _has_hvg_col:
                print("  Computing top variable genes...")
                t0 = time.time()
                import scipy.sparse as _sp
                X = ADATA.X
                n_sub = min(50000, X.shape[0])
                if n_sub < X.shape[0]:
                    rng = np.random.default_rng(42)
                    X_sub = X[rng.choice(X.shape[0], n_sub, replace=False)]
                else:
                    X_sub = X
                if _sp.issparse(X_sub):
                    mean = np.asarray(X_sub.mean(axis=0)).ravel()
                    sq_mean = np.asarray(X_sub.multiply(X_sub).mean(axis=0)).ravel()
                    var = sq_mean - mean ** 2
                else:
                    var = np.var(np.asarray(X_sub), axis=0)
                del X_sub
                n_hvg = min(2000, len(var))
                hvg_idx = np.argsort(var)[-n_hvg:][::-1]
                HVG_NAMES = [str(ADATA.var_names[i]) for i in hvg_idx]
                print(f"  Computed top {n_hvg:,} variable genes ({time.time()-t0:.1f}s, {n_sub:,}-cell subsample)")

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
                # Save to cache
                try:
                    _sp.save_npz(str(csc_cache_path), X_CSC, compressed=False)
                    cache_mb = os.path.getsize(str(csc_cache_path)) / 1e6
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
    umap_from_cache = len(dims_needed) == 0  # True if ALL dims were cache hits
    load_elapsed = time.time() - t0
    LOAD_SETTINGS.clear()
    LOAD_SETTINGS["backed"] = BACKED
    LOAD_SETTINGS["quickmap"] = use_fast
    LOAD_SETTINGS["quickmap_n"] = fast_umap_subsample if use_fast else 0
    LOAD_SETTINGS["umap_cached"] = umap_from_cache
    LOAD_SETTINGS["load_time"] = round(load_elapsed, 1)

    # Load marker database for cell type annotation
    load_marker_db()

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

    if "connectivities" in adata.obsp:
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


def get_gene_expression(gene_name):
    """Get normalized 0-1 expression for a single gene."""
    import scipy.sparse as sp

    idx = GENE_INDEX.get(gene_name)
    if idx is None:
        return None

    if X_CSC is not None:
        # Fast path: full CSC column slice — O(nnz in column)
        col = X_CSC[:, idx].toarray().ravel().astype(np.float32)
    elif BACKED:
        # Backed mode: read in row batches
        n = ADATA.n_obs
        col = np.empty(n, dtype=np.float32)
        batch = 50000
        for bi in range(0, n, batch):
            be = min(bi + batch, n)
            chunk = ADATA.X[bi:be, idx]
            if sp.issparse(chunk):
                chunk = chunk.toarray().ravel()
            else:
                chunk = np.asarray(chunk).ravel()
            col[bi:be] = chunk
    else:
        # Dense or fallback
        col = ADATA.X[:, idx]
        if sp.issparse(col):
            col = col.toarray().ravel()
        else:
            col = np.asarray(col).ravel()
        col = col.astype(np.float32)

    # Normalize to 0-1
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

    header = {
        "version": 3,
        "n_cells": N_CELLS,
        "n_dims": max(avail),
        "available_dims": avail,
        "n_genes": 0,  # genes loaded on demand
        "n_clusters": len(CLUSTER_NAMES),
        "gene_names": GENE_NAMES_LIST[:500],  # send first 500 for search
        "total_genes": len(GENE_NAMES_LIST),
        "cluster_names": CLUSTER_NAMES,
        "metadata_columns": {},
        "load_settings": dict(LOAD_SETTINGS),
        "hvg_genes": HVG_NAMES[:200],  # top HVGs for quick-access UI
    }

    # Include cached annotations if available
    for annot_key, suffix in [("cached_markers", ".annot_markers.json"), ("cached_llm", ".annot_llm.json")]:
        if LOADED_PATH:
            ap = Path(LOADED_PATH).with_suffix(suffix)
            if ap.exists():
                try:
                    with open(ap) as f:
                        cached = json.load(f)
                    if cached.get("cluster_names") == CLUSTER_NAMES and cached.get("n_cells") == N_CELLS:
                        header[annot_key] = cached.get("results", [])
                except Exception:
                    pass
    header_bytes = json.dumps(header, separators=(",", ":")).encode("utf-8")

    buf = io.BytesIO()
    buf.write(b"SCRN")
    # Pad header to 4-byte alignment
    while len(header_bytes) % 4 != 0:
        header_bytes += b" "
    buf.write(struct.pack("<I", len(header_bytes)))
    buf.write(header_bytes)

    for nd in avail:
        umap = UMAP_2D if nd == 2 else UMAP_3D
        for d in range(nd):
            buf.write(umap[:, d].tobytes())

    buf.write(CLUSTER_IDS.tobytes())

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

    sub = ADATA[indices]
    if BACKED:
        try:
            sub = sub.to_memory()
        except AttributeError:
            sub = sub.copy()
    else:
        sub = sub.copy()
    sub.obs["deg_group"] = labels

    # Normalize if raw counts
    sample = sub.X[:min(100, sub.n_obs), :min(100, sub.n_vars)]
    if sp.issparse(sample):
        sample = sample.toarray()
    sample = np.asarray(sample)
    is_raw = np.allclose(sample, sample.astype(int)) and sample.max() > 5
    if is_raw:
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
    """Export a subset of cells as h5ad, return file path."""
    import tempfile
    sub = ADATA[cell_indices]
    if BACKED:
        try:
            sub = sub.to_memory()
        except AttributeError:
            sub = sub.copy()
    else:
        sub = sub.copy()
    safe_name = "".join(c if c.isalnum() or c in "._-" else "_" for c in group_name)
    fname = f"subset_{safe_name}_{len(cell_indices)}cells.h5ad"
    fpath = os.path.join(tempfile.gettempdir(), fname)
    sub.write_h5ad(fpath)
    print(f"  Export: {len(cell_indices)} cells -> {fpath}")
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
        elif path == "/api/progress":
            self._serve_progress()
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
        elif path == "/api/export":
            self._serve_export(req)
        elif path == "/api/load":
            self._serve_load(req)
        elif path == "/api/annotate":
            self._serve_annotate(req)
        elif path == "/api/annotate-llm":
            self._serve_annotate_llm(req)
        else:
            self.send_error(404)

    def _serve_html(self):
        html_path = Path(__file__).parent / "fleek.html"
        if not html_path.exists():
            self.send_error(500, "fleek.html not found next to server script")
            return
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
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

    def _serve_progress(self):
        p = dict(PROGRESS)
        p["csc_building"] = CSC_BUILDING
        p["csc_cached"] = CSC_CACHED
        p["csc_time"] = CSC_TIME
        p["hvg_ready"] = len(HVG_NAMES) > 0
        if len(HVG_NAMES) > 0:
            p["hvg_genes"] = HVG_NAMES[:200]
        data = json.dumps(p).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(data)

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
            matches = [g for g in GENE_NAMES_LIST if q in g.upper()][:50]
        data = json.dumps(matches).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

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
        """Load a new h5ad file from path (kept for compatibility)."""
        try:
            fpath = req.get("path", "")
            if not fpath or not os.path.exists(fpath):
                raise FileNotFoundError(f"File not found: {fpath}")
            if not fpath.endswith(".h5ad"):
                raise ValueError("File must be .h5ad")
            load_and_prepare(fpath, max_cells=0, n_dims_list=[2, 3])
            result = {"ok": True, "n_cells": N_CELLS, "n_genes": len(GENE_NAMES_LIST)}
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

    def _serve_annotate_llm(self, req):
        """Run cell type annotation using Claude LLM."""
        try:
            if not ANTHROPIC_API_KEY:
                raise ValueError("No API key configured. Start server with --api-key sk-ant-... or set ANTHROPIC_API_KEY env var.")
            top_n = req.get("top_n", 30)
            tissue = req.get("tissue", "")
            result = annotate_clusters_llm(top_n=top_n, tissue_hint=tissue)
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

            # Start processing in background thread
            def _bg():
                try:
                    with PROCESSING_LOCK:
                        load_and_prepare(fpath, max_cells=0, n_dims_list=[2, 3],
                                         fast_umap=fast, fast_umap_subsample=fast_sub,
                                         backed=backed)
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
    parser.add_argument("--api-key", type=str, default="",
                        help="Anthropic API key for Claude cell type annotation (or set ANTHROPIC_API_KEY)")
    args = parser.parse_args()

    if args.api_key:
        global ANTHROPIC_API_KEY
        ANTHROPIC_API_KEY = args.api_key

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
    print(f"{'='*50}\n")

    # Try to open browser
    try:
        import webbrowser
        webbrowser.open(url)
    except:
        pass

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()