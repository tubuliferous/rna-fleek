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
                     fast_umap_subsample=50000, backed="auto"):
    """Load h5ad, subsample if needed, compute UMAP(s).
    
    backed: "auto" (use file size vs RAM), "on" (force backed), "off" (force in-memory)
    """
    global ADATA, UMAP_2D, UMAP_3D, CLUSTER_IDS, CLUSTER_NAMES, GENE_NAMES_LIST, N_CELLS, BACKED

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
        N_CELLS = 0
        gc.collect()

    _progress("Loading h5ad file...", 5)
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

    if gene_name not in ADATA.var_names:
        return None

    idx = list(ADATA.var_names).index(gene_name)

    if BACKED:
        # Backed mode: read column in chunks to avoid full-matrix scan
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
    }
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
    """Run DEG between two cell groups with Cohen's d."""
    import scanpy as sc
    import scipy.sparse as sp

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
            sub = sub.copy()  # older anndata
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

    # Cohen's d
    print(f"  DEG: computing Cohen's d...")
    X_a = sub.X[:na]
    X_b = sub.X[na:]
    if sp.issparse(X_a):
        X_a = X_a.toarray()
        X_b = X_b.toarray()
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
    print(f"  DEG: done in {elapsed:.1f}s — {len(clean):,} genes, {n_sig:,} significant")

    return {
        "genes": clean, "test": method, "n_a": na, "n_b": nb,
        "n_total_genes": n_all, "elapsed": round(elapsed, 1),
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
        data = json.dumps(PROGRESS).encode("utf-8")
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
            with open(fpath, "rb") as f:
                data = f.read()
            os.remove(fpath)
            self.send_response(200)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Disposition", f'attachment; filename="{fname}"')
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
            backed = self.headers.get("X-Backed", "auto")  # "auto", "on", "off"

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
    parser.add_argument("--backed", type=str, default="auto", choices=["auto", "on", "off"],
                        help="Backed mode: auto (detect), on (force disk), off (force RAM)")
    args = parser.parse_args()

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
