#!/usr/bin/env python3
"""
Inspect an h5ad file to discover where counts, metadata, and layers are stored.
Useful for figuring out what's available for pseudo-bulk DEG analysis.

Usage:
    python inspect_h5ad.py path/to/data.h5ad
"""

import sys
import numpy as np

def inspect(path):
    import anndata as ad

    print(f"\nLoading: {path}")
    adata = ad.read_h5ad(path, backed="r")
    print(f"Cells: {adata.n_obs:,}  Genes: {adata.n_vars:,}\n")

    # === X ===
    print("=" * 50)
    print("X (main expression matrix)")
    print("=" * 50)
    print(f"  Shape: {adata.X.shape}")
    print(f"  Type:  {type(adata.X).__name__}")
    try:
        chunk = adata.X[:10, :10]
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        is_int = np.allclose(chunk, chunk.astype(int))
        vmin, vmax = chunk.min(), chunk.max()
        print(f"  Range: {vmin:.4f} — {vmax:.4f}")
        print(f"  Looks like: {'RAW COUNTS (integers)' if is_int else 'NORMALIZED (floats)'}")
        print(f"  Sample:\n    {chunk[0, :8]}")
    except Exception as e:
        print(f"  Error reading X: {e}")

    # === raw ===
    print(f"\n{'=' * 50}")
    print("raw (adata.raw)")
    print("=" * 50)
    if adata.raw is not None:
        print(f"  Shape: {adata.raw.X.shape}")
        print(f"  Genes: {adata.raw.n_vars:,}")
        try:
            chunk = adata.raw.X[:10, :10]
            if hasattr(chunk, "toarray"):
                chunk = chunk.toarray()
            is_int = np.allclose(chunk, chunk.astype(int))
            vmin, vmax = chunk.min(), chunk.max()
            print(f"  Range: {vmin:.4f} — {vmax:.4f}")
            print(f"  Looks like: {'RAW COUNTS (integers)' if is_int else 'NORMALIZED (floats)'}")
            print(f"  Sample:\n    {chunk[0, :8]}")
        except Exception as e:
            print(f"  Error reading raw: {e}")
    else:
        print("  Not present")

    # === layers ===
    print(f"\n{'=' * 50}")
    print("Layers")
    print("=" * 50)
    if len(adata.layers) == 0:
        print("  None")
    else:
        for name in adata.layers:
            try:
                chunk = adata.layers[name][:10, :10]
                if hasattr(chunk, "toarray"):
                    chunk = chunk.toarray()
                is_int = np.allclose(chunk, chunk.astype(int))
                vmin, vmax = chunk.min(), chunk.max()
                label = "RAW COUNTS" if is_int else "NORMALIZED"
                print(f"  '{name}': {label} (range {vmin:.4f}—{vmax:.4f})")
                print(f"    Sample: {chunk[0, :8]}")
            except Exception as e:
                print(f"  '{name}': Error — {e}")

    # === obs (cell metadata) ===
    print(f"\n{'=' * 50}")
    print("Obs columns (cell metadata)")
    print("=" * 50)
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        n_unique = adata.obs[col].nunique()
        tag = ""
        if n_unique <= 30:
            counts = adata.obs[col].value_counts().head(8)
            vals = ", ".join(f"{k}({v})" for k, v in counts.items())
            tag = f" — {vals}"
            if 2 <= n_unique <= 20:
                tag += "  ◄ POTENTIAL CONDITION/SAMPLE"
        print(f"  {col} [{dtype}]: {n_unique} unique{tag}")

    # === var (gene metadata) ===
    print(f"\n{'=' * 50}")
    print("Var columns (gene metadata)")
    print("=" * 50)
    for col in adata.var.columns:
        print(f"  {col} [{adata.var[col].dtype}]")
    print(f"  Gene name examples: {list(adata.var_names[:5])}")

    # === uns ===
    print(f"\n{'=' * 50}")
    print("Uns keys (unstructured metadata)")
    print("=" * 50)
    for key in sorted(adata.uns.keys()):
        val = adata.uns[key]
        vtype = type(val).__name__
        if isinstance(val, str) and len(val) < 100:
            print(f"  {key}: \"{val}\"")
        elif isinstance(val, (int, float, bool)):
            print(f"  {key}: {val}")
        elif isinstance(val, dict):
            print(f"  {key}: dict with {len(val)} keys — {list(val.keys())[:5]}")
        elif isinstance(val, np.ndarray):
            print(f"  {key}: array {val.shape}")
        else:
            print(f"  {key}: {vtype}")

    # === obsm (embeddings) ===
    print(f"\n{'=' * 50}")
    print("Obsm (embeddings)")
    print("=" * 50)
    if len(adata.obsm) == 0:
        print("  None")
    else:
        for name in adata.obsm:
            shape = adata.obsm[name].shape
            print(f"  {name}: {shape}")

    print(f"\n{'=' * 50}")
    print("SUMMARY")
    print("=" * 50)
    # Guess where counts are
    candidates = []
    try:
        chunk = adata.X[:100, :100]
        if hasattr(chunk, "toarray"): chunk = chunk.toarray()
        if np.allclose(chunk, chunk.astype(int)):
            candidates.append("adata.X")
    except: pass
    if adata.raw is not None:
        try:
            chunk = adata.raw.X[:100, :100]
            if hasattr(chunk, "toarray"): chunk = chunk.toarray()
            if np.allclose(chunk, chunk.astype(int)):
                candidates.append("adata.raw.X")
        except: pass
    for name in adata.layers:
        try:
            chunk = adata.layers[name][:100, :100]
            if hasattr(chunk, "toarray"): chunk = chunk.toarray()
            if np.allclose(chunk, chunk.astype(int)):
                candidates.append(f"adata.layers['{name}']")
        except: pass

    if candidates:
        print(f"  Raw counts found in: {', '.join(candidates)}")
    else:
        print("  ⚠ No raw counts detected — all matrices appear normalized")
        print("    Pseudo-bulk DEG requires raw counts. Check if your pipeline")
        print("    stores them elsewhere or if they need to be regenerated.")

    # Guess sample/condition columns
    good_cols = []
    for col in adata.obs.columns:
        n = adata.obs[col].nunique()
        if 2 <= n <= 20:
            good_cols.append((col, n))
    if good_cols:
        print(f"  Likely sample/condition columns:")
        for col, n in good_cols:
            print(f"    {col} ({n} values)")
    else:
        print("  ⚠ No obvious sample/condition columns found (2-20 unique values)")

    print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_h5ad.py path/to/data.h5ad")
        sys.exit(1)
    inspect(sys.argv[1])
