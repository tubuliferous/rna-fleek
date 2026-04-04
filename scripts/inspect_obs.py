#!/usr/bin/env python3
"""Quick inspector for obs metadata columns in an h5ad file."""
import sys
import anndata as ad
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python inspect_obs.py file.h5ad")
    sys.exit(1)

path = sys.argv[1]
print(f"\nLoading {path} (obs only)...")
adata = ad.read_h5ad(path, backed="r")  # backed mode — don't load expression matrix

n_cells = adata.n_obs
print(f"{n_cells:,} cells  ·  {adata.n_vars:,} genes\n")

# Keywords that suggest experimental metadata
KEYWORDS = {"condition", "sample", "treatment", "batch", "timepoint", "time",
            "donor", "patient", "subject", "replicate", "group", "stim",
            "stimulation", "genotype", "disease", "status", "phase", "sex",
            "age", "tissue", "organ", "library", "lane", "pool", "site"}

categorical, numerical, other = [], [], []

for col in adata.obs.columns:
    series = adata.obs[col]
    dtype = series.dtype
    n_unique = series.nunique()
    flag = "  ★" if any(kw in col.lower() for kw in KEYWORDS) else ""

    if hasattr(dtype, "categories") or dtype == object or str(dtype) == "category":
        cats = list(series.astype("category").cat.categories[:8])
        preview = ", ".join(str(c) for c in cats)
        if n_unique > 8:
            preview += f", … (+{n_unique - 8} more)"
        categorical.append((col, n_unique, preview, flag))
    elif np.issubdtype(dtype, np.number):
        mn, mx = series.min(), series.max()
        numerical.append((col, dtype, f"{mn:.3g} – {mx:.3g}", n_unique, flag))
    else:
        other.append((col, str(dtype), n_unique, flag))

def header(title):
    print(f"\n{'─'*60}")
    print(f"  {title}")
    print(f"{'─'*60}")

header("CATEGORICAL / STRING columns  (★ = likely experimental metadata)")
if categorical:
    for col, n_u, preview, flag in categorical:
        print(f"  {col}{flag}")
        print(f"    {n_u} unique values: {preview}")
else:
    print("  (none)")

header("NUMERICAL columns")
if numerical:
    for col, dtype, rng, n_u, flag in numerical:
        print(f"  {col}{flag}  [{dtype}]  range: {rng}  ({n_u} unique)")
else:
    print("  (none)")

if other:
    header("OTHER columns")
    for col, dtype, n_u, flag in other:
        print(f"  {col}{flag}  [{dtype}]  {n_u} unique")

print(f"\n{'─'*60}")
print(f"  Total obs columns: {len(adata.obs.columns)}")
print(f"{'─'*60}\n")
