#!/usr/bin/env python3
"""
Download a tissue slice from CELLxGENE Census and save as a local h5ad file.
Run this ONCE, then use preprocess_census.py --h5ad on the saved file.

Usage:
    # Start small to test (fast):
    python download_census.py --tissue liver

    # Big immune atlas:
    python download_census.py --tissue blood

    # Cap the download to a manageable size:
    python download_census.py --tissue blood --max-cells 2000000

Then preprocess from local file (fast, no network):
    python preprocess_census.py --h5ad liver_census.h5ad --max-cells 500000
"""

import argparse
import os
import sys
import time
import threading


def monitor_memory(stop_event, interval=10):
    """Background thread that prints memory usage periodically."""
    try:
        import psutil
        proc = psutil.Process()
        while not stop_event.is_set():
            mem = proc.memory_info().rss / 1e9
            print(f"    [monitor] RAM usage: {mem:.1f} GB", flush=True)
            stop_event.wait(interval)
    except ImportError:
        # psutil not installed, fall back to simple timer
        t0 = time.time()
        while not stop_event.is_set():
            elapsed = time.time() - t0
            mins = int(elapsed // 60)
            secs = int(elapsed % 60)
            print(f"    [monitor] Still downloading... {mins}m {secs}s elapsed", flush=True)
            stop_event.wait(interval)


def main():
    parser = argparse.ArgumentParser(description="Download Census tissue data as h5ad")
    parser.add_argument("--tissue", type=str, required=True,
                        help="tissue_general value (blood, lung, brain, liver, heart, kidney, etc.)")
    parser.add_argument("--max-cells", type=int, default=0,
                        help="Cap cell count (0 = download all). Recommended for blood/brain.")
    parser.add_argument("--output", type=str, default=None,
                        help="Output h5ad path (default: <tissue>_census.h5ad)")
    parser.add_argument("--census-version", type=str, default="stable")
    args = parser.parse_args()

    out = args.output or f"{args.tissue}_census.h5ad"

    import cellxgene_census
    import numpy as np

    print(f"Opening Census (version: {args.census_version})...")
    census = cellxgene_census.open_soma(census_version=args.census_version)

    obs_filter = f"tissue_general == '{args.tissue}' and is_primary_data == True"

    # Count first
    print(f"Counting cells for tissue_general='{args.tissue}'...")
    obs_df = cellxgene_census.get_obs(
        census, "Homo sapiens", value_filter=obs_filter,
        column_names=["soma_joinid"],
    )
    total = len(obs_df)
    print(f"Found {total:,} cells")

    # Estimate file size (very rough: ~3-5 KB per cell for sparse h5ad)
    est_gb = total * 4e-6  # ~4 KB/cell is a reasonable average
    print(f"Estimated h5ad size: ~{est_gb:.0f} GB")
    print(f"Estimated RAM needed: ~{est_gb * 1.5:.0f} GB (1.5x file size)")

    # Subsample if requested
    download_n = total
    obs_coords = None
    if args.max_cells > 0 and total > args.max_cells:
        download_n = args.max_cells
        print(f"Capping download to {download_n:,} cells (use --max-cells 0 for all)")
        rng = np.random.default_rng(42)
        keep = rng.choice(total, size=download_n, replace=False)
        keep.sort()
        obs_coords = obs_df.iloc[keep]["soma_joinid"].to_numpy()
    elif total > 3_000_000:
        print(f"\n⚠️  That's a LOT of cells. This will use ~{est_gb * 1.5:.0f} GB RAM.")
        print(f"    Consider: python download_census.py --tissue {args.tissue} --max-cells 2000000")
        print(f"    Proceeding with full download in 5 seconds... (Ctrl+C to cancel)\n")
        time.sleep(5)

    print(f"Downloading {download_n:,} cells × ~60k genes...")
    print(f"This will take a while. Progress updates every 30 seconds.\n")

    # Start background monitor
    stop_event = threading.Event()
    monitor = threading.Thread(target=monitor_memory, args=(stop_event, 30), daemon=True)
    monitor.start()

    t0 = time.time()
    try:
        if obs_coords is not None:
            adata = cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_coords=obs_coords,
                obs_column_names=[
                    "cell_type", "tissue", "disease", "assay",
                    "donor_id", "sex", "development_stage",
                    "self_reported_ethnicity", "suspension_type",
                ],
            )
        else:
            adata = cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_value_filter=obs_filter,
                obs_column_names=[
                    "cell_type", "tissue", "disease", "assay",
                    "donor_id", "sex", "development_stage",
                    "self_reported_ethnicity", "suspension_type",
                ],
            )
    finally:
        stop_event.set()
        monitor.join(timeout=2)

    elapsed = time.time() - t0
    census.close()

    print(f"\nDownloaded {adata.n_obs:,} cells × {adata.n_vars:,} genes in {elapsed/60:.1f} min")
    print(f"Saving to {out}...")

    t1 = time.time()
    adata.write_h5ad(out)
    save_time = time.time() - t1

    size_mb = os.path.getsize(out) / 1e6
    print(f"Saved {size_mb/1000:.1f} GB → {out} ({save_time:.0f}s to write)")
    print(f"\nNow preprocess with:")
    print(f"  python preprocess_census.py --h5ad {out} --max-cells 1000000")


if __name__ == "__main__":
    main()
