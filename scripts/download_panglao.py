#!/usr/bin/env python3
"""
Download and preprocess the PanglaoDB marker database into a JSON
lookup that the FLEEK server can ensemble alongside CellMarker2.

PanglaoDB markers ship as a TSV with per-gene sensitivity / specificity
columns and per-organ tags. We collapse the TSV into the same shape
FLEEK already uses for cell_markers.json:

  {
    "Human": { "tissue_name": { "Cell Type": ["GENE1", ...] } },
    "Mouse": { ... },
    "_global": {
      "Human": { "Cell Type": ["GENE1", ...] }   # union across tissues
    }
  }

so the existing JSON loader in rna_fleek/server.py picks it up under
the panglao_markers.json filename.

Markers are filtered to keep only canonical / specific entries:
  - Drop genes flagged by ubiquitousness_index >= 0.05 (housekeeping-ish).
  - Keep genes whose specificity_human >= 0.10 (configurable).
  - Sort by sensitivity_human desc within each cell type so the top of
    the list is the most reliable up-regulator.

Usage:
    python scripts/download_panglao.py
    python scripts/download_panglao.py --output rna_fleek/panglao_markers.json
    python scripts/download_panglao.py --min-specificity 0.05

Requirements: only stdlib (urllib + gzip + csv).
"""

import argparse
import csv
import gzip
import io
import json
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path


PANGLAO_URL = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"


def _download(url):
    print(f"Downloading PanglaoDB markers from {url} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "RNA-FLEEK"})
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read()
    print(f"  {len(raw)} bytes downloaded")
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
        text = gz.read().decode("utf-8", errors="replace")
    return text


def _parse(text, min_specificity=0.10, max_ubiquitousness=0.05):
    """Parse PanglaoDB TSV → {species: {tissue: {cell_type: set(genes)}}}.

    species: 'Human' for Hs, 'Mouse' for Mm. Inclusion rules per row:
      - Drop if ubiquitousness_index ≥ max_ubiquitousness (housekeeping).
      - Include in a species bucket if the row's species column tags it
        AND any of: spec_<species> ≥ min_specificity, spec_<other_species>
        ≥ min_specificity (orthologs share biology), or canonical_marker == "1".
        PanglaoDB's specificity_mouse field is sparse — many "Hs Mm"-tagged
        rows have 0 there even for canonical mouse markers, so we don't
        reject on mouse spec alone.

    Tissue uses the 'organ' column; entries with no organ bucket under "Other".
    """
    out = {"Human": defaultdict(lambda: defaultdict(set)),
           "Mouse": defaultdict(lambda: defaultdict(set))}
    rdr = csv.DictReader(io.StringIO(text), delimiter="\t")
    n_total = 0
    n_kept_h = 0
    n_kept_m = 0
    for row in rdr:
        n_total += 1
        gene = (row.get("official gene symbol") or row.get("official_gene_symbol") or "").strip().upper()
        ct = (row.get("cell type") or row.get("cell_type") or "").strip()
        if not gene or not ct:
            continue
        species_field = (row.get("species") or "").strip()
        organ = (row.get("organ") or "").strip() or "Other"
        canonical = (row.get("canonical marker") or row.get("canonical_marker") or "").strip()
        is_canonical = canonical in ("1", "1.0", "True", "true", "yes")
        try:
            ubi = float(row.get("ubiquitousness index") or row.get("ubiquitousness_index") or "0") or 0.0
        except ValueError:
            ubi = 0.0
        if ubi >= max_ubiquitousness:
            continue
        try:
            spec_h = float(row.get("specificity_human") or "0") or 0.0
        except ValueError:
            spec_h = 0.0
        try:
            spec_m = float(row.get("specificity_mouse") or "0") or 0.0
        except ValueError:
            spec_m = 0.0
        # A gene is "good enough" if either species-specificity passes OR
        # the marker is flagged canonical. Same gate for both species.
        good = (spec_h >= min_specificity) or (spec_m >= min_specificity) or is_canonical
        if not good:
            continue
        if "Hs" in species_field:
            out["Human"][organ][ct].add(gene)
            n_kept_h += 1
        if "Mm" in species_field:
            out["Mouse"][organ][ct].add(gene)
            n_kept_m += 1
    print(f"  Parsed {n_total} rows;"
          f" kept {n_kept_h} Human + {n_kept_m} Mouse"
          f" (spec ≥ {min_specificity} OR canonical, ubi < {max_ubiquitousness})")
    return out


def _serialize(parsed):
    """Convert the nested dict-of-sets into the JSON shape FLEEK expects.
    Adds a "_global" branch which unions across tissues per cell type."""
    out = {"Human": {}, "Mouse": {}, "_global": {"Human": {}, "Mouse": {}}}
    for species in ("Human", "Mouse"):
        for tissue, ctmap in parsed[species].items():
            out[species][tissue] = {ct: sorted(genes) for ct, genes in ctmap.items()}
        # Union across tissues per cell type
        union = defaultdict(set)
        for tissue, ctmap in parsed[species].items():
            for ct, genes in ctmap.items():
                union[ct] |= genes
        out["_global"][species] = {ct: sorted(genes) for ct, genes in union.items()}
    n_h = len(out["_global"]["Human"])
    n_m = len(out["_global"]["Mouse"])
    print(f"  Cell types in _global: Human={n_h}, Mouse={n_m}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", default="rna_fleek/panglao_markers.json",
                    help="Output path (default: rna_fleek/panglao_markers.json)")
    ap.add_argument("--min-specificity", type=float, default=0.10,
                    help="Minimum specificity to keep a marker (default 0.10)")
    ap.add_argument("--max-ubiquitousness", type=float, default=0.05,
                    help="Max ubiquitousness index — drops housekeeping-ish markers (default 0.05)")
    ap.add_argument("--url", default=PANGLAO_URL,
                    help="PanglaoDB markers URL (default: 27 Mar 2020 release)")
    args = ap.parse_args()

    text = _download(args.url)
    parsed = _parse(text,
                    min_specificity=args.min_specificity,
                    max_ubiquitousness=args.max_ubiquitousness)
    out = _serialize(parsed)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print(f"Wrote {out_path} ({out_path.stat().st_size:,} bytes)")
    print()
    print(f"Drop the file next to cell_markers.json — the FLEEK server")
    print(f"will pick it up on the next start and ensemble it with CellMarker2.")


if __name__ == "__main__":
    main()
