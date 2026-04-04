#!/usr/bin/env python3
"""
Download Gene Ontology annotations and build a compact lookup JSON for RNA-FLEEK.

Usage:
    python download_go.py                    # human (default)
    python download_go.py --organism mouse
    python download_go.py --organism all     # all common organisms
    python download_go.py --list             # show available organisms

Output:
    go_human.json, go_mouse.json, etc. — place alongside fleek_server.py

Data sources:
    - GO OBO (ontology structure): http://purl.obolibrary.org/obo/go.obo
    - GAF (gene associations): https://current.geneontology.org/annotations/
"""

import argparse
import gzip
import json
import os
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

# Available organisms and their GAF file names on geneontology.org
ORGANISMS = {
    "human":       {"gaf": "goa_human.gaf.gz",            "label": "Homo sapiens"},
    "mouse":       {"gaf": "mgi.gaf.gz",                  "label": "Mus musculus"},
    "rat":         {"gaf": "rgd.gaf.gz",                  "label": "Rattus norvegicus"},
    "zebrafish":   {"gaf": "zfin.gaf.gz",                 "label": "Danio rerio"},
    "fly":         {"gaf": "fb.gaf.gz",                   "label": "Drosophila melanogaster"},
    "worm":        {"gaf": "wb.gaf.gz",                   "label": "Caenorhabditis elegans"},
    "yeast":       {"gaf": "sgd.gaf.gz",                  "label": "Saccharomyces cerevisiae"},
    "arabidopsis": {"gaf": "tair.gaf.gz",                 "label": "Arabidopsis thaliana"},
    "chicken":     {"gaf": "goa_chicken.gaf.gz",          "label": "Gallus gallus"},
    "pig":         {"gaf": "goa_pig.gaf.gz",              "label": "Sus scrofa"},
    "xenopus":     {"gaf": "xenbase.gaf.gz",              "label": "Xenopus"},
}

GAF_BASE_URL = "https://current.geneontology.org/annotations/"
OBO_URL = "http://purl.obolibrary.org/obo/go.obo"


def download(url, dest, label=""):
    """Download a file with progress indicator."""
    print(f"  Downloading {label or url}...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "RNA-FLEEK/1.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            data = bytearray()
            while True:
                chunk = resp.read(65536)
                if not chunk:
                    break
                data += chunk
                if total > 0:
                    pct = len(data) * 100 // total
                    print(f"\r  {len(data)//1024}KB / {total//1024}KB ({pct}%)", end="", flush=True)
            print()
        with open(dest, "wb") as f:
            f.write(data)
        return True
    except Exception as e:
        print(f"  Error downloading {url}: {e}")
        return False


def parse_obo(obo_path):
    """Parse GO OBO file into {go_id: {name, namespace}}."""
    terms = {}
    current_id = None
    current_name = None
    current_ns = None
    current_alt_ids = []
    is_obsolete = False

    with open(obo_path, "r") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                # Save previous term
                if current_id and current_name and not is_obsolete:
                    entry = {"name": current_name, "namespace": current_ns or ""}
                    terms[current_id] = entry
                    for alt in current_alt_ids:
                        terms[alt] = entry
                current_id = None
                current_name = None
                current_ns = None
                current_alt_ids = []
                is_obsolete = False
            elif line.startswith("id: GO:"):
                current_id = line[4:]
            elif line.startswith("name: "):
                current_name = line[6:]
            elif line.startswith("namespace: "):
                current_ns = line[11:]
            elif line.startswith("alt_id: GO:"):
                current_alt_ids.append(line[8:])
            elif line == "is_obsolete: true":
                is_obsolete = True

    # Save last term
    if current_id and current_name and not is_obsolete:
        entry = {"name": current_name, "namespace": current_ns or ""}
        terms[current_id] = entry
        for alt in current_alt_ids:
            terms[alt] = entry

    return terms


def parse_gaf(gaf_path, go_terms):
    """Parse GAF file into gene-to-GO and GO-to-genes mappings with synonyms.

    GAF format columns (tab-separated):
        0: DB
        1: DB Object ID
        2: DB Object Symbol (gene symbol)
        3: Qualifier
        4: GO ID
        5-9: various
        10: DB Object Synonyms (pipe-separated)
        ...
    """
    gene_to_terms = defaultdict(set)
    term_to_genes = defaultdict(set)
    synonyms = {}  # synonym -> canonical symbol

    opener = gzip.open if gaf_path.endswith(".gz") else open

    with opener(gaf_path, "rt") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue

            symbol = parts[2].strip()
            go_id = parts[4].strip()
            syns_raw = parts[10].strip()

            if not symbol or not go_id:
                continue
            # Skip if GO term is obsolete / not in our OBO
            if go_id not in go_terms:
                continue

            gene_to_terms[symbol].add(go_id)
            term_to_genes[go_id].add(symbol)

            # Parse synonyms
            if syns_raw:
                for syn in syns_raw.split("|"):
                    syn = syn.strip()
                    if syn and syn != symbol:
                        synonyms[syn] = symbol

    return gene_to_terms, term_to_genes, synonyms


def build_go_json(organism_key, obo_path, output_dir):
    """Download GAF for an organism, parse, and build the lookup JSON."""
    org = ORGANISMS[organism_key]
    gaf_url = GAF_BASE_URL + org["gaf"]
    gaf_local = output_dir / org["gaf"]

    print(f"\n{'='*50}")
    print(f"  Building GO database for: {org['label']} ({organism_key})")
    print(f"{'='*50}")

    # Download GAF
    if not gaf_local.exists():
        if not download(gaf_url, str(gaf_local), org["gaf"]):
            return False
    else:
        print(f"  Using cached: {gaf_local.name}")

    # Parse OBO (shared across organisms)
    print("  Parsing GO ontology...")
    go_terms = parse_obo(str(obo_path))
    print(f"  {len(go_terms):,} GO terms loaded")

    # Parse GAF
    print(f"  Parsing {org['gaf']}...")
    gene_to_terms, term_to_genes, synonyms = parse_gaf(str(gaf_local), go_terms)
    print(f"  {len(gene_to_terms):,} genes, {len(term_to_genes):,} terms with annotations, {len(synonyms):,} synonyms")

    # Build compact output
    terms_out = {}
    for go_id, genes in term_to_genes.items():
        t = go_terms.get(go_id)
        if not t:
            continue
        # Abbreviate namespace
        ns = t["namespace"]
        ns_short = "BP" if "biological" in ns else "MF" if "molecular" in ns else "CC" if "cellular" in ns else ns
        terms_out[go_id] = {
            "name": t["name"],
            "ns": ns_short,
            "genes": sorted(genes)
        }

    gene_to_terms_out = {}
    for gene, tids in gene_to_terms.items():
        gene_to_terms_out[gene] = sorted(tids)

    result = {
        "organism": organism_key,
        "label": org["label"],
        "source": org["gaf"],
        "n_terms": len(terms_out),
        "n_genes": len(gene_to_terms_out),
        "n_synonyms": len(synonyms),
        "terms": terms_out,
        "gene_to_terms": gene_to_terms_out,
        "synonyms": synonyms
    }

    out_path = output_dir / f"go_{organism_key}.json"
    with open(out_path, "w") as f:
        json.dump(result, f, separators=(",", ":"))

    size_mb = out_path.stat().st_size / 1024 / 1024
    print(f"  Written: {out_path.name} ({size_mb:.1f} MB)")
    print(f"  {len(terms_out):,} annotated terms, {len(gene_to_terms_out):,} genes, {len(synonyms):,} synonyms")

    # Clean up GAF
    if gaf_local.exists():
        gaf_local.unlink()
        print(f"  Cleaned up: {gaf_local.name}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Download Gene Ontology annotations for RNA-FLEEK",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python download_go.py                     # human (default)
  python download_go.py --organism mouse
  python download_go.py --organism all      # all organisms
  python download_go.py --list              # show available"""
    )
    parser.add_argument("--organism", type=str, default="human",
                        help="Organism to download (default: human). Use 'all' for all.")
    parser.add_argument("--output-dir", type=str, default=str(Path(__file__).parent.parent / "data"),
                        help="Output directory for JSON files (default: ../data/)")
    parser.add_argument("--list", action="store_true",
                        help="List available organisms and exit")
    args = parser.parse_args()

    if args.list:
        print("Available organisms:")
        for key, org in sorted(ORGANISMS.items()):
            print(f"  {key:15s} {org['label']}")
        sys.exit(0)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download OBO (shared)
    obo_path = output_dir / "go.obo"
    if not obo_path.exists():
        if not download(OBO_URL, str(obo_path), "GO ontology (go.obo)"):
            print("Failed to download GO ontology. Exiting.")
            sys.exit(1)
    else:
        print(f"Using cached: {obo_path.name}")

    # Build for requested organism(s)
    if args.organism == "all":
        targets = list(ORGANISMS.keys())
    else:
        if args.organism not in ORGANISMS:
            print(f"Unknown organism: {args.organism}")
            print(f"Available: {', '.join(sorted(ORGANISMS.keys()))}")
            sys.exit(1)
        targets = [args.organism]

    success = 0
    for org in targets:
        if build_go_json(org, obo_path, output_dir):
            success += 1

    # Clean up OBO
    if obo_path.exists():
        obo_path.unlink()
        print(f"\nCleaned up: go.obo")

    print(f"\nDone. Built {success}/{len(targets)} GO databases.")


if __name__ == "__main__":
    main()
