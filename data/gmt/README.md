# Pathway / gene-set databases for the Pathways panel

The Pathways panel runs over-representation analysis (ORA) against gene-set
collections in [GMT format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).

GO Biological Process / Molecular Function / Cellular Component already ship
with FLEEK (loaded from `data/go_<organism>.json`). Other databases are
opt-in: drop a `.gmt` file in `data/gmt/<organism>/` and the server will
discover it on dataset load.

## Directory layout

```
data/
  gmt/
    human/
      hallmark.gmt
      reactome.gmt
    mouse/
      hallmark.gmt
      reactome.gmt
```

The file stem becomes the database id in the panel dropdown
(`hallmark` → "Hallmark", `reactome` → "Reactome", etc.).

## Recommended sources (license-clean)

- **MSigDB Hallmark (50 sets, CC-BY-4.0)**
  Download from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
  Look for "h.all.v*.symbols.gmt" (human) and the matching mouse file.
  Rename to `hallmark.gmt` and place in `data/gmt/human/` or `data/gmt/mouse/`.

- **Reactome (~2500 sets, CC-BY-4.0)**
  Download from: https://reactome.org/download-data
  Look for "ReactomePathways.gmt" (organism-specific symbol-based).

- **MSigDB C2 / C5 / C7** are also CC-BY-4.0 and follow the same workflow;
  drop them in with whatever filename you like.

## What about KEGG?

KEGG's pathway gene sets are released under a license that restricts
redistribution. FLEEK doesn't bundle them. Academic users can obtain a
KEGG GMT through their institution's MSigDB access (KEGG_LEGACY collection)
or build one from the KEGG REST API for personal use.

## File format

One pathway per line, tab-separated:

```
PATHWAY_ID<TAB>DESCRIPTION_OR_URL<TAB>GENE1<TAB>GENE2<TAB>...
```

If column 2 starts with `http://` or `https://` the panel treats it as a
"View source" link. Otherwise it's shown as a free-form description in the
expanded row detail.

Gene name matching is case-insensitive against the dataset's gene list, so
a human-style GMT (uppercase symbols) works on a mouse dataset (Title-case
symbols) — the server's loader handles the canonical form mapping.
