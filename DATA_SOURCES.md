# Bundled data and their sources

The repository contains the annotation database that the program needs
(`semanticsbml3/data`, about 100 MB; `semanticsbml3/data_2011`, about 54 MB) and
some test models. These are **not** covered by the GPL of the program code; each
part keeps the terms of its source. The licences below are stated to the best of
our knowledge — **check them before making this repository public**, in
particular KEGG.

| Data | Source | Terms (as we understand them) |
|---|---|---|
| ChEBI names, relations, formulas | EMBL-EBI, `chebi.obo` (release 257 in `data/`, 2011 state in `data_2011/`) | CC BY 4.0 |
| Gene Ontology | geneontology.org, `go-basic.obo` (2026-07-26 / 2011) | CC BY 4.0 |
| Systems Biology Ontology | github.com/EBI-BioModels/SBO | Artistic License 2.0 |
| identifiers.org namespaces | registry.identifiers.org | CC0 |
| MIRIAM registry (2011 export) | EMBL-EBI | predecessor of identifiers.org |
| KEGG Compound and Reaction: identifiers, names, formulas, reaction participants (2011), names of newer entries | KEGG ftp server (2011), KEGG web service name lists (2026) | **KEGG is free for academic use only; redistribution is restricted.** A public repository should probably not contain these tables. |
| NCBI Taxonomy (small subset, 2011) | NCBI | public domain |
| InterPro (2011) | EMBL-EBI | CC0 today |
| SGD cross references (2011) | Saccharomyces Genome Database | CC BY 4.0 today |
| Test models `BIOMD0000000061`, `…64`, `…239` | BioModels | CC0 |

The tables in `data/db` mix these sources (one entry can carry a ChEBI and a KEGG
identifier), so they cannot simply be separated by file. `data/database_sources.txt`
records what the last update brought in. `semanticsbml update-database --no-kegg`
updates without contacting KEGG.

Files in `semanticsbml3/tests/reference_py2/` are outputs of the original
program, kept for comparison.
