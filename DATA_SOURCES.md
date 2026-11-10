# Data sources

## Built on your computer (not in the repository)

`semanticsbml build-database` downloads the following files into
`semanticsbml3/data/downloads` and builds the annotation database in
`semanticsbml3/data/db`. Neither folder is part of the repository; each user
fetches the data from the providers. The terms below are stated to the best of
our knowledge; check them at the source.

| Data | Source | Terms |
|---|---|---|
| ChEBI (names, synonyms, relations, formulas, cross references) | EMBL-EBI, `https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz` | CC BY 4.0 |
| Gene Ontology | `https://purl.obolibrary.org/obo/go/go-basic.obo` | CC BY 4.0 |
| Systems Biology Ontology | `https://github.com/EBI-BioModels/SBO` | Artistic License 2.0 |
| identifiers.org registry | `https://registry.api.identifiers.org/resolutionApi/getResolverDataset` | CC0 |

KEGG identifiers occur in the database only as cross references published by ChEBI
and the Gene Ontology. Nothing is downloaded from KEGG unless `--with-kegg` is given;
KEGG data are free for academic use only and may not be redistributed.

## In the repository

| Data | Source | Terms |
|---|---|---|
| `data/ncbi_taxonomy_small_*.dmp`: 1500 organisms of the NCBI taxonomy | extract made by the original semanticSBML project (2009) | NCBI data are public domain |
| `data/MIRIAM_resources.xml`: MIRIAM registry, export of February 2011 | EMBL-EBI | predecessor of identifiers.org |
| `data/listofresources.xml`, `libsbmlresources.xml`, `elementscore.txt`, `qualifiers.csv`, `config` | original semanticSBML project | GPL 3, like the code |
| `data/bimatrix.csv`, `simmatrix.csv`, `biomodels_probabilities`: statistics of the annotations in BioModels (2010–2011), used only by library functions not reachable from the command line | original semanticSBML project | GPL 3 |
| Test models `tests/models/*.xml` | BioModels | CC0 |
| `tests/reference_py2/`: outputs of the original program, kept for comparison | this port | GPL 3 |
