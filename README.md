# semanticSBML 3

Check, annotate, merge and compare SBML models using their semantic
(MIRIAM / identifiers.org) annotations.

This is a port to current Python (3.10 or newer) of the core of semanticSBML,
developed 2005–2011 at the Max Planck Institute for Molecular Genetics and
Humboldt University Berlin (Falko Krause, Marvin Schulz, Jannis Uhlendorf,
Wolfram Liebermeister, Anselm Helbig, Timo Lubitz, Stefan Stiller). It is based
on the last state of the SourceForge repository (trunk, revision 3402, version 1.1.0).

Licence: GPL 3 (see `COPYING`).

## Installation

```
cd semanticsbml3
python3 -m venv .venv
.venv/bin/pip install -e .
```

This installs the command `semanticsbml` into `.venv/bin`. The pictures of
`semanticsbml view` need [Graphviz](https://graphviz.org) (`sudo apt install graphviz`).

## Usage

```
semanticsbml check        MODEL.xml            # report missing/duplicate annotations, atom balance
semanticsbml view         MODEL.xml -o m.svg   # picture of the reaction network (.svg .png .pdf .dot)
semanticsbml annotations  MODEL.xml            # table of all annotations (--missing, --tsv)
semanticsbml search       "glucose"            # look up a name in the annotation database
semanticsbml annotate     MODEL.xml -o OUT.xml --auto
semanticsbml annotate     MODEL.xml -o OUT.xml --add GLCi ChEBI CHEBI:17234 is
semanticsbml annotate     MODEL.xml -o OUT.xml --remove GLCi ChEBI CHEBI:17234
semanticsbml annotate     MODEL.xml -o OUT.xml --copy-from OTHER.xml
semanticsbml merge        A.xml B.xml -o MERGED.xml
semanticsbml match        A.xml B.xml          # corresponding elements, by similarity
semanticsbml similarity   A.xml B.xml C.xml    # table of model similarities
semanticsbml similarity   --query A.xml --dir FOLDER   # rank the models of a folder
```

`semanticsbml COMMAND --help` lists the options of each command.

**Annotate.** `--auto` annotates every element without annotation whose name
occurs exactly once in the database. `search` shows what the database offers for a
name; `--add ELEMENT DATABASE IDENTIFIER [QUALIFIER]` adds one annotation by hand.

**Merge.** Elements whose annotations agree are joined; everything else is taken
over from each model. Where joined elements differ (values, formulas, units,
reaction participants) the model named first wins. A report lists the joined
elements and the differences, and says whether the result is valid SBML. Models
that contradict each other (e.g. a species set by a rule in one model and by
reactions in the other) give an invalid result that has to be repaired by hand;
changing the order of the models often helps.

**Match** finds corresponding elements more tolerantly than merge: it also uses
relations between database entries (D-glucose *is a* glucose) and the
neighbourhood of an element in the network. It only reports, it does not merge.

**Similarity** compares whole models (0 = nothing in common, 1 = same annotations).

## The annotation database

The folder `data/` holds the annotation database: names and relations from ChEBI,
KEGG Compound, KEGG Reaction, Gene Ontology, NCBI Taxonomy (small), SBO, InterPro
and SGD, and the list of databases known to MIRIAM / identifiers.org. It is needed
for names, for `search`, `--auto`, `check` (atom balance), `match` and `similarity`.

```
semanticsbml database           # which versions are in use
semanticsbml update-database    # fetch current ChEBI, Gene Ontology, SBO, KEGG names, identifiers.org
```

`update-database` downloads about 100 MB into `data/downloads` and takes a few
minutes. It keeps the existing tables and replaces the names and relations of
ChEBI, Gene Ontology and SBO by the current ones; new entries are added. From
KEGG only the two lists of names are fetched (KEGG is free for academic use only;
`--no-kegg` leaves it out). Taxonomy, InterPro and SGD stay as they were in 2011.
What an update cannot bring back: chemical formulas and reaction equations of KEGG
entries added after 2011 are taken from ChEBI and from the reaction names where possible.

`data_2011/` is the last database of the original tool (February 2011), untouched.
The tests use it, and `semanticsbml --data data_2011 ...` works with it.
Another data folder can also be chosen with the environment variable `SEMANTICSBML_DATA`.

## Tests

```
.venv/bin/pip install pytest
.venv/bin/python -m pytest
```

## Layout

| folder | content |
|---|---|
| `semanticSBML/` | check, annotate, merge, view, similarity of elements; `cli.py` is the command line tool |
| `libSBAnnotation/` | annotation database, MIRIAM resources, SBO |
| `libModelSimilarity/` | similarity and clustering of whole models |
| `data/` | annotation database (updated) and resource files |
| `data_2011/` | annotation database as of February 2011 |
| `tests/` | tests and example models (BioModels 61, 64, 239) |

See `CHANGES.md` for what differs from the original program.
