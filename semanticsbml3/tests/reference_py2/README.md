Results of the original semanticSBML (trunk r3402) run under Python 2.7.18 with
libSBML 5.18.0, numpy 1.16.6 and the February 2011 database, on BioModels 61, 64
and 239 (November 2026). Produced by `make_reference.py MODELDIR OUTDIR annotations
check merge match similarity`, which runs unchanged under the original and the port.
`tests/test_against_original.py` compares the port with these files.

`annotations.json` lists element, qualifier, database and identifier; the names the
original program showed for them (taken from its database, partly from KEGG) were
removed.

Running the comparison needs the database of February 2011, which is not in the
repository: `SEMANTICSBML_DATA=<folder> python -m pytest tests/test_against_original.py`.
