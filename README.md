# semanticSBML 3

Check, annotate, merge and compare SBML models using their semantic
(MIRIAM / identifiers.org) annotations — from the command line or in the browser.

This is a port to current Python of the core of **semanticSBML**, developed
2005–2011 at the Max Planck Institute for Molecular Genetics and Humboldt
University Berlin (Falko Krause, Marvin Schulz, Jannis Uhlendorf, Wolfram
Liebermeister, Anselm Helbig, Timo Lubitz, Stefan Stiller). The original was
written for Python 2 and no longer runs on current systems. The port is based on
the last state of the SourceForge repository (trunk, revision 3402, version 1.1.0)
and was made in November 2026 with the help of an AI coding assistant (Claude Code).

## What it does

| | |
|---|---|
| **Annotate** | list, search, add and remove annotations; annotate elements automatically by name |
| **Check** | missing and duplicate annotations, atom balance of reactions |
| **Merge** | join models: elements with the same annotations become one, with a report of all differences |
| **Match** | find corresponding elements of two models, using relations between database entries |
| **Similarity** | compare whole models by their annotations |
| **View** | picture of the reaction network |

Not included: parameter balancing and SBtab (developed further elsewhere), the
old graphical and web interfaces, building models from KEGG.

## Quick start

Needs Python 3.10 or newer; for network pictures also [Graphviz](https://graphviz.org).

```
cd semanticsbml3
python3 -m venv .venv
.venv/bin/pip install -e .

.venv/bin/semanticsbml gui                          # browser interface
.venv/bin/semanticsbml merge A.xml B.xml -o AB.xml  # command line
.venv/bin/semanticsbml --help
```

## Where to find what

| | |
|---|---|
| [`semanticsbml3/README.md`](semanticsbml3/README.md) | installation, all commands, browser interface, annotation database |
| [`semanticsbml3/CHANGES.md`](semanticsbml3/CHANGES.md) | every difference from the original program, and the comparison with it |
| [`DATA_SOURCES.md`](DATA_SOURCES.md) | where the bundled annotation data come from, and their licences |
| [`LICENSE`](LICENSE) | GNU General Public License, version 3 (the licence of the original) |

## State

Alpha. Annotation, merging and similarity were compared with the original program
run under Python 2.7 on three BioModels models (results in
`semanticsbml3/tests/reference_py2/`); the differences are intended corrections
and are listed in `CHANGES.md`. 56 automated tests (`.venv/bin/python -m pytest`
inside `semanticsbml3`, after `.venv/bin/pip install pytest`).

## References

- Krause F, Uhlendorf J, Lubitz T, Schulz M, Klipp E, Liebermeister W (2010).
  Annotation and merging of SBML models with semanticSBML. *Bioinformatics* 26(3):421–422.
- Schulz M, Krause F, Le Novère N, Klipp E, Liebermeister W (2011).
  Retrieval, alignment, and clustering of computational models based on semantic
  annotations. *Molecular Systems Biology* 7:512.
