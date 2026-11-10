# Differences from the original semanticSBML (trunk r3402, Python 2)

The aim of the port was to keep the algorithms as they were. This file lists
everything that is not a pure translation to Python 3.

## Not ported

- the graphical interface (PyQt4) and the web interface (web2py); replaced by a new browser interface
- parameter balancing, kineticizer, SBtab (developed further elsewhere)
- building models from KEGG (`kegg2sbml`): the KEGG web service it used was shut down
- libTI / libTI2 / libPyless (simulation, identifiability)
- the original way of building the annotation database: most download addresses
  no longer exist. The code (`libSBAnnotation/buildDB.py`) is kept (its taxonomy
  readers are used); the database is now built by `semanticsbml build-database`, see below.
  A missing database gives an error message instead of an attempt to download.

## New

- one command `semanticsbml` with sub-commands (`semanticSBML/cli.py`)
- browser interface `semanticsbml gui` (`semanticSBML/webgui.py`, `webgui.html`): a new
  page, not a port of the PyQt or web2py interfaces; it needs no additional software.
  The server listens on 127.0.0.1 only and refuses requests from other web pages.
- annotations in the identifiers.org form (`http://identifiers.org/chebi/CHEBI:17234`,
  `https://identifiers.org/CHEBI:17234`) are read; new annotations are written
  in this form (set `uri_style = urn` in section `[annotate]` of `data/config`
  for the old `urn:miriam:` form)
- annotations from databases unknown to the 2011 MIRIAM table are kept instead
  of being rejected
- qualifiers added to SBML since 2011: `hasTaxon`, `isInstanceOf`, `hasInstance`
  (with low scores in `data/elementscore.txt` and `data/qualifiers.csv`)
- data are read from `./data` (or `$SEMANTICSBML_DATA`), not from `~/.semanticSBML`;
  nothing is written to the home directory
- `libModelSimilarity/api.py`: simple functions for model similarity
- `semanticsbml build-database` (`libSBAnnotation/update.py`) builds the annotation
  database from current public resources, without KEGG:
  - ChEBI, Gene Ontology and SBO from their current ontology files, the identifiers.org
    registry, and the small NCBI taxonomy extract of the original project (`data/`).
    Of ChEBI the manually curated entries ("3 star"), those linked to a KEGG compound
    and all their parents are taken (`--all-chebi`: all).
  - KEGG compound and reaction identifiers are known only as links given by ChEBI and
    the Gene Ontology; there are no KEGG names, formulas or reaction equations.
    InterPro and SGD, part of the 2011 database, are not included.
  - Former identifiers (`alt_id`, `replaced_by`) are linked to their successor the way the
    original did for ChEBI ("successor is_a former entry"); now also for GO.
  - `--update` applies the same files to an existing database instead (names and
    relations of ChEBI, GO and SBO are replaced, new entries added, links made earlier
    are never rearranged); `--with-kegg` additionally takes the KEGG name lists.
    This is how the database of 2011 was brought up to date in November 2026.
  - The database is no longer distributed with the program.
  - Formulas from ChEBI are kept in a new table (`formula`); KEGG formulas of old
    databases (`kegg_formula`) are still used.
- database tables build their indexes when first needed (start-up 4 s instead of 15 s
  with the updated database)
- tests (`tests/`)

## Changes needed for current libSBML

- libSBML objects can be used as dictionary keys again (`semanticSBML/libsbml_compat.py`)
- renaming a rule or initial assignment: libSBML now ignores `setId` for them;
  the variable is set instead. (Without this, rules were lost when merging.)
- removing the last annotation of a qualifier left an empty group, which makes
  current libSBML crash on writing; the groups are rebuilt instead

## Corrections of the original behaviour

Merge
- A better matching element now replaces a weaker one. Before, the pair was
  dissolved and the weaker match was restored later (e.g. ADP was paired with a
  pool "ATP + ADP" although the other model contained ADP).
- "First model wins": all merged elements are created first, then formulas are
  translated. Before, formulas could refer to identifiers that were replaced later.
- Unit definitions are kept by identifier. Before, two definitions with the same
  content but different identifiers were reduced to one, leaving references to
  the other one dangling (invalid SBML).
- Function definitions used by the kinetic law of a merged reaction are carried
  over to the merged model.
- In `updateMathElementIDs` an assignment `math.setName = ...` was meant to be a call.
- The state kept between merges (identifiers in use) is reset for every merge.
- The order of elements in the output is the same in every run.
- The command line merge script of the original (`semanticsbml-automerge`) passed
  the wrong kind of object to the merger and could not run; `semanticsbml merge` replaces it.

Check
- The atom balance check took formulas only from KEGG Compound annotations; ChEBI
  annotations are used as well now (with the database of 2011 this gives the same
  results for the test models; with a database built today many more reactions can
  be checked).

Search
- The search for a name answered from a cache of names displayed before, so the
  result depended on what had been looked up earlier in the same run (e.g. a name
  with two database entries was reported with one). The cache is no longer used.
  The merge code relied on this cache when comparing two annotations by name; the
  decision it took in effect (same identifier) is now written out and unchanged.
- An error in the database code that only showed when entries are added
  (`AnnotationObject.items`) was fixed.

View
- `semanticsbml-view` used an undefined variable and wrong arguments; fixed.
- no temporary file `test.dot` is written to the working directory.

Similarity
- `2 / (height1 + height2)` in `similarity.py` was an integer division under
  Python 2 (result 0 or 1); it is now a true division, as the formula intends.
  **Element similarities for related (not identical) database entries therefore
  differ from values computed with the original program.**
  Set `semanticSBML.similarity.HEIGHT_INTEGER_DIVISION = True` to get the old numbers.
- The height of a database entry is computed breadth first instead of following
  every path (same values; the recursive version was very slow under Python 3).
- The walk through the "is_a" hierarchy for the information content did not
  terminate when the hierarchy contained circles (the 2011 database does); guarded.
- removed SciPy functions replaced by their NumPy equivalents.

## Known limitations

- Merging converts all models to SBML Level 2 Version 4, as the original did.
  Events, constraints, species types, compartment types and algebraic rules are ignored.
- `match` with similarity propagation (`GlobalSigmaSimProp`) is still in the code
  but too slow for use (dense matrices of size (n·m)²); the command line uses
  annotation propagation, the default of the original.

## Comparison with the original program

In November 2026 the original code (trunk r3402) was run under Python 2.7.18 with
libSBML 5.18.0 and the same 2011 database, on BioModels 61, 64 and 239. (libSBML
5.18 is from 2019, not the version of 2011, so this is the original code but not
exactly the original installation.) The results are stored in
`tests/reference_py2/` and compared by `tests/test_against_original.py`.

- Reading annotations: identical (213 + 116 annotations). Where an entry has several
  names or a reaction several products, the name shown can differ in order or choice
  (e.g. "ATPase activity, coupled" / "…, uncoupled"); the original picked arbitrarily.
- Check: identical numbers of errors and warnings in every test.
- Merge of 61 and 64: both find 27 groups; 24 are the same. The three that differ:
  ATP and ADP were paired with the pooled species "sum of AXP" and "High energy
  phosphates" in the original (correction listed above), and extracellular glucose
  of model 64 was paired with "GlcX0" instead of "GlcX" of model 61 (both carry the
  same annotation; a tie). The merged file of the original was not valid SBML
  (kinetic laws referring to replaced identifiers such as `ATP_0`; 7 instead of 9
  unit definitions); the file of the port is valid.
- Element matching (`match`): identical matches and scores when the integer division
  is restored. With the true division one weak match (vinGlc/vGLT, 0.53) falls below
  0.5 and F6P/F6P changes from 0.71 to 0.64.
- (These comparisons need the database of 2011, which is not distributed; with a
  database built today results differ slightly, e.g. model similarities by about 0.005.)
- Model similarity: identical to six digits without the information content update.
  With it (what `semanticsbml similarity` does) values differ by less than 0.001.
  The original could do this update only when loading a folder of BioModels files,
  and there it stops with "maximum recursion depth exceeded" on the 2011 database
  (circles in the hierarchy), so this path cannot be compared.
