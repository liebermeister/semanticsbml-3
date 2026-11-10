"""The annotation database built from ChEBI, Gene Ontology, SBO (semanticsbml build-database)."""
import os
import subprocess
import sys

import pytest

from conftest import DATA, DATABASE_2011, ROOT, model_path

built = pytest.mark.skipif(DATABASE_2011, reason='tests for a database built with build-database')


def run(*args):
    result = subprocess.run([sys.executable, '-m', 'semanticSBML.cli'] + list(args), cwd=ROOT,
                            env=dict(os.environ, SEMANTICSBML_DATA=DATA), capture_output=True, text=True)
    assert result.returncode in (0, 1), result.stderr
    return result.stdout


@built
def test_database_command_names_sources():
    out = run('database')
    assert 'ChEBI release' in out and 'Gene Ontology' in out


@built
def test_entry_added_after_2011_is_known():
    assert 'CHEBI:456216' in run('search', 'ADP(3-)', '--exact')        # ADP(3-), not in the 2011 database


def test_kegg_identifiers_understood_through_chebi():
    out = run('search', 'glucose', '--exact')
    assert 'CHEBI:17234' in out and 'C00293' in out                      # ChEBI names the KEGG compound


def test_merge_and_similarity_with_current_database(tmp_path):
    hynne = model_path(os.path.join('current', 'BIOMD0000000061.xml'))
    teusink = model_path(os.path.join('current', 'BIOMD0000000064.xml'))
    merged = str(tmp_path / 'merged.xml')
    out = run('merge', hynne, teusink, '-o', merged)
    assert 'The merged model is valid SBML.' in out
    assert '36 species, 31 reactions' in out
    out = run('similarity', hynne, teusink, '--tsv')
    value = float(out.strip().splitlines()[1].split('\t')[2])
    assert 0.5 < value < 0.7


@built
def test_no_kegg_data():
    """without --with-kegg nothing comes from KEGG: KEGG identifiers only as links of ChEBI and GO entries"""
    from libSBAnnotation import update
    with open(update.info_path(), encoding='utf-8') as f:
        assert 'KEGG Compound: ' not in f.read()
    folder = os.path.join(DATA, 'db')
    with open(os.path.join(folder, 'internal_db.version_1.8.external_item_kegg_formula')) as f:
        assert f.read().strip() == ''
    with open(os.path.join(folder, 'internal_db.version_1.8.takes_part')) as f:
        assert f.read().strip() == ''                                     # no reaction equations
    with open(os.path.join(folder, 'internal_db.version_1.8.string_repr_from_database')) as f:
        sources = set(line.split('\t')[2] for line in f if line.strip())
    assert sources <= {'chebi', 'go', 'sbo', 'taxonomy'}


@pytest.mark.nodb
def test_parse_obo(tmp_path):
    from libSBAnnotation import update
    path = tmp_path / 'x.obo'
    path.write_text('format-version: 1.2\ndata-version: 7\n\n[Term]\nid: CHEBI:1\nname: thing\n'
                    'synonym: "a \\"quoted\\" name" EXACT []\nxref: kegg.compound:C00001 {source="kegg"}\n'
                    'is_a: CHEBI:2 ! parent\nrelationship: RO:0000087 CHEBI:3 ! has role x\n'
                    'property_value: chemrof:generalized_empirical_formula "H2O" xsd:string\n\n'
                    '[Term]\nid: CHEBI:9\nname: old\nis_obsolete: true\nreplaced_by: CHEBI:1\n\n'
                    '[Typedef]\nid: RO:0000087\nname: has role\n')
    header, terms = update.parse_obo(str(path))
    assert header['data-version'] == '7' and len(terms) == 2
    t = terms[0]
    assert t['synonyms'] == ['a "quoted" name'] and t['xrefs'] == [('kegg.compound', 'C00001')]
    assert t['is_a'] == ['CHEBI:2'] and t['relationships'] == [('RO:0000087', 'CHEBI:3')] and t['formula'] == 'H2O'
    assert terms[1]['obsolete'] and terms[1]['replaced_by'] == ['CHEBI:1']
