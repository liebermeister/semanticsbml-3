import libsbml
import pytest

import libModelSimilarity.api as api
from semanticSBML import cli

from conftest import model_path


def test_model_similarity_matrix(hynne, teusink, teusink_2011):
    other = model_path('BIOMD0000000239.xml')          # glucose-stimulated insulin secretion
    names, matrix = api.similarity_matrix([hynne, teusink, teusink_2011, other],
                                          ['hynne', 'teusink', 'teusink_2011', 'other'])
    sim = dict(((a, b), matrix[i][j]) for i, a in enumerate(names) for j, b in enumerate(names))
    for name in names:
        assert sim[(name, name)] == pytest.approx(1.0)
    for a in names:
        for b in names:
            assert 0. <= sim[(a, b)] <= 1.0 + 1e-9
            assert sim[(a, b)] == pytest.approx(sim[(b, a)])
    # the same model with urn:miriam and with identifiers.org annotations
    assert sim[('teusink', 'teusink_2011')] > 0.95
    # two models of yeast glycolysis are much closer to each other than to another model
    assert sim[('hynne', 'teusink')] > 0.5
    assert sim[('hynne', 'other')] < 0.3
    assert sim[('teusink', 'other')] < 0.3


def test_similarity_is_repeatable(hynne, teusink):
    # the similarity code keeps loaded models in global variables; a second call must not be influenced
    first = api.similarity_matrix([hynne, teusink])[1]
    second = api.similarity_matrix([hynne, teusink])[1]
    assert first == second


def test_match_elements(hynne, teusink):
    documents = [libsbml.readSBML(hynne), libsbml.readSBML(teusink)]
    matches = cli.match_documents(documents, 0.6)
    pairs = dict((frozenset((mid, i) for mid, t, i, n in elements), score) for score, elements in matches)
    for a, b in [('ATP', 'ATP'), ('Pyr', 'PYR'), ('Glc', 'GLCi'), ('vHK', 'vGLK'), ('vPK', 'vPYK')]:
        assert pairs[frozenset([(1, a), (2, b)])] > 0.9
    # every element is matched at most once
    used = [(mid, i) for score, elements in matches for mid, t, i, n in elements]
    assert len(used) == len(set(used))
    # only elements of the same kind are matched
    for score, elements in matches:
        assert len(set(t for mid, t, i, n in elements)) == 1


def test_cli_similarity_and_match(hynne, teusink, capsys):
    assert cli.main(['similarity', hynne, teusink, '--tsv']) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert len(lines) == 3 and float(lines[1].split('\t')[2]) > 0.5

    assert cli.main(['match', hynne, teusink, '--tsv']) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].split('\t') == ['similarity', 'type', 'model_1', 'model_2']
    assert any(l.split('\t')[1:] == ['species', 'ATP', 'ATP'] for l in lines)
