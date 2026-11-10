"""Compare the port with results of the original program (see tests/reference_py2/README.md)."""
import json
import os
import warnings

import libsbml
import pytest

from conftest import model_path

REF = os.path.join(os.path.dirname(__file__), 'reference_py2')
FILES = ['BIOMD0000000061.xml', 'BIOMD0000000064.xml', 'BIOMD0000000239.xml']


def ref(name):
    with open(os.path.join(REF, name + '.json')) as f:
        return json.load(f)


def test_annotations_as_original():
    from semanticSBML import annotate
    for name, rows in ref('annotations').items():
        doc = libsbml.readSBML(model_path(name))
        mine = set()
        for ea in annotate.ModelElementsAnnotations(doc.getModel(), True).getElementAnnotations():
            for a in ea.getAnnotations():
                mine.add((str(ea.id), str(a.qualifier), str(a.db), str(a.id)))
        assert mine == set(tuple(r[:4]) for r in rows)


def test_check_as_original():
    from semanticSBML.Model import Model
    from semanticSBML.sbmlcheck import Checker
    for name, rows in ref('check').items():
        mine = sorted([str(r.name), r.errors, r.warnings] for r in Checker(Model(model_path(name))).run())
        assert mine == [r[:3] for r in rows]


def test_merge_matches_as_original():
    """same groups as the original, except the documented correction (ATP/ADP no longer
    paired with pooled species) and a tie between two glucose species of model 61"""
    from semanticSBML import merge
    docs = [libsbml.readSBML(model_path(FILES[0])), libsbml.readSBML(model_path(FILES[1]))]
    merger = merge.Merger(docs, lambda *a: None, lambda *a: None)
    mine = set(tuple(sorted(str(e.id) for e in t.getElements())) for t in merger.getTuples())
    orig = set(tuple(sorted(e[1] for e in m)) for m in ref('merge')['matches'])
    assert orig - mine == {('ADP_0', 'P_1'), ('ATP_0', 'SUM_P_1'), ('GLCo_1', 'GlcX0_0')}
    assert mine - orig == {('ADP_0', 'ADP_1'), ('ATP_0', 'ATP_1'), ('GLCo_1', 'GlcX_0')}


def test_model_similarity_as_original():
    """without the information content update the numbers are those of the original"""
    from semanticSBML import annotate
    from libModelSimilarity import api, cluster
    api.reset()
    coll = cluster.DatabaseAnnotation()
    for name in FILES:
        doc = libsbml.readSBML(model_path(name))
        mea = annotate.ModelElementsAnnotations(doc.getModel(), True)
        coll.model_annotations.append(cluster.ModelAnnotation(mea, name))
    ms = coll.model_annotations
    expected = ref('similarity')
    for i in range(3):
        for j in range(i, 3):
            key = '%s~%s' % (ms[i].model_name, ms[j].model_name)
            assert float(ms[i].similarity(ms[j])) == pytest.approx(expected[key], abs=2e-6)
    api.reset()


def test_match_as_original_with_integer_division(monkeypatch):
    """element matching equals the original when the Python 2 integer division is restored"""
    from semanticSBML import similarity
    monkeypatch.setattr(similarity, 'HEIGHT_INTEGER_DIVISION', True)
    docs = [libsbml.readSBML(model_path(FILES[0])), libsbml.readSBML(model_path(FILES[1]))]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        matcher = similarity.ModelMatcher(docs, [1, 2], 1, True)
        matcher.pick_propagation_type('local')
        mine = {}
        for match in matcher.match(0.5):
            ids = [d['element_id'] for d in match]
            mids = [d['mid'] for d in match]
            key = tuple(sorted('%s:%s' % (d['mid'], d['element_id']) for d in match))
            mine[key] = float(matcher.global_esm.get_similarity(mids[0], mids[1], ids[0], ids[1]))
    orig = {tuple(r[1:]): r[0] for r in ref('match')['local']}
    assert set(mine) == set(orig)
    for key in orig:
        assert mine[key] == pytest.approx(orig[key], abs=2e-6)
