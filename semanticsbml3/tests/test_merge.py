import contextlib
import io

import libsbml
import pytest

from semanticSBML import cli


def merge(*paths):
    documents = [libsbml.readSBML(p) for p in paths]
    with contextlib.redirect_stdout(io.StringIO()):
        merged, report = cli.merge_documents(documents)
    text = libsbml.writeSBMLToString(merged)
    return libsbml.readSBMLFromString(text), report, text


def matched(report):
    """set of frozensets with the ids (including the model suffix) of joined elements"""
    return set(frozenset(i for mid, i, name in m['elements']) for m in report['matches'])


@pytest.mark.parametrize('first, second', [('hynne', 'teusink'), ('hynne_2011', 'teusink_2011')])
def test_merge_two_glycolysis_models(first, second, request):
    document, report, text = merge(request.getfixturevalue(first), request.getfixturevalue(second))
    model = document.getModel()
    assert cli.validation_errors(document) == []
    # 25 + 26 species, 15 of them are joined; 24 + 17 reactions, 10 are joined
    assert model.getNumSpecies() == 36
    assert model.getNumReactions() == 31
    assert model.getNumCompartments() == 2
    pairs = matched(report)
    for a, b in [('ATP_0', 'ATP_1'), ('ADP_0', 'ADP_1'), ('AMP_0', 'AMP_1'), ('Pyr_0', 'PYR_1'),
                 ('Glc_0', 'GLCi_1'), ('vPFK_0', 'vPFK_1'), ('vHK_0', 'vGLK_1'), ('cytosol_0', 'cytosol_1')]:
        assert frozenset([a, b]) in pairs
    # glucose 6-phosphate carries different ChEBI identifiers in the two models: not joined
    assert not any('G6P_0' in pair for pair in pairs)


def test_merge_is_reproducible(hynne, teusink):
    assert merge(hynne, teusink)[2] == merge(hynne, teusink)[2]


def test_merge_model_with_itself(teusink):
    original_document = libsbml.readSBML(teusink)
    original = original_document.getModel()
    document, report, text = merge(teusink, teusink)
    model = document.getModel()
    assert cli.validation_errors(document) == []
    assert model.getNumSpecies() == original.getNumSpecies()
    assert model.getNumReactions() == original.getNumReactions()
    assert model.getNumRules() == original.getNumRules()
    # the kinetic law of phosphofructokinase uses three function definitions
    assert model.getNumFunctionDefinitions() == original.getNumFunctionDefinitions()


def test_merge_keeps_kinetics_of_first_model(hynne, teusink):
    first_document = libsbml.readSBML(hynne)
    first = first_document.getModel()
    document, report, text = merge(hynne, teusink)
    merged_law = document.getModel().getReaction('vPFK_0').getKineticLaw()
    original_law = first.getReaction('vPFK').getKineticLaw()
    assert merged_law.getNumParameters() == original_law.getNumParameters()
    assert document.getModel().getSpecies('ATP_0').getInitialConcentration() == \
        first.getSpecies('ATP').getInitialConcentration()


def test_contradicting_models_are_reported(hynne, teusink):
    # in Teusink's model ATP, ADP and AMP are set by rules, in Hynne's by reactions;
    # with Teusink first the rules win and collide with Hynne's reactions
    document, report, text = merge(teusink, hynne)
    assert any('rule' in message for message in cli.validation_errors(document))


def test_merge_three_models(hynne, teusink):
    document, report, text = merge(hynne, teusink, teusink)
    assert document.getModel().getNumSpecies() == 36
    assert any(len(m['elements']) == 3 for m in report['matches'])


def test_cli_merge(hynne, teusink, tmp_path, capsys):
    out = tmp_path / 'merged.xml'
    assert cli.main(['merge', hynne, teusink, '-o', str(out)]) == 0
    report = capsys.readouterr().out
    assert 'Matched elements (27)' in report
    assert 'The merged model is valid SBML.' in report
    written = libsbml.readSBML(str(out))
    assert written.getModel().getNumSpecies() == 36
