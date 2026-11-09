import shutil

import libsbml
import pytest

from semanticSBML import cli, sbml2dot
from semanticSBML.Model import Model
from semanticSBML.sbmlcheck import Checker


def test_check_reports(hynne_2011):
    results = dict((r.name, r) for r in Checker(Model(hynne_2011)).run())
    assert set(results) == {'Atom number balance', 'Duplicate elements', 'Annotations'}
    assert all(r.errors == 0 for r in results.values())
    # the glycerol transport reaction carries the volume ratio 59 as stoichiometry
    assert any('vdifGlyc' in m and 'atom number conservation' in m
               for m in results['Atom number balance'].warning_messages)
    # four reactions and the species P have no annotation
    assert results['Annotations'].warnings == 5


def test_check_cli_exit_code(hynne_2011, capsys):
    assert cli.main(['check', hynne_2011]) == 1          # warnings, no errors
    assert 'warnings have been issued' in capsys.readouterr().out


def test_broken_file_gives_readable_error(tmp_path, capsys):
    broken = tmp_path / 'broken.xml'
    broken.write_text('<sbml><model>')
    assert cli.main(['check', str(broken)]) == 2
    assert 'Could not read' in capsys.readouterr().err


def test_dot_graph(teusink):
    dot = sbml2dot.sbml2dot_string(libsbml.readSBML(teusink))
    assert dot.startswith('digraph')
    assert 'Hexokinase' in dot or 'vGLK' in dot


@pytest.mark.skipif(shutil.which('dot') is None, reason='graphviz is not installed')
def test_svg_picture(teusink, tmp_path):
    out = tmp_path / 'model.svg'
    assert cli.main(['view', teusink, '-o', str(out)]) == 0
    assert '<svg' in out.read_text()
