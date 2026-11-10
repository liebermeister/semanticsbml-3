"""Browser interface: the functions behind the page, and the web server itself."""
import json
import threading
import urllib.error
import urllib.request
from http.server import HTTPServer

import pytest

from semanticSBML import webgui


def load(workspace, path):
    with open(path, encoding='utf-8') as f:
        return workspace.add_model(path, f.read())


@pytest.fixture
def workspace(hynne, teusink):
    w = webgui.Workspace()
    load(w, hynne)
    load(w, teusink)
    return w


def test_state_and_elements(workspace):
    models = workspace.state()['models']
    assert [m['name'] for m in models] == ['BIOMD0000000061', 'BIOMD0000000064']
    assert models[0]['species'] == 25
    rows = workspace.elements('m1')
    glucose = [r for r in rows if r['id'] == 'Glc'][0]
    assert ('ChEBI', 'CHEBI:17234') in [(a['db'], a['id']) for a in glucose['annotations']]


def test_add_and_remove_annotation(workspace):
    rows = workspace.elements('m1')
    target = [r for r in rows if not r['annotations'] and r['type'] == 'species'][0]
    workspace.add_annotation('m1', target['n'], 'ChEBI', 'CHEBI:15343', 'isVersionOf')
    after = workspace.elements('m1')[target['n']]['annotations']
    assert [(a['qualifier'], a['db'], a['id']) for a in after] == [('isVersionOf', 'ChEBI', 'CHEBI:15343')]
    with pytest.raises(webgui.GuiError):
        workspace.add_annotation('m1', target['n'], 'ChEBI', 'CHEBI:15343', 'is')       # already there
    with pytest.raises(webgui.GuiError):
        workspace.add_annotation('m1', target['n'], 'ChEBI', 'nonsense', 'is')
    workspace.remove_annotation('m1', target['n'], 'ChEBI', 'CHEBI:15343')
    assert workspace.elements('m1')[target['n']]['annotations'] == []
    assert 'CHEBI:15343' in workspace.sbml('m1')        # other elements still carry it


def test_search_prefers_chebi(workspace):
    rows = workspace.search('glucose', exact=True, element_type='species')
    assert (rows[0]['db'], rows[0]['id']) == ('ChEBI', 'CHEBI:17234')
    assert workspace.search('', exact=True) == []


def test_merge_match_similarity(workspace):
    report = workspace.merge(['m1', 'm2'])
    assert report['errors'] == [] and len(report['matches']) == 27
    assert report['counts']['species'] == 36
    assert len(workspace.state()['models']) == 3
    with pytest.raises(webgui.GuiError):
        workspace.merge(['m1'])
    matches = workspace.match(['m1', 'm2'], 0.6)['matches']
    assert any(set(e['id'] for e in m['elements']) == {'ATP'} for m in matches)
    similarity = workspace.similarity(['m1', 'm2'])
    assert similarity['matrix'][0][0] == 1.0 and 0.5 < similarity['matrix'][0][1] < 0.7


def test_bad_file_is_refused(workspace):
    with pytest.raises(webgui.GuiError):
        workspace.add_model('x.xml', '<notsbml/>')


def test_server(workspace):
    server = HTTPServer(('127.0.0.1', 0), None)
    port = server.server_address[1]
    server.RequestHandlerClass = webgui.make_handler(workspace, {'127.0.0.1:%d' % port})
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = 'http://127.0.0.1:%d' % port

    def get(path, headers=None, data=None):
        request = urllib.request.Request(base + path, data=data, headers=headers or {})
        try:
            with urllib.request.urlopen(request, timeout=30) as r:
                return r.status, r.read().decode('utf-8')
        except urllib.error.HTTPError as e:
            return e.code, e.read().decode('utf-8')

    try:
        status, page = get('/')
        assert status == 200 and '<title>semanticSBML</title>' in page
        status, body = get('/api/state')
        assert status == 200 and len(json.loads(body)['result']['models']) == 2
        assert get('/api/state', {'Host': 'evil.example'})[0] == 403                    # other web sites are refused
        assert get('/api/state', {'Origin': 'http://evil.example'})[0] == 403
        assert get('/api/merge', {'Content-Type': 'text/plain'}, b'{}')[0] == 415        # forms of other pages too
        status, body = get('/api/models/m1/download')
        assert status == 200 and '<sbml' in body
        assert get('/api/models/nothing/elements')[0] == 400
        assert get('/somewhere/else')[0] == 404
    finally:
        server.shutdown()
        server.server_close()
