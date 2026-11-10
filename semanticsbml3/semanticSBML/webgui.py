"""
Browser interface of semanticSBML (added in the Python 3 port).

    semanticsbml gui

starts a small web server on this computer only (127.0.0.1) and opens the page
in the browser.  Models are kept in memory while the program runs; nothing is
written to disk and nothing is sent to other computers.

The page is the file webgui.html next to this module; it talks to the server
through the JSON interface implemented by the class Workspace below.
"""

import contextlib
import io
import json
import os
import re
import sys
import threading
import urllib.parse
import warnings
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer

MAX_UPLOAD = 100 * 1024 * 1024
# databases shown first in search results
PREFERRED = ['ChEBI', 'KEGG Compound', 'KEGG Reaction', 'Gene Ontology', 'Systems Biology Ontology', 'Taxonomy',
             'Enzyme Nomenclature', 'UniProt']
ELEMENT_TYPES = ['species', 'reaction', 'compartment', 'parameter', 'rule', 'event', 'model']


class GuiError(Exception):
    """an error that is shown to the user as it is"""


def _search_type(element_type):
    t = (element_type or '').lower()
    if 'rule' in t:
        return 'rule'
    for known in ELEMENT_TYPES:
        if known in t:
            return known
    return None


class Workspace:
    """the models loaded in the browser session and everything that can be done with them"""

    def __init__(self):
        self.models = {}            # key -> {'name': str, 'document': libsbml.SBMLDocument}
        self.order = []
        self.counter = 0

    # -- models -------------------------------------------------------------
    def _get(self, key):
        if key not in self.models:
            raise GuiError('This model is no longer loaded.')
        return self.models[key]

    def _elements(self, document):
        from . import cli
        return cli.element_annotations(document)

    def add_document(self, name, document):
        self.counter += 1
        key = 'm%d' % self.counter
        names = set(m['name'] for m in self.models.values())
        base, n = name, 2
        while name in names:
            name = '%s (%d)' % (base, n)
            n += 1
        self.models[key] = {'name': name, 'document': document}
        self.order.append(key)
        return key

    def add_model(self, name, text):
        import libsbml
        document = libsbml.readSBMLFromString(text)
        fatal = [document.getError(i) for i in range(document.getNumErrors())
                 if document.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
        if fatal or document.getModel() is None:
            messages = '; '.join('line %s: %s' % (e.getLine(), e.getMessage().strip().split('\n')[-1]) for e in fatal[:3])
            raise GuiError('"%s" could not be read as SBML. %s' % (name, messages))
        name = re.sub(r'\.(xml|sbml)$', '', os.path.basename(name or 'model'), flags=re.I) or 'model'
        return self.add_document(name, document)

    def remove_model(self, key):
        self._get(key)
        del self.models[key]
        self.order.remove(key)

    def summary(self, key):
        entry = self._get(key)
        model = entry['document'].getModel()
        elements = self._elements(entry['document'])
        return {'key': key, 'name': entry['name'], 'title': model.getName() or model.getId() or '',
                'level': '%d.%d' % (entry['document'].getLevel(), entry['document'].getVersion()),
                'species': model.getNumSpecies(), 'reactions': model.getNumReactions(),
                'compartments': model.getNumCompartments(), 'parameters': model.getNumParameters(),
                'elements': len(elements), 'unannotated': sum(1 for e in elements if not e.getAnnotations())}

    def state(self):
        from libSBAnnotation import update
        database = 'annotation database of February 2011'
        if os.path.exists(update.info_path()):
            with open(update.info_path(), encoding='utf-8') as f:
                lines = [line.strip() for line in f if line.startswith(('Built', 'Updated'))]
            if lines:
                database = 'annotation database, %s' % lines[-1].lower()
        return {'models': [self.summary(key) for key in self.order], 'database': database}

    def sbml(self, key):
        import libsbml
        return libsbml.writeSBMLToString(self._get(key)['document'])

    # -- annotations --------------------------------------------------------
    def _annotation(self, a):
        uri = str(a)
        try:
            name = a.getName() or ''
        except Exception:
            name = ''
        return {'qualifier': a.qualifier, 'db': a.db, 'id': a.id, 'name': name,
                'url': uri.replace('http://', 'https://') if uri.startswith('http') else ''}

    def elements(self, key):
        rows = []
        for number, e in enumerate(self._elements(self._get(key)['document'])):
            rows.append({'n': number, 'type': e.type, 'id': e.id or '', 'name': e.name or '',
                         'annotations': [self._annotation(a) for a in e.getAnnotations()],
                         'problem': str(e.getReadErrors() or '')})
        return rows

    def _element(self, key, number):
        elements = self._elements(self._get(key)['document'])
        if not isinstance(number, int) or not 0 <= number < len(elements):
            raise GuiError('Element not found; please reload the model.')
        return elements[number]

    def qualifiers(self):
        from . import annotate
        return {'bio': list(annotate.bio_qualifiers), 'model': list(annotate.model_qualifiers)}

    def add_annotation(self, key, number, db, identifier, qualifier='is'):
        from . import annotate
        e = self._element(key, number)
        qualifier_type = 'bio'
        if qualifier not in annotate.bio_qualifiers:
            if qualifier not in annotate.model_qualifiers:
                raise GuiError('Unknown qualifier "%s".' % qualifier)
            qualifier_type = 'model'
        try:
            a = annotate.Annotation(db, identifier.strip(), qualifier_type, qualifier)
        except annotate.AnnotationError as error:
            raise GuiError(re.sub(r'<[^>]+>', ' ', str(error)).strip())
        same = [x for x in e.getAnnotations() if x.db == a.db and x.id == a.id]
        if same:
            raise GuiError('%s already has the annotation %s %s (as "%s"). Remove it first to change the relation.'
                           % (e.id or e.type, a.db, a.id, same[0].qualifier))
        e.addAnnotation(a)

    def remove_annotation(self, key, number, db, identifier, qualifier=None):
        e = self._element(key, number)
        hits = [a for a in e.getAnnotations()
                if a.db == db and a.id == identifier and (qualifier is None or a.qualifier == qualifier)]
        if not hits:
            raise GuiError('Annotation not found.')
        for a in hits:
            e.remAnnotation(a)

    def clear_annotations(self, key, number):
        self._element(key, number).unsetAnnotations()

    def auto_annotate(self, key):
        from . import annotate
        done = []
        failed = 0
        for e in self._elements(self._get(key)['document']):
            if e.isAnnotated() or e.type == 'model':
                continue
            try:
                e.addAnnotationAutomatic()
            except annotate.AnnotationError:
                failed += 1
                continue
            done.append({'type': e.type, 'id': e.id, 'name': e.name or '',
                         'annotations': [self._annotation(a) for a in e.getAnnotations()]})
        return {'annotated': done, 'left': failed}

    def copy_annotations(self, key, source_key, overwrite=False):
        from . import cli
        if key == source_key:
            raise GuiError('Choose a different model to copy from.')
        return {'copied': cli._copy_annotations(self._get(key)['document'], self._get(source_key)['document'], overwrite)}

    def search(self, query, exact=False, element_type=None, precision=0.8):
        import libSBAnnotation.annotate
        import libSBAnnotation.config
        import libSBAnnotation.database_new
        query = (query or '').strip().lower()
        if not query:
            return []
        db = libSBAnnotation.database_new.Meta_DB()
        groups = db.searchEntries(query) if exact else db.searchEntries(query, precision)
        allowed = None
        kind = _search_type(element_type)
        if kind:
            allowed = libSBAnnotation.config.Config().get('annotate', 'annotate_%s_with' % kind).split(',')
        rows = []
        for number, group in enumerate(groups[:200]):
            for database, identifier in sorted(group.items(), key=lambda x: (PREFERRED.index(x[0]) if x[0] in PREFERRED else len(PREFERRED), x)):
                if allowed is not None and database not in allowed:
                    continue
                try:
                    a = libSBAnnotation.annotate.Annotation(database, identifier, qualifier='is')
                except libSBAnnotation.annotate.AnnotationError:
                    continue
                row = self._annotation(a)
                row['group'] = number + 1
                rows.append(row)
        return rows

    # -- check, picture -----------------------------------------------------
    def check(self, key):
        import tempfile
        from .Model import Model
        from .sbmlcheck import Checker
        with tempfile.TemporaryDirectory() as folder:           # the checker reads a file
            path = os.path.join(folder, 'model.xml')
            with open(path, 'w', encoding='utf-8') as f:
                f.write(self.sbml(key))
            results = []
            for r in Checker(Model(path)).run():
                results.append({'name': r.name, 'errors': [str(m) for m in r.error_messages],
                                'warnings': [str(m) for m in r.warning_messages]})
        return results

    def picture(self, key, ids=False, max_participations=1000):
        import libSBAnnotation.config
        from . import sbml2dot
        document = self._get(key)['document']
        if document.getModel().getNumSpecies() > 400:
            raise GuiError('The model has more than 400 species; the picture would not be readable.')
        config = libSBAnnotation.config.Config()
        config.set('sbml2dot', 'colors', '2')
        config.set('sbml2dot', 'view_names', '1' if ids else '0')
        config.set('sbml2dot', 'max_number', '500')
        try:
            return sbml2dot.sbml2image(document, 'svg', max_participations)
        except FileNotFoundError:
            raise GuiError('Graphviz is not installed (needed for pictures): sudo apt install graphviz')

    # -- merge, match, similarity --------------------------------------------
    def merge(self, keys, name=''):
        import libsbml
        from . import cli
        if len(keys) < 2 or len(set(keys)) != len(keys):
            raise GuiError('Choose at least two different models to merge.')
        entries = [self._get(key) for key in keys]
        documents = [libsbml.readSBMLFromString(libsbml.writeSBMLToString(e['document'])) for e in entries]
        with contextlib.redirect_stdout(io.StringIO()):
            title = name or ' - '.join(d.getModel().getName() or d.getModel().getId() or 'model' for d in documents)
            merged, report = cli.merge_documents(documents, None, title)
            if not merged.getModel().isSetId():
                merged.getModel().setId('merged_model')

        def origin(new_id, mid):
            suffix = '_%s' % mid
            return new_id[:-len(suffix)] if new_id.endswith(suffix) else new_id

        matches = []
        for match in sorted(report['matches'], key=lambda x: (x['type'], x['merged_id'])):
            matches.append({'type': match['type'], 'merged_id': match['merged_id'], 'conflicts': match['conflicts'],
                            'elements': [{'model': mid + 1, 'id': origin(i, mid), 'name': n or ''}
                                         for mid, i, n in match['elements']]})
        unmatched = [{'model': mid + 1, 'type': t, 'id': origin(i, mid), 'name': n or ''}
                     for mid, t, i, n in report['unmatched']]
        text = libsbml.writeSBMLToString(merged)
        merged = libsbml.readSBMLFromString(text)
        errors = sorted(set(cli.validation_errors(libsbml.readSBMLFromString(text))))
        key = self.add_document('merged ' + ' + '.join(e['name'] for e in entries), merged)
        model = merged.getModel()
        return {'key': key, 'models': [e['name'] for e in entries], 'matches': matches, 'unmatched': unmatched,
                'warnings': [str(w) for w in report['warnings']], 'errors': errors,
                'counts': {'compartments': model.getNumCompartments(), 'species': model.getNumSpecies(),
                           'reactions': model.getNumReactions(), 'parameters': model.getNumParameters(),
                           'rules': model.getNumRules()}}

    def match(self, keys, minimum=0.5, propagation=True):
        from . import cli
        if len(keys) != 2 or keys[0] == keys[1]:
            raise GuiError('Choose two different models to compare.')
        entries = [self._get(key) for key in keys]
        rows = []
        for score, elements in cli.match_documents([e['document'] for e in entries], float(minimum), bool(propagation)):
            rows.append({'score': round(float(score), 4), 'type': elements[0][1],
                         'elements': [{'model': mid, 'id': i, 'name': n or ''} for mid, t, i, n in elements]})
        return {'models': [e['name'] for e in entries], 'matches': rows}

    def similarity(self, keys):
        import tempfile
        import libModelSimilarity.api as api
        if len(keys) < 2:
            raise GuiError('Choose at least two models.')
        entries = [self._get(key) for key in keys]
        with tempfile.TemporaryDirectory() as folder, warnings.catch_warnings():
            warnings.simplefilter('ignore')
            paths = []
            for number, key in enumerate(keys):
                paths.append(os.path.join(folder, 'model%d.xml' % number))
                with open(paths[-1], 'w', encoding='utf-8') as f:
                    f.write(self.sbml(key))
            collection = api.load_models(paths, [e['name'] for e in entries])
            models = collection.model_annotations
            names = [m.model_name for m in models]
            matrix = [[round(float(a.similarity(b)), 4) for b in models] for a in models]
            order = list(range(len(names)))
            if len(models) > 2:
                try:
                    clustered = api.cluster_order(collection)
                    order = sorted(order, key=lambda i: clustered.index(names[i]))
                except Exception:
                    pass
            api.reset()
        return {'names': [names[i] for i in order], 'matrix': [[matrix[i][j] for j in order] for i in order]}


# ---------------------------------------------------------------------------
# web server
# ---------------------------------------------------------------------------

def make_handler(workspace, allowed_hosts):
    page_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'webgui.html')
    lock = threading.Lock()

    class Handler(BaseHTTPRequestHandler):
        server_version = 'semanticSBML'

        def log_message(self, format, *args):
            if os.environ.get('SEMANTICSBML_DEBUG'):
                BaseHTTPRequestHandler.log_message(self, format, *args)

        # -- answers --
        def _send(self, status, body, content_type, extra=None):
            if isinstance(body, str):
                body = body.encode('utf-8')
            self.send_response(status)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.send_header('Content-Security-Policy',
                             "default-src 'self'; style-src 'self' 'unsafe-inline'; script-src 'self' 'unsafe-inline'; "
                             "img-src 'self' data: blob:")
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.end_headers()
            self.wfile.write(body)

        def _json(self, data, status=200):
            self._send(status, json.dumps(data), 'application/json; charset=utf-8')

        def _error(self, status, message):
            self._json({'error': message}, status)

        # -- protection: only pages of this server may use it --
        def _allowed(self):
            host = (self.headers.get('Host') or '').lower()
            if host not in allowed_hosts:
                return False
            origin = self.headers.get('Origin')
            if origin and urllib.parse.urlparse(origin).netloc.lower() not in allowed_hosts:
                return False
            return True

        def _run(self, function, *args, **kwargs):
            try:
                with lock, warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                    warnings.simplefilter('ignore')
                    return function(*args, **kwargs), None
            except GuiError as e:
                return None, (400, str(e))
            except Exception as e:
                if os.environ.get('SEMANTICSBML_DEBUG'):
                    import traceback
                    traceback.print_exc(file=sys.stderr)
                text = re.sub(r'<[^>]+>', ' ', str(e)).strip() or e.__class__.__name__
                return None, (500, text)

        def _answer(self, function, *args, **kwargs):
            result, error = self._run(function, *args, **kwargs)
            if error:
                self._error(*error)
            else:
                self._json({'result': result})

        def do_GET(self):
            if not self._allowed():
                return self._error(403, 'Request refused.')
            url = urllib.parse.urlparse(self.path)
            query = dict(urllib.parse.parse_qsl(url.query))
            parts = [urllib.parse.unquote(p) for p in url.path.strip('/').split('/')]
            if url.path in ('/', '/index.html'):
                with open(page_path, 'rb') as f:
                    return self._send(200, f.read(), 'text/html; charset=utf-8')
            if parts[:1] != ['api']:
                return self._error(404, 'Not found.')
            w = workspace
            if parts == ['api', 'state']:
                return self._answer(w.state)
            if parts == ['api', 'qualifiers']:
                return self._answer(w.qualifiers)
            if parts == ['api', 'search']:
                return self._answer(w.search, query.get('q', ''), query.get('exact') == '1', query.get('type'))
            if len(parts) == 4 and parts[1] == 'models':
                key, what = parts[2], parts[3]
                if what == 'elements':
                    return self._answer(w.elements, key)
                if what == 'check':
                    return self._answer(w.check, key)
                if what == 'picture':
                    try:
                        limit = int(query.get('max', '1000'))
                    except ValueError:
                        limit = 1000
                    return self._answer(w.picture, key, query.get('ids') == '1', limit)
                if what == 'download':
                    text, error = self._run(w.sbml, key)
                    if error:
                        return self._error(*error)
                    name = re.sub(r'[^A-Za-z0-9._-]+', '_', w.models[key]['name']).strip('_') or 'model'
                    return self._send(200, text, 'application/xml; charset=utf-8',
                                      {'Content-Disposition': 'attachment; filename="%s.xml"' % name})
            return self._error(404, 'Not found.')

        def do_POST(self):
            if not self._allowed():
                return self._error(403, 'Request refused.')
            if (self.headers.get('Content-Type') or '').split(';')[0].strip().lower() != 'application/json':
                return self._error(415, 'JSON expected.')
            try:
                length = int(self.headers.get('Content-Length') or 0)
            except ValueError:
                length = -1
            if length < 0 or length > MAX_UPLOAD:
                return self._error(413, 'The file is too large.')
            try:
                data = json.loads(self.rfile.read(length).decode('utf-8') or '{}')
            except (ValueError, UnicodeDecodeError):
                return self._error(400, 'Could not read the request.')
            if not isinstance(data, dict):
                return self._error(400, 'Could not read the request.')
            parts = [urllib.parse.unquote(p) for p in urllib.parse.urlparse(self.path).path.strip('/').split('/')]
            w = workspace
            g = data.get
            if parts == ['api', 'models']:
                return self._answer(w.add_model, str(g('name', 'model')), str(g('content', '')))
            if parts == ['api', 'merge']:
                return self._answer(w.merge, list(g('models', [])), str(g('name', '')))
            if parts == ['api', 'match']:
                return self._answer(w.match, list(g('models', [])), g('minimum', 0.5), g('propagation', True))
            if parts == ['api', 'similarity']:
                return self._answer(w.similarity, list(g('models', [])))
            if len(parts) == 4 and parts[1] == 'models':
                key, what = parts[2], parts[3]
                if what == 'remove':
                    return self._answer(w.remove_model, key)
                if what == 'add-annotation':
                    return self._answer(w.add_annotation, key, g('element'), str(g('db', '')), str(g('id', '')),
                                        str(g('qualifier', 'is')))
                if what == 'remove-annotation':
                    return self._answer(w.remove_annotation, key, g('element'), str(g('db', '')), str(g('id', '')),
                                        g('qualifier'))
                if what == 'clear-annotations':
                    return self._answer(w.clear_annotations, key, g('element'))
                if what == 'auto-annotate':
                    return self._answer(w.auto_annotate, key)
                if what == 'copy-annotations':
                    return self._answer(w.copy_annotations, key, str(g('source', '')), bool(g('overwrite', False)))
            return self._error(404, 'Not found.')

    return Handler


def serve(port=8765, open_browser=True, files=()):
    """start the browser interface; returns when the user presses Ctrl-C"""
    log = sys.stderr.write
    log('Loading the annotation database ...\n')
    import libSBAnnotation.database_new
    from . import annotate          # noqa: F401  (loads the resource tables)
    libSBAnnotation.database_new.Meta_DB()
    workspace = Workspace()
    for path in files:
        with open(path, encoding='utf-8', errors='replace') as f:
            workspace.add_model(os.path.basename(path), f.read())
    server = None
    for candidate in [port] + list(range(port + 1, port + 20)):
        try:
            hosts = set(['127.0.0.1:%d' % candidate, 'localhost:%d' % candidate])
            server = HTTPServer(('127.0.0.1', candidate), make_handler(workspace, hosts))
            port = candidate
            break
        except OSError:
            continue
    if server is None:
        raise GuiError('No free port found near %d; choose another one with --port.' % port)
    url = 'http://127.0.0.1:%d/' % port
    log('semanticSBML is running at  %s\n(only reachable from this computer; stop it with Ctrl-C)\n' % url)
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        log('\nstopped\n')
    finally:
        server.server_close()
    return 0
