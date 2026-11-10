"""
Command line interface of semanticSBML.

    semanticsbml check        MODEL.xml ...
    semanticsbml view         MODEL.xml -o picture.svg
    semanticsbml annotations  MODEL.xml
    semanticsbml search       "glucose"
    semanticsbml annotate     MODEL.xml -o OUT.xml --auto
    semanticsbml merge        A.xml B.xml -o MERGED.xml
    semanticsbml match        A.xml B.xml
    semanticsbml similarity   A.xml B.xml C.xml

Run "semanticsbml COMMAND --help" for the options of each command.
"""

import argparse
import os
import sys


class CliError(Exception):
    pass


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def read_document(path):
    """read an SBML file, stop with a readable message if it is broken"""
    import libsbml
    if not os.path.isfile(path):
        raise CliError('File not found: %s' % path)
    document = libsbml.readSBML(path)
    fatal = [document.getError(i) for i in range(document.getNumErrors())
             if document.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
    if fatal or document.getModel() is None:
        messages = '\n'.join('  line %s: %s' % (e.getLine(), e.getMessage().strip()) for e in fatal[:10])
        raise CliError('Could not read %s as SBML\n%s' % (path, messages))
    return document


def write_document(document, path):
    import libsbml
    if not libsbml.writeSBMLToFile(document, path):
        raise CliError('Could not write %s' % path)


def make_box(string):
    line = '+' + '-' * (len(string) + 2) + '+'
    return '%s\n| %s |\n%s' % (line, string, line)


def element_annotations(document):
    """list of ElementAnnotations of all annotatable elements, in file order (model first)"""
    from . import annotate
    elements = annotate.ModelElementsAnnotations(document.getModel(), True).getElementAnnotations()
    return list(reversed(elements))


def find_element(elements, element_id):
    """find an element by its id (or 'model' for the model itself)"""
    hits = [e for e in elements if e.id == element_id]
    if not hits and element_id == 'model':
        hits = [e for e in elements if e.type == 'model']
    if not hits:
        hits = [e for e in elements if e.metaid == element_id]
    if not hits:
        raise CliError('No element with id "%s" in the model' % element_id)
    if len(hits) > 1:
        raise CliError('Id "%s" is ambiguous (%s)' % (element_id, ', '.join(e.type for e in hits)))
    return hits[0]


def print_table(rows, out=sys.stdout):
    if not rows:
        return
    widths = [max(len(str(row[i])) for row in rows) for i in range(len(rows[0]) - 1)]
    for row in rows:
        cells = [str(c).ljust(w) for c, w in zip(row, widths)] + [str(row[-1])]
        out.write('  '.join(cells).rstrip() + '\n')


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------

def cmd_check(args):
    from .Model import Model
    from .sbmlcheck import Checker

    total_errors = total_warnings = 0
    for f in args.files:
        read_document(f)        # readable error message for broken files
        file_errors = file_warnings = 0
        if not args.quiet:
            print(make_box('CHECKING FILE %s' % f))
        for r in Checker(Model(f)).run():
            file_errors += r.errors
            file_warnings += r.warnings
            if not args.quiet and (r.errors or r.warnings or args.verbose):
                print('%s returned %d errors and %d warnings.' % (r.name, r.errors, r.warnings))
                for m in r.error_messages:
                    print(' - ERROR: ' + m)
                for m in r.warning_messages:
                    print(' - WARNING: ' + m)
                print()
        if not args.quiet:
            print('%d errors have been found in this file, %d warnings have been issued.\n' % (file_errors, file_warnings))
        total_errors += file_errors
        total_warnings += file_warnings
    if not args.quiet and len(args.files) > 1:
        print('A total of %d errors and %d warnings have been found in all files.' % (total_errors, total_warnings))
    return (1 if total_warnings else 0) + (2 if total_errors else 0)


# ---------------------------------------------------------------------------
# view
# ---------------------------------------------------------------------------

def cmd_view(args):
    import libSBAnnotation.config
    from . import sbml2dot

    document = read_document(args.file)
    config = libSBAnnotation.config.Config()
    config.set('sbml2dot', 'colors', '0' if args.black_and_white else '2')
    config.set('sbml2dot', 'view_names', '1' if args.ids else '0')
    config.set('sbml2dot', 'max_number', str(args.max_species))

    out = args.outfile
    fmt = (os.path.splitext(out)[1][1:].lower() if out else '') or 'dot'
    if fmt == 'dot':
        data = sbml2dot.sbml2dot_string(document, args.max_participations)
    else:
        data = sbml2dot.sbml2image(document, fmt, args.max_participations)
    if not out:
        sys.stdout.write(data)
    elif isinstance(data, bytes):
        with open(out, 'wb') as f:
            f.write(data)
    else:
        with open(out, 'w', encoding='utf-8') as f:
            f.write(data)
    return 0


# ---------------------------------------------------------------------------
# annotations (list)
# ---------------------------------------------------------------------------

def cmd_annotations(args):
    document = read_document(args.file)
    elements = element_annotations(document)
    rows = []
    missing = 0
    for e in elements:
        annotations = e.getAnnotations()
        if not annotations:
            missing += 1
            if not args.annotated:
                rows.append((e.type, e.id, e.name, '-', '-', '-', ''))
            continue
        if args.missing:
            continue
        for a in annotations:
            rows.append((e.type, e.id, e.name, a.qualifier, a.db, a.id, a.getName()))
    if args.tsv:
        print('\t'.join(('type', 'id', 'name', 'qualifier', 'database', 'identifier', 'identifier_name')))
        for row in rows:
            print('\t'.join(str(c) for c in row))
    else:
        print_table([('TYPE', 'ID', 'NAME', 'QUALIFIER', 'DATABASE', 'IDENTIFIER', 'MEANING')] + rows)
        print('\n%d of %d elements have no annotation.' % (missing, len(elements)))
    for e in elements:
        if e.getReadErrors():
            sys.stderr.write('Problem in %s %s: %s\n' % (e.type, e.id, e.getReadErrors()))
    return 0


# ---------------------------------------------------------------------------
# search the annotation database
# ---------------------------------------------------------------------------

def cmd_search(args):
    import libSBAnnotation.database_new
    import libSBAnnotation.annotate
    import libSBAnnotation.config

    db = libSBAnnotation.database_new.Meta_DB()
    query = args.query.lower()
    if args.exact:
        groups = db.searchEntries(query)
    else:
        groups = db.searchEntries(query, args.precision)
    allowed = None
    if args.type:
        allowed = libSBAnnotation.config.Config().get('annotate', 'annotate_%s_with' % args.type).split(',')
    rows = []
    for number, group in enumerate(groups):
        for database, identifier in sorted(group.items()):
            if allowed is not None and database not in allowed:
                continue
            try:
                a = libSBAnnotation.annotate.Annotation(database, identifier, qualifier='is')
                rows.append((number + 1, a.getName(), database, identifier, str(a)))
            except libSBAnnotation.annotate.AnnotationError:
                pass
    if not rows:
        print('Nothing found for "%s".' % args.query)
        return 1
    print_table([('#', 'NAME', 'DATABASE', 'IDENTIFIER', 'URI')] + rows)
    return 0


# ---------------------------------------------------------------------------
# annotate (modify)
# ---------------------------------------------------------------------------

def _copy_annotations(document, source_document, overwrite):
    """copy annotations from elements of another model that have the same id"""
    from . import annotate
    source = dict(((e.type, e.id), e) for e in element_annotations(source_document))
    copied = 0
    for e in element_annotations(document):
        other = source.get((e.type, e.id))
        if e.type == 'model' and other is None:
            others = [x for (t, i), x in source.items() if t == 'model']
            other = others[0] if others else None
        if other is None or not other.getAnnotations():
            continue
        if e.isAnnotated():
            if not overwrite:
                continue
            e.unsetAnnotations()
        for a in other.getAnnotations():
            e.addAnnotation(annotate.Annotation(a.db, a.id, a.qualifier_type, a.qualifier))
        copied += 1
    return copied


def cmd_annotate(args):
    from . import annotate

    document = read_document(args.file)
    elements = element_annotations(document)
    log = sys.stderr.write

    for element_id in args.clear or []:
        find_element(elements, element_id).unsetAnnotations()
        log('cleared annotations of %s\n' % element_id)
    if args.clear_all:
        for e in elements:
            e.unsetAnnotations()
        log('cleared all annotations\n')

    for element_id, database, identifier in args.remove or []:
        e = find_element(elements, element_id)
        wanted = annotate.Annotation(database, identifier, qualifier='is')
        hits = [a for a in e.getAnnotations() if a.db == wanted.db and a.id == wanted.id]
        if not hits:
            raise CliError('%s has no annotation %s %s' % (element_id, database, identifier))
        for a in hits:
            e.remAnnotation(a)
        log('removed %s %s from %s\n' % (wanted.db, wanted.id, element_id))

    if args.copy_from:
        n = _copy_annotations(document, read_document(args.copy_from), args.overwrite)
        elements = element_annotations(document)
        log('copied annotations to %d elements from %s\n' % (n, args.copy_from))

    for item in args.add or []:
        if len(item) not in (3, 4):
            raise CliError('--add needs ELEMENT DATABASE IDENTIFIER [QUALIFIER]')
        element_id, database, identifier = item[:3]
        qualifier = item[3] if len(item) == 4 else 'is'
        e = find_element(elements, element_id)
        qualifier_type = 'bio'
        if qualifier not in annotate.bio_qualifiers:
            if qualifier in annotate.model_qualifiers:
                qualifier_type = 'model'
            else:
                raise CliError('Unknown qualifier "%s". Known: %s' % (qualifier, ', '.join(annotate.bio_qualifiers)))
        a = annotate.Annotation(database, identifier, qualifier_type, qualifier)
        if any(x.db == a.db and x.id == a.id and x.qualifier == a.qualifier for x in e.getAnnotations()):
            log('%s already has %s %s\n' % (element_id, a.db, a.id))
            continue
        e.addAnnotation(a)
        log('added %s %s (%s) to %s  [%s]\n' % (a.db, a.id, qualifier, element_id, a.getName()))

    if args.auto:
        done = failed = 0
        for e in elements:
            if e.isAnnotated() or e.type == 'model':
                continue
            try:
                e.addAnnotationAutomatic()
            except annotate.AnnotationError:
                failed += 1
                if args.verbose:
                    log('  no unique database entry for %s "%s"\n' % (e.type, e.name or e.id))
                continue
            done += 1
            log('  %s %s "%s" -> %s\n' % (e.type, e.id, e.name,
                                           ', '.join('%s %s' % (a.db, a.id) for a in e.getAnnotations())))
        log('automatic annotation: %d elements annotated, %d left without annotation\n' % (done, failed))

    if args.outfile:
        write_document(document, args.outfile)
        log('written %s\n' % args.outfile)
    else:
        import libsbml
        sys.stdout.write(libsbml.writeSBMLToString(document))
    return 0


# ---------------------------------------------------------------------------
# merge
# ---------------------------------------------------------------------------

def validation_errors(document):
    """error messages of the libSBML consistency check"""
    import libsbml
    document.checkConsistency()
    messages = []
    for i in range(document.getNumErrors()):
        e = document.getError(i)
        if e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            messages.append(e.getMessage().strip().split('\n')[-1].strip())
    return messages


def merge_documents(documents, status=None, name=''):
    """
    merge libSBML documents; if elements are in conflict the earlier model wins.
    @return: (merged libsbml document, report dictionary)
    """
    from . import merge

    warnings = []
    merger = merge.Merger(documents, status or (lambda *a: None),
                          lambda title, text='': warnings.append(text or title))
    matches = []
    for t in merger.getTuples():
        me = t.getMergedElement()
        conflicts = [label for label, flag in (('annotation', me.bioEntity_inConflict),
                                               ('quantity/unit/participants', me.bioQuantity_inConflict),
                                               ('location', me.bioQLocation_inConflict),
                                               ('values/formula', me.modelStatement_inConflict)) if flag]
        matches.append({'type': t.type,
                        'elements': [(e.merge_model.mid, e.id, e.name) for e in t.getElements()],
                        'conflicts': conflicts})
    merger.resolveUsingModel(merger.merge_models)
    for match, t in zip(matches, merger.getTuples()):
        match['merged_id'] = t.getMergedElement().getId()
    document = merger.createSBMLmodel(name)
    unmatched = [(e.merge_model.mid, e.type, e.id, e.name) for e in merger.getNonTupleElements()]
    report = {'matches': matches, 'unmatched': unmatched,
              'warnings': warnings + list(getattr(merger, 'warnings', []))}
    return document, report


def cmd_merge(args):
    import io
    import contextlib
    import libsbml

    if len(args.files) < 2:
        raise CliError('merge needs at least two models')
    documents = [read_document(f) for f in args.files]
    log = sys.stderr.write
    status = (lambda message, *a: log(message + '\n')) if args.verbose else None
    with contextlib.redirect_stdout(io.StringIO()):       # the merge code prints debugging output
        name = args.name or ' - '.join(d.getModel().getName() or d.getModel().getId() or 'model' for d in documents)
        merged, report = merge_documents(documents, status, name)
        if not merged.getModel().isSetId():
            merged.getModel().setId('merged_model')

    def origin_id(new_id, mid):                            # ids got the suffix _<model number> for merging
        suffix = '_%s' % mid
        return new_id[:-len(suffix)] if new_id.endswith(suffix) else new_id

    out = sys.stderr if not args.outfile else sys.stdout
    out.write('Models:\n')
    for i, f in enumerate(args.files):
        out.write('  [%d] %s\n' % (i + 1, f))
    out.write('\nMatched elements (%d):\n' % len(report['matches']))
    rows = []
    for match in sorted(report['matches'], key=lambda x: (x['type'], x['merged_id'])):
        elements = ' = '.join('[%d] %s%s' % (mid + 1, origin_id(i, mid), (' "%s"' % n) if n else '')
                              for mid, i, n in match['elements'])
        rows.append((match['type'], match['merged_id'], elements,
                     ('CONFLICT: ' + ', '.join(match['conflicts'])) if match['conflicts'] else ''))
    print_table(rows, out)
    conflicts = [m for m in report['matches'] if m['conflicts']]
    if conflicts:
        out.write('\n%d matched elements differ between the models; the values of the first model in the list were used.\n' % len(conflicts))
    if args.unmatched:
        out.write('\nElements without partner (%d), taken over unchanged:\n' % len(report['unmatched']))
        print_table([(t, '[%d] %s' % (mid + 1, origin_id(i, mid)), n) for mid, t, i, n in report['unmatched']], out)
    else:
        out.write('\n%d elements without partner were taken over unchanged (list them with --unmatched).\n' % len(report['unmatched']))
    for w in report['warnings']:
        out.write('\nWarning: %s\n' % w)

    text = libsbml.writeSBMLToString(merged)
    errors = validation_errors(libsbml.readSBMLFromString(text))
    model = merged.getModel()
    out.write('\nMerged model: %d compartments, %d species, %d reactions, %d parameters, %d rules.\n' % (
        model.getNumCompartments(), model.getNumSpecies(), model.getNumReactions(),
        model.getNumParameters(), model.getNumRules()))
    if errors:
        out.write('\nThe merged model is NOT valid SBML (%d problems). They usually arise where the models\n'
                  'describe the same element in contradicting ways and need to be repaired by hand:\n' % len(errors))
        for message in sorted(set(errors)):
            out.write('  - %s\n' % message)
    else:
        out.write('The merged model is valid SBML.\n')

    if args.outfile:
        with open(args.outfile, 'w', encoding='utf-8') as f:
            f.write(text)
        out.write('Written to %s\n' % args.outfile)
    else:
        sys.stdout.write(text)
    return 1 if errors else 0


# ---------------------------------------------------------------------------
# similarity of whole models
# ---------------------------------------------------------------------------

def _model_files(args_files, directory):
    files = list(args_files or [])
    if directory:
        if not os.path.isdir(directory):
            raise CliError('Not a directory: %s' % directory)
        files += sorted(os.path.join(directory, f) for f in os.listdir(directory)
                        if f.lower().endswith(('.xml', '.sbml')))
    for f in files:
        if not os.path.isfile(f):
            raise CliError('File not found: %s' % f)
    return files


def cmd_similarity(args):
    import warnings
    import libModelSimilarity.api as api

    files = _model_files(args.files, args.dir)
    if args.query:
        files = [args.query] + [f for f in files if os.path.abspath(f) != os.path.abspath(args.query)]
    if len(files) < 2:
        raise CliError('similarity needs at least two models')
    for f in files:
        read_document(f)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        collection = api.load_models(files)
        models = collection.model_annotations
        names = [m.model_name for m in models]

        if args.query:
            # rank all models by their similarity to the query model
            ranking = sorted(((models[0].similarity(m), m.model_name) for m in models[1:]), reverse=True)
            if args.tsv:
                print('similarity\tmodel')
                for sim, name in ranking:
                    print('%.4f\t%s' % (sim, name))
            else:
                print('Similarity to %s:' % names[0])
                print_table([('%.3f' % sim, name) for sim, name in ranking])
            return 0

        if args.cluster and len(models) > 2:
            order = api.cluster_order(collection)
            models = sorted(models, key=lambda m: order.index(m.model_name))
            names = [m.model_name for m in models]
        matrix = [[a.similarity(b) for b in models] for a in models]
    if args.tsv:
        print('\t'.join([''] + names))
        for name, row in zip(names, matrix):
            print('\t'.join([name] + ['%.4f' % x for x in row]))
    else:
        width = max(len(n) for n in names)
        for number, name in enumerate(names):
            print('[%d] %s' % (number + 1, name))
        print()
        print(' ' * 5 + ''.join(('[%d]' % (i + 1)).rjust(7) for i in range(len(names))))
        for number, row in enumerate(matrix):
            print(('[%d]' % (number + 1)).ljust(5) + ''.join(('%.3f' % x).rjust(7) for x in row))
    return 0


# ---------------------------------------------------------------------------
# matching elements of two or more models by similarity
# ---------------------------------------------------------------------------

def match_documents(documents, minimum_similarity=0.5, propagation=True, neighbours=1):
    """
    find corresponding elements in several models
    @return: list of matches (score, [(model number, element type, id, name), ...]), best first
    """
    import warnings
    from . import similarity

    model_ids = list(range(1, len(documents) + 1))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')            # elements without annotation give 0/0
        matcher = similarity.ModelMatcher(documents, model_ids, neighbours, True)
        if not propagation:
            matcher.pick_propagation_type('local')
        names = {}
        for mid, mea in zip(model_ids, matcher.meas):
            for ea in mea.getElementAnnotations():
                names[(mid, ea.id)] = ea.name
        result = []
        for match in matcher.match(minimum_similarity):
            score = matcher.get_similarity([d['mid'] for d in match], [d['element_id'] for d in match])
            pair_scores = []
            for i in range(len(match)):
                for j in range(i + 1, len(match)):
                    pair_scores.append(matcher.global_esm.get_similarity(
                        match[i]['mid'], match[j]['mid'], match[i]['element_id'], match[j]['element_id']))
            elements = sorted((d['mid'], d['element_type'], d['element_id'], names.get((d['mid'], d['element_id']), ''))
                              for d in match)
            result.append((min(pair_scores) if pair_scores else float(score), elements))
    result.sort(key=lambda x: (-x[0], x[1]))
    return result


def cmd_match(args):
    if len(args.files) < 2:
        raise CliError('match needs at least two models')
    documents = [read_document(f) for f in args.files]
    matches = match_documents(documents, args.min, not args.local, args.neighbours)
    if args.tsv:
        print('\t'.join(['similarity', 'type'] + ['model_%d' % (i + 1) for i in range(len(documents))]))
        for score, elements in matches:
            ids = dict((mid, i) for mid, t, i, n in elements)
            print('\t'.join(['%.4f' % score, elements[0][1]] + [ids.get(i + 1, '') for i in range(len(documents))]))
        return 0
    for i, f in enumerate(args.files):
        print('[%d] %s' % (i + 1, f))
    print('\nCorresponding elements (similarity at least %s): %d\n' % (args.min, len(matches)))
    rows = []
    for score, elements in matches:
        rows.append(('%.3f' % score, elements[0][1],
                     ' = '.join('[%d] %s%s' % (mid, i, (' "%s"' % n) if n else '') for mid, t, i, n in elements)))
    print_table(rows)
    return 0


# ---------------------------------------------------------------------------
# argument parser
# ---------------------------------------------------------------------------
# database
# ---------------------------------------------------------------------------

def cmd_update_database(args):
    from libSBAnnotation import update
    try:
        report = update.run(download_dir=args.download_dir, use_existing=args.use_existing,
                            kegg=not args.no_kegg, all_chebi=args.all_chebi, only=args.only)
    except update.UpdateError as e:
        raise CliError(str(e))
    print()
    for line in report:
        print(line)
    return 0


def cmd_database(args):
    from libSBAnnotation import config, update
    data = config.data_dir()
    print('data folder: %s' % data)
    info = os.path.join(data, update.INFO_FILE)
    if os.path.exists(info):
        with open(info, encoding='utf-8') as f:
            print(f.read().rstrip())
    else:
        print('annotation database of February 2011 (never updated)')
    return 0


# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        prog='semanticsbml',
        description='Check, annotate, merge and compare SBML models using their semantic (MIRIAM) annotations.')
    parser.add_argument('--data', metavar='DIR',
                        help='directory with the annotation database (default: the "data" folder of the installation, '
                             'or the environment variable SEMANTICSBML_DATA)')
    parser.add_argument('--progress', action='store_true', help='show progress messages while loading the database')
    sub = parser.add_subparsers(dest='command', metavar='COMMAND')

    p = sub.add_parser('check', help='check annotated models for completeness and consistency')
    p.add_argument('files', nargs='+', metavar='MODEL.xml')
    p.add_argument('-v', '--verbose', action='store_true', help='also list the checks that found nothing')
    p.add_argument('-q', '--quiet', action='store_true', help='no output, only the exit code (1: warnings, 2: errors)')
    p.set_defaults(func=cmd_check)

    p = sub.add_parser('view', help='draw the reaction network of a model (needs graphviz)')
    p.add_argument('file', metavar='MODEL.xml')
    p.add_argument('-o', '--outfile', help='output file; the ending decides the format (.svg .png .pdf .dot). '
                                           'Without -o the graphviz dot code is printed')
    p.add_argument('--ids', action='store_true', help='label nodes with the ids of the model instead of the names of the annotations')
    p.add_argument('--black-and-white', action='store_true')
    p.add_argument('--max-participations', type=int, default=1000, metavar='N',
                   help='leave out species that take part in more than N reactions (e.g. ATP, water)')
    p.add_argument('--max-species', type=int, default=500, metavar='N', help='refuse to draw models with more species (default 500)')
    p.set_defaults(func=cmd_view)

    p = sub.add_parser('annotations', help='list the annotations of a model')
    p.add_argument('file', metavar='MODEL.xml')
    p.add_argument('--missing', action='store_true', help='show only elements without annotation')
    p.add_argument('--annotated', action='store_true', help='show only elements with annotation')
    p.add_argument('--tsv', action='store_true', help='tab separated output')
    p.set_defaults(func=cmd_annotations)

    p = sub.add_parser('search', help='search the annotation database for a name')
    p.add_argument('query')
    p.add_argument('--exact', action='store_true', help='only exact matches of the name')
    p.add_argument('--precision', type=float, default=0.8, help='similarity of names, between 0 and 1 (default 0.8)')
    p.add_argument('--type', choices=['species', 'reaction', 'compartment', 'parameter', 'rule', 'event', 'model'],
                   help='only databases that are appropriate for this kind of SBML element')
    p.set_defaults(func=cmd_search)

    p = sub.add_parser('annotate', help='add or remove annotations',
                       description='Add or remove annotations. The steps are carried out in this order: '
                                   'clear, remove, copy-from, add, auto.')
    p.add_argument('file', metavar='MODEL.xml')
    p.add_argument('-o', '--outfile', help='file for the annotated model (default: print to screen)')
    p.add_argument('--auto', action='store_true',
                   help='annotate all elements without annotation whose name is found exactly once in the database')
    p.add_argument('--add', nargs='+', action='append', metavar='X',
                   help='--add ELEMENT DATABASE IDENTIFIER [QUALIFIER]   e.g. --add GLCi ChEBI CHEBI:17234 is')
    p.add_argument('--remove', nargs=3, action='append', metavar=('ELEMENT', 'DATABASE', 'IDENTIFIER'))
    p.add_argument('--clear', action='append', metavar='ELEMENT', help='remove all annotations of this element')
    p.add_argument('--clear-all', action='store_true', help='remove all annotations of the model')
    p.add_argument('--copy-from', metavar='OTHER.xml', help='copy annotations from elements with the same id in another model')
    p.add_argument('--overwrite', action='store_true', help='with --copy-from: replace existing annotations')
    p.add_argument('-v', '--verbose', action='store_true')
    p.set_defaults(func=cmd_annotate)

    p = sub.add_parser('merge', help='merge annotated models',
                       description='Merge two or more models. Elements with matching annotations are joined. '
                                   'Where joined elements differ (values, formulas, units), the model that comes '
                                   'first in the list wins. A report is printed.')
    p.add_argument('files', nargs='+', metavar='MODEL.xml')
    p.add_argument('-o', '--outfile', help='file for the merged model (default: print to screen)')
    p.add_argument('--name', help='name of the merged model')
    p.add_argument('--unmatched', action='store_true', help='also list the elements that were not matched')
    p.add_argument('-v', '--verbose', action='store_true')
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser('similarity', help='compare whole models by their annotations',
                       description='Similarity (0 to 1) between annotated models. Without --query a table of all '
                                   'pairs is printed; with --query the other models are ranked by their similarity '
                                   'to that model. How much an annotation counts depends on how often it occurs '
                                   'in the given set of models.')
    p.add_argument('files', nargs='*', metavar='MODEL.xml')
    p.add_argument('--dir', metavar='FOLDER', help='use all .xml files in this folder')
    p.add_argument('--query', metavar='MODEL.xml', help='rank the other models by similarity to this one')
    p.add_argument('--cluster', action='store_true', help='order the table by hierarchical clustering')
    p.add_argument('--tsv', action='store_true', help='tab separated output')
    p.set_defaults(func=cmd_similarity)

    p = sub.add_parser('match', help='find corresponding elements of models by the similarity of their annotations',
                       description='Find corresponding elements in two or more models. Unlike "merge", which needs '
                                   'identical annotations, this also uses relations between database entries '
                                   '(e.g. D-glucose is a glucose) and the neighbourhood of an element in the network.')
    p.add_argument('files', nargs='+', metavar='MODEL.xml')
    p.add_argument('--min', type=float, default=0.5, metavar='S', help='smallest similarity that is reported (default 0.5)')
    p.add_argument('--local', action='store_true',
                   help='compare only the annotations of the elements themselves, not those of their neighbours in the network')
    p.add_argument('--neighbours', type=int, default=1, metavar='N',
                   help='number of relation steps between database entries that are followed (default 1)')
    p.add_argument('--tsv', action='store_true', help='tab separated output')
    p.set_defaults(func=cmd_match)

    p = sub.add_parser('database', help='show which version of the annotation database is in use')
    p.set_defaults(func=cmd_database)

    p = sub.add_parser('update-database', help='download current ChEBI, Gene Ontology, SBO, KEGG names and '
                                               'identifiers.org namespaces and update the annotation database',
                       description='Updates the annotation database in the data folder (about 100 MB of downloads, '
                                   'several minutes). Names and relations of ChEBI, Gene Ontology and SBO are '
                                   'replaced by the current ones and new entries are added. From KEGG only the '
                                   'lists of compound and reaction names are fetched (two requests); KEGG is free '
                                   'for academic use only, use --no-kegg otherwise.')
    p.add_argument('--no-kegg', action='store_true', help='do not contact the KEGG web service')
    p.add_argument('--all-chebi', action='store_true',
                   help='take over all ChEBI entries (about 220000) instead of the manually curated ones and those '
                        'with a KEGG identifier; makes every command slower to start')
    p.add_argument('--only', nargs='+', metavar='PART', choices=['chebi', 'go', 'sbo', 'kegg', 'registry'],
                   help='update only these parts: chebi go sbo kegg registry')
    p.add_argument('--download-dir', metavar='DIR', help='folder for the downloaded files (default: data/downloads)')
    p.add_argument('--use-existing', action='store_true', help='use files in the download folder instead of downloading again')
    p.set_defaults(func=cmd_update_database)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 0
    if args.data:
        os.environ['SEMANTICSBML_DATA'] = args.data
    if not args.progress:
        os.environ.setdefault('SEMANTICSBML_QUIET', '1')
    try:
        return args.func(args) or 0
    except CliError as e:
        sys.stderr.write('Error: %s\n' % e)
        return 2
    except BrokenPipeError:
        return 0
    except Exception as e:
        # errors of the library (AnnotationError, MergeError, DatabaseError ...)
        if os.environ.get('SEMANTICSBML_DEBUG'):
            raise
        sys.stderr.write('Error: %s\n(set the environment variable SEMANTICSBML_DEBUG=1 for details)\n' % e)
        return 2


if __name__ == '__main__':
    sys.exit(main())
