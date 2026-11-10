# Runs the same basic operations under the original (Python 2) and the port (Python 3)
# usage: python ref.py MODELDIR OUTDIR section [section ...]
from __future__ import print_function
import sys, os, json, time, traceback, warnings
warnings.simplefilter('ignore')
import libsbml

MODELDIR, OUTDIR = sys.argv[1], sys.argv[2]
SECTIONS = sys.argv[3:]
M61 = os.path.join(MODELDIR, 'BIOMD0000000061.xml')
M64 = os.path.join(MODELDIR, 'BIOMD0000000064.xml')
M239 = os.path.join(MODELDIR, 'BIOMD0000000239.xml')


def save(name, obj):
    with open(os.path.join(OUTDIR, name + '.json'), 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=True)


def sec_annotations():
    from semanticSBML import annotate
    out = {}
    for path in (M61, M64, M239):
        doc = libsbml.readSBML(path)
        rows = []
        for ea in annotate.ModelElementsAnnotations(doc.getModel(), True).getElementAnnotations():
            for a in ea.getAnnotations():
                try:
                    name = a.getName()
                except Exception as e:
                    name = 'ERR'
                rows.append([str(ea.id), str(a.qualifier), str(a.db), str(a.id), '%s' % (name,)])
        out[os.path.basename(path)] = sorted(rows)
    return out


def sec_merge():
    from semanticSBML import merge
    docs = [libsbml.readSBML(M61), libsbml.readSBML(M64)]
    merger = merge.Merger(docs, lambda *a: None, lambda *a: None)
    matches = []
    for t in merger.getTuples():
        matches.append(sorted([str(e.merge_model.mid), str(e.id), str(e.name)] for e in t.getElements()))
    out = {'matches': sorted(matches)}
    merger.resolveUsingModel(merger.merge_models)
    doc = merger.createSBMLmodel('merged')
    libsbml.writeSBML(doc, os.path.join(OUTDIR, 'merged.xml'))
    m = doc.getModel()
    out['counts'] = {'species': m.getNumSpecies(), 'reactions': m.getNumReactions(),
                     'parameters': m.getNumParameters(), 'compartments': m.getNumCompartments(),
                     'rules': m.getNumRules(), 'units': m.getNumUnitDefinitions(),
                     'functions': m.getNumFunctionDefinitions()}
    out['species'] = sorted(str(m.getSpecies(i).getName() or m.getSpecies(i).getId()) for i in range(m.getNumSpecies()))
    doc2 = libsbml.readSBML(os.path.join(OUTDIR, 'merged.xml'))
    doc2.checkConsistency()
    errs = []
    for i in range(doc2.getNumErrors()):
        e = doc2.getError(i)
        if e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            errs.append(e.getMessage().strip().split('\n')[-1].strip())
    out['validation_errors'] = errs
    return out


def sec_match():
    from semanticSBML import similarity
    docs = [libsbml.readSBML(M61), libsbml.readSBML(M64)]
    out = {}
    for label, local in (('local', True), ('propagation', False)):
        t = time.time()
        matcher = similarity.ModelMatcher(docs, [1, 2], 1, True)
        if local:
            matcher.pick_propagation_type('local')
        res = []
        for match in matcher.match(0.5):
            ids = [d['element_id'] for d in match]
            mids = [d['mid'] for d in match]
            s = matcher.global_esm.get_similarity(mids[0], mids[1], ids[0], ids[1]) if len(match) == 2 else -1
            res.append([round(float(s), 6)] + sorted('%s:%s' % (d['mid'], d['element_id']) for d in match))
        out[label] = sorted(res)
        out[label + '_seconds'] = round(time.time() - t, 1)
        save('match', out)
    return out


def sec_similarity():
    from semanticSBML import annotate
    from libModelSimilarity import cluster
    coll = cluster.DatabaseAnnotation()
    for path in (M61, M64, M239):
        doc = libsbml.readSBML(path)
        mea = annotate.ModelElementsAnnotations(doc.getModel(), True)
        coll.model_annotations.append(cluster.ModelAnnotation(mea, os.path.basename(path)))
    if hasattr(coll, 'update_information_content'):
        coll.update_information_content()
    ms = coll.model_annotations
    out = {}
    for i in range(3):
        for j in range(i, 3):
            out['%s~%s' % (ms[i].model_name, ms[j].model_name)] = round(float(ms[i].similarity(ms[j])), 6)
    return out


def sec_check():
    from semanticSBML.Model import Model
    from semanticSBML.sbmlcheck import Checker
    out = {}
    for path in (M61, M64):
        rows = []
        for r in Checker(Model(path)).run():
            rows.append([str(r.name), r.errors, r.warnings, sorted(str(m) for m in r.error_messages),
                         sorted(str(m) for m in r.warning_messages)])
        out[os.path.basename(path)] = sorted(rows)
    return out


for s in SECTIONS:
    t = time.time()
    try:
        save(s, globals()['sec_' + s]())
        print(s, 'ok', round(time.time() - t, 1), 's')
    except Exception:
        print(s, 'FAILED')
        traceback.print_exc()
    sys.stdout.flush()
