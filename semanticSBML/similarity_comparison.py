#!/usr/bin/env python

import math
import semanticSBML.similarity
import semanticSBML.annotate

k0=0
k1=1

class similarityDictionary:
    def __init__(self, esm, names):
        self.esm = esm
        [self.names0, self.names1] = names
    def get_similarity(self, id1, id2):
        return self.esm.get_similarity(k0, k1, id1, id2)
    def get_matrix(self, decorated=False, as_str=False, normalized=False):
        # get max entries
        if normalized:
            max_values0 = {}
            max_values1 = {}
            for name0 in self.names0:
                maxval = max([self.get_similarity(name0, name1) for name1 in self.names1])
                if maxval == 0.:
                    maxval = 1.
                max_values0[name0] = math.sqrt(maxval)
            for name1 in self.names1:
                maxval = max([self.get_similarity(name0, name1) for name0 in self.names0])
                if maxval == 0.:
                    maxval = 1.
                max_values1[name1] = math.sqrt(maxval)
        m = []
        if decorated:
            m.append(['']+self.names1)
        for name0 in self.names0:
            l = []
            if decorated:
                l.append(name0)
            for name1 in self.names1:
                sim = self.get_similarity(name0, name1)
                if normalized:
                    sim = sim / max_values0[name0] / max_values1[name1]
                l.append(sim)
            m.append(l)
        if as_str:
            for i1 in range(len(m)):
                for i2 in range(len(m[i1])):
                    m[i1][i2] = str(m[i1][i2])
        return m
    def diff(self, other):
        new = semanticSBML.similarity.ElementSimilarityMatrix()
        for name0 in self.names0:
            for name1 in self.names1:
                new_sim = self.get_similarity(name0, name1) - other.get_similarity(name0, name1)
                if new_sim != 0.:
                    new.add_similarity(k0, k1, name0, name1, new_sim)
        return similarityDictionary(new, [self.names0, self.names1])
    def get_greedy_matches(self):
        new = semanticSBML.similarity.ElementSimilarityMatrix()
        black0 = []
        black1 = []
        while True:
            pop = self.esm.pop_similarity()
            if pop is None:
                break
            [sim, foo, bar, id0, id1] = pop
            if id0 not in black0 and id1 not in black1:
                new.add_similarity(k0, k1, id0, id1, 1.)
                black0.append(id0)
                black1.append(id1)
        self.esm.rebuild_stack()
        return similarityDictionary(new, [self.names0, self.names1])
    def get_size(self):
        size = 0
        while True:
            pop = self.esm.pop_similarity()
            if pop is None:
                break
            [sim, foo, bar, id0, id1] = pop
            size += abs(sim)
        self.esm.rebuild_stack()
        return size
        
        
        

class similarityComparison:
    def __init__(self, SBMLDocuments, max_distance_neighbours=1):
        self.SBMLDocuments = SBMLDocuments
        self.meas = [semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True) for doc in self.SBMLDocuments]
        self.names = [[ea.id for ea in reversed(mea.getElementAnnotations())] for mea in self.meas]
        self.sigma = semanticSBML.similarity.SigmaMatrix(self.SBMLDocuments, [k0, k1], max_distance_neighbours)
        self.global_sigma = semanticSBML.similarity.GlobalSigma(self.SBMLDocuments, [k0, k1], self.meas, self.sigma)
        self.global_sigma_sim_prop = semanticSBML.similarity.GlobalSigmaSimProp(self.SBMLDocuments, [k0, k1], self.meas, self.sigma)

        self.local_dict = similarityDictionary(self.sigma.esm, self.names)
        self.global_dict = similarityDictionary(self.global_sigma.global_esm, self.names)
        self.global_dict_sim_prop = similarityDictionary(self.global_sigma_sim_prop.global_esm, self.names)
        self.diff_dict = self.global_dict.diff(self.local_dict)

        self.local_match = self.local_dict.get_greedy_matches()
        self.global_match = self.global_dict.get_greedy_matches()
        self.global_match_sim_prop = self.global_dict_sim_prop.get_greedy_matches()
        self.diff_match = self.global_match.diff(self.local_match)
        

def triangle():
    import libsbml
    sc = similarityComparison([libsbml.readSBML('/home/schulzma/Desktop/triangle.xml')]*2)
    gvs = sc.global_sigma.global_vector[0]
    all_indices = []
    maxes = []
    for gv in gvs:
        if len([x for x in gv.tolist()[0] if x != 0.])>0:
            vec = gv.tolist()[0]
            maxes.append(max(vec))
            for i in range(len(vec)):
                if vec[i] != 0. and i not in all_indices:
                    all_indices.append(i)
    all_indices = sorted(all_indices)
    all_max = max(maxes)
    for gv in gvs:
        if len([x for x in gv.tolist()[0] if x != 0.])>0:
            vec = gv.tolist()[0]
            print([str(int(vec[i]/all_max*256)) for i in all_indices])
    print([ea.id for ea in sc.meas[0].getElementAnnotations()])
            

def kegg():
    import os, libsbml, re
    import libSBAnnotation.database_new
    d=libSBAnnotation.database_new.Meta_DB()
    reaction_re = re.compile('R\\d{5}')
    ec_re = re.compile('E\\d+ \\d+ \\d+ \\d+')
    folder = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/hsa/'
    files = sorted([x for x in os.listdir(folder) if x.endswith('.xml')])
    none_found = []
    name2annotation = {}
    for f in files:
        print(f)
        doc = libsbml.readSBML(folder+f)
        mea = semanticSBML.annotate.ModelElementsAnnotations(doc.getModel())
        for ea in mea.getElementAnnotations():
            if ea.type in 'model compartment'.split():
                continue
            if reaction_re.match(ea.id):
                ea.addAnnotation(semanticSBML.annotate.Annotation(db='KEGG Reaction', id=ea.id, qualifier_type='bio', qualifier='is'))
                raid = d.database_instance.__getAbstractIdByIdDict__({'KEGG Reaction':ea.id})
                paids = [x[0] for x in d.database_instance.takes_part.findAllRows({'reaction_id':raid})]
                for paid in paids:
                    ao = d.database_instance.getAnnotationObjectByAbstractId(paid)
                    if ao is None:
                        continue
                    ao = ao.toMiriamCompliantDatabase2id()
                    names = d.database_instance.__getAllNamesByAbstractId__(paid)
                    for name in names:
                        for char in '()':
                            name = name.replace(char, ' ')
                        name2annotation[name] = ao
            else:
                name = ea.id
                name = name.replace('_space_',' ')
                name = name.replace('_plus_','+')
                name = name.replace('_minus_','-')
                for s in 'alpha beta gamma'.split():
                    name = name.replace("_%s_"%s, s)
                name = name.strip()

                if ec_re.match(name):
                    name = name.replace(' ', '.')[1:]
                    ea.addAnnotation(semanticSBML.annotate.Annotation(db='Enzyme Nomenclature', id=name, qualifier_type='bio', qualifier='is'))
                else:
                    import libSBAnnotation.mydifflib
                    new_names = libSBAnnotation.mydifflib.get_close_matches(name.lower(), list(name2annotation.keys()), 1, .9)
                    if len(new_names) > 0:
                        ao = name2annotation[new_names[0]]
                        for db in ao:
                            try:
                                ea.addAnnotation(semanticSBML.annotate.Annotation(db=db, id=ao[db], qualifier_type='bio', qualifier='is'))
                            except semanticSBML.annotate.AnnotationCreationError:
                                continue
                    else:
                        if name not in none_found:
                            none_found.append(name)
        break
    print(none_found)

def kegg2():
    import os, libsbml, re
    import libSBAnnotation.database_new
    d=libSBAnnotation.database_new.Meta_DB()
    reaction_re = re.compile('R\\d{5}')
    ec_re = re.compile('E\\d+ (\\d+|-) (\\d+|-) (\\d+|-)')
    folder = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/hsa/'
    files = sorted([x for x in os.listdir(folder) if x.endswith('.xml') and not x.endswith('.xml.xml')])
    none_found = []
    s2aid = {}
    s2finalaid = {}
    for f in files:
        print(f)
        doc = libsbml.readSBML(folder+f)
        mea = semanticSBML.annotate.ModelElementsAnnotations(doc.getModel())
        for ea in mea.getElementAnnotations():
            if ea.type in 'model compartment'.split():
                continue
            if reaction_re.match(ea.id):
                #ea.addAnnotation(semanticSBML.annotate.Annotation(db='KEGG Reaction', id=ea.id, qualifier_type='bio', qualifier='is'))
                raid = d.database_instance.__getAbstractIdByIdDict__({'KEGG Reaction':ea.id})
                paids = [x[0] for x in d.database_instance.takes_part.findAllRows({'reaction_id':raid})]
                participants = [s.getSpecies() for s in ea.libsbml_element.getListOfReactants()] + [s.getSpecies() for s in ea.libsbml_element.getListOfProducts()]
                for participant in participants:
                    if participant not in s2aid:
                        s2aid[participant] = set(paids)
                    else:
                        s2aid[participant] = s2aid[participant].intersection(set(paids))
            else:
                name = ea.id.replace('_space_',' ').replace('_plus_','+').replace('_minus_','-')
                for s in 'alpha beta gamma'.split():
                    name = name.replace("_%s_"%s, s)
                name = name.strip()

                if ec_re.match(name):
                    name = name.replace(' ', '.')[1:]
                    #ea.addAnnotation(semanticSBML.annotate.Annotation(db='Enzyme Nomenclature', id=name, qualifier_type='bio', qualifier='is'))
                else:
                    if name in s2finalaid:
                        continue
                    if ea.id not in s2aid:
                        none_found.append(name)
                        continue
                    aid = None
                    aids = list(s2aid[ea.id])
                    if len(aids) == 1:
                        aid = aids[0]
                    else:
                        bestaid = None
                        bestname = None
                        bestratio = 0.
                        for thisaid in aids:
                            for othername in d.database_instance.__getAllNamesByAbstractId__(thisaid):
                                for brace in '()':
                                    othername = othername.replace(brace, ' ')
                                othername = othername.strip()
                                if othername == name:
                                    aid = thisaid
                                    s2finalaid[name] = thisaid
                        if aid is not None:
                            ao = d.database_instance.getAnnotationObjectByAbstractId(aid)
                            if ao is None:
                                continue
                            ao = ao.toMiriamCompliantDatabase2id()
                            for db in ao:
                                try:
                                    ea.addAnnotation(semanticSBML.annotate.Annotation(db=db, id=ao[db], qualifier_type='bio', qualifier='is'))
                                except semanticSBML.annotate.AnnotationCreationError:
                                    continue
                        else:
                            if name not in none_found:
                                none_found.append(name)
        libsbml.writeSBML(doc, folder+f+'.xml')
    print(none_found)

def kegg_new():
    import os, libsbml, re
    import libSBAnnotation.database_new
    d=libSBAnnotation.database_new.Meta_DB()
    file_re = re.compile('sbml_hsa\\d{5}\\.xml')
    reaction_re = re.compile('R\\d{5}')
    compound_re = re.compile('C\\d{5}')
    ec_re = re.compile('\\d+\\.\\d+\\.(\\d+|-)\\.(\\d+|-)')
    folder = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/hsa_new/'
    files = sorted([x for x in os.listdir(folder) if file_re.match(x)])
    for f in files:
        print(f)
        doc = libsbml.readSBML(folder+f)
        mea = semanticSBML.annotate.ModelElementsAnnotations(doc.getModel())
        for ea in mea.getElementAnnotations():
            if ea.type in 'model compartment'.split():
                continue
            if reaction_re.match(ea.id):
                for rid in ea.id.split(' '):
                    ea.addAnnotation(semanticSBML.annotate.Annotation(db='KEGG Reaction', id=rid, qualifier_type='bio', qualifier='is'))
            else:
                # compound
                for (regex,db) in [(compound_re, 'KEGG Compound'), (ec_re, 'Enzyme Nomenclature')]:
                    s = regex.search(ea.libsbml_element.getName())
                    if s is not None:
                        _id = s.group(0)
                        ea.addAnnotation(semanticSBML.annotate.Annotation(db=db, id=_id, qualifier_type='bio', qualifier='is'))
        libsbml.writeSBML(doc, folder+f+'.xml')

def main():
    import libsbml
    models = [9,11]#[9, 10, 11, 14, 26, 27, 84]
    for model0 in models:
        for model1 in models:
            if model1 > model0:
                sbml_documents = [libsbml.readSBML('/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/BIOMD'+str(f).rjust(10, '0')+'.xml') for f in [model0, model1]]
                print('======', model0, model1)
                sc = similarityComparison(sbml_documents)
