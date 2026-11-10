#!/usr/bin/env python

import numpy
import os
import scipy
import libsbml
import semanticSBML.annotate
import semanticSBML.similarity
import libSBAnnotation.cofactors
import libSBAnnotation.database_new
import libSBAnnotation.config
database = libSBAnnotation.database_new.Meta_DB()

def vis_aid(aid):
    ao = database.database_instance.getAnnotationObjectByAbstractId(aid)
    name = database.database_instance.__getNameByAbstractId__(aid)
    return name

class SparseVector:
    def __init__(self):
        self.d = {}
        self.l = None
    def fromSAID2val(self, SAID2val, blacklist=[]):
        self.d = {}
        for said in SAID2val:
            if said in blacklist:
                continue
            try:
                self.d[int(said)] = SAID2val[said]
            except ValueError:
                pass
    def length(self):
        if self.l is None:
            l = 0.
            for val in list(self.d.values()):
                l += numpy.power(val,2)
            self.l = numpy.sqrt(l)
        return self.l
    def angle(self, other):
        cross_prod = 0
        for key in list(self.d.keys()):
            if key in other.d:
                cross_prod += self.d[key]*other.d[key]
        try:
            return cross_prod / self.length() / other.length()
        except ZeroDivisionError:
            return 0.
    def int_keys(self):
        return list(self.d.keys())
    def tostr(self):
        vals = list(set(self.d.values()))
        vals.sort()
        vals.reverse()
        lines = []
        for val in vals:
            lines.append("%s:%s"%(val,[vis_aid(aid) for aid in self.d if self.d[aid] == val]))
        return '\n'.join(lines)
    def __repr__(self):
        return self.tostr()
    def __str__(self):
        return self.tostr()

class AnnotationLookaheadDatabase:
    __shared_state = {}
    def __init__(self, limited_propagation=None):
        self.__dict__ = self.__shared_state
        if len(list(self.__dict__.keys())) == 0:
            self.limited_propagation = limited_propagation
            self.documents = []
            self.filenames = []
            self.vector_and_annotation = []
            self.hash = None
            self.blacklist = []#[str(x) for x in libSBAnnotation.cofactors.Cofactors().get_cofactor_aids()]
            self.filename = libSBAnnotation.config.Config().getpath('database','database_dir')+os.sep+'annotation_lookahead.txt'
            if not self.load(self.filename):
                self.load_from_path('/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/')
                self.dump(self.filename)
    def delete(self):
        self.__class__.__shared_state = {}
        del self
    def create_hash(self):
        self.hash = {}
        for i in range(len(self.vector_and_annotation)):
            v = self.vector_and_annotation[i][0]
            for key in v.int_keys():
                if key not in self.hash:
                    self.hash[key] = []
                self.hash[key].append(i)
    def find_annotations(self, vector):
        if self.hash is None:
            self.create_hash()
        index_list = set([])
        for key in vector.int_keys():
            if key in self.hash:
                for index in self.hash[key]:
                    index_list.add(index)
        anno2val = {}
        for index in list(index_list):
            [other_vector, anno] = self.vector_and_annotation[index]
            sim = vector.angle(other_vector)
            if anno in anno2val:
                anno2val[anno] = max(sim, anno2val[anno])
            else:
                anno2val[anno] = sim
        return anno2val
    def load_from_path(self, directory):
        files = sorted(os.listdir(directory))
        for f in files:
            if f.endswith('.xml'):
                self.documents.append(libsbml.readSBML(directory+os.sep+f))
                self.filenames.append(f)
        for i in range(len(self.documents)):
            doc = self.documents[i]
            meas = [semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True)]
            sigma = semanticSBML.similarity.SigmaMatrix([doc], [self.filenames[i]], 0, use_full_database=False)
            global_sigma_anno_prop = semanticSBML.similarity.GlobalSigma([doc], [self.filenames[i]], meas, sigma, self.limited_propagation)
            elements = global_sigma_anno_prop._get_model_elements(doc)
            doc_global_vector = global_sigma_anno_prop.global_vector[0]
            saids = sigma.keys
            for element,global_vector in zip(elements,doc_global_vector):
                local_vector = sigma._element_annotation_to_vector(semanticSBML.annotate.ElementAnnotations(element, True))
                SAID2val = {}
                for j in range(len(sigma.keys)):
                    key = sigma.keys[j]
                    lv = local_vector[0,j]
                    gv = global_vector[0,j]
                    if gv != 0. and lv == 0.:
                        SAID2val[key] = gv
                vector = SparseVector()
                vector.fromSAID2val(SAID2val, self.blacklist)
                for j in range(len(sigma.keys)):
                    lv = local_vector[0,j]
                    key = sigma.keys[j]
                    if lv != 0.:
                        try:
                            self.vector_and_annotation.append([vector, int(key)])
                        except ValueError:
                            pass
            print('Done with',self.filenames[i])
        self.create_hash()
    def dump(self, filename):
        import pickle
        f = open(filename, 'w')
        pickle.dump([self.vector_and_annotation, self.hash], f)
        f.close()
    def load(self, filename):
        import pickle
        try:
            f = open(filename, 'r')
        except:
            return False
        loaded = pickle.load(f)
        self.vector_and_annotation = loaded[0]
        self.hash = loaded[1]
        f.close()
        return True

class AnnotationLookahead:
    def __init__(self, document, model_element_annotations=None, annotation_lookahead_database=None, limited_propagation=None):
        self.filenames = ['blabla']
        self.document = document
        self.mea = model_element_annotations
        self.limited_propagation = limited_propagation
        if self.mea is None:
            self.mea = semanticSBML.annotate.ModelElementsAnnotations(document.getModel(), True)
        if annotation_lookahead_database is None:
            annotation_lookahead_database = AnnotationLookaheadDatabase(limited_propagation=self.limited_propagation)
        self.sigma = semanticSBML.similarity.SigmaMatrix([self.document], self.filenames, 0, use_full_database=False)
        self.global_sigma_anno_prop = semanticSBML.similarity.GlobalSigma([self.document], self.filenames, [self.mea], self.sigma, limited_propagation=self.limited_propagation)
        self.element2aid2value = {}

        elements = self.global_sigma_anno_prop._get_model_elements(self.document)
        doc_global_vector = self.global_sigma_anno_prop.global_vector[0]
        self.element_id2vector = {}
        for element,global_vector in zip(elements,doc_global_vector):
            self.element2aid2value[element.getId()] = {}
            SAID2val = {}
            for j in range(len(self.sigma.keys)):
                if global_vector[0,j] > 10e-10:
                    SAID2val[self.sigma.keys[j]] = global_vector[0,j]
            vector = SparseVector()
            vector.fromSAID2val(SAID2val, annotation_lookahead_database.blacklist)
            self.element2aid2value[element.getId()] = annotation_lookahead_database.find_annotations(vector)
            self.element_id2vector[element.getId()] = vector
    def aid2type(self,aid):
        aodk = list(database.getAnnotationObjectByAbstractId(aid).dictionary.keys())
        if 'kegg_reaction' in aodk or 'ec' in aodk:
            return ['reaction']
        elif 'kegg_compound' in aodk or 'chebi' in aodk:
            return ['species']
        elif 'go' in aodk:
            return ['compartment', 'species', 'reaction']
        return []
        
    def predict_aid(self, element_id, element_type):
        try:
            aid2value = self.element2aid2value[element_id]
        except KeyError:
            # maybe a model has been inserted
            return None
        return [[aid, aid2value[aid]] for aid in aid2value if element_type in self.aid2type(aid) and not numpy.isnan(aid2value[aid])]

    def predict_grouped_db_id(self, element_id, element_type):
        aid_values = self.predict_aid(element_id, element_type)
        # sorting out double aids
        new_aids = []
        for aid_value in aid_values:
            if aid_value[0] not in new_aids:
                new_aids.append(aid_value[0])
        # converting each aid to Miriam compliant database 2 id
        grouped = []
        for aid in new_aids:
            ao = database.database_instance.getAnnotationObjectByAbstractId(aid)
            if ao is not None:
                grouped.append(ao.toMiriamCompliantDatabase2id())
        return grouped



def main():
    lp = 1
    import libSBAnnotation.database_new
    database = libSBAnnotation.database_new.Meta_DB()

    ald = AnnotationLookaheadDatabase(limited_propagation=lp)
    #ald.load_from_path('/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/')
    
    sbmlfilename = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/BIOMD0000000061.xml'
    doc = libsbml.readSBML(sbmlfilename)
    mea = semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True)
    id2urns = {}
    id2annotations = {}
    for ea in mea.getElementAnnotations():
        if ea.libsbml_element.getId().startswith('vPFK'):
            id2urns[ea.libsbml_element.getId()] = []
            id2annotations[ea.libsbml_element.getId()] = []
            for anno in ea.getAnnotations():
                id2urns[ea.libsbml_element.getId()].append(anno.db+':'+anno.id)
                id2annotations[ea.libsbml_element.getId()].append(anno)
            # unset the annotations
            ea.annotations = []
    al = AnnotationLookahead(doc, mea, ald, limited_propagation=lp)
    el_aid_val = []
    positions = []
    for ea in mea.getElementAnnotations():
        if ea.libsbml_element.getId() in id2urns:
            aid_and_value = al.predict_aid(ea.libsbml_element.getId(), ea.type)
            if aid_and_value is not None:
                print(ea.libsbml_element.getId(), id2urns[ea.libsbml_element.getId()])
                sorted_aid_and_value = sorted(aid_and_value, key=lambda l:l[1], reverse=True)[:20]
                for aid,value in sorted_aid_and_value:
                    print(aid, value , database.getAnnotationObjectByAbstractId(aid), database.database_instance.__getNameByAbstractId__(aid))
                    for vector in [x[0] for x in ald.vector_and_annotation if int(x[1])==aid ]:
                        print(' ',vector.tostr())
                    print() 
                    el_aid_val.append([ea.libsbml_element.getId(), aid, value])
                print()
                these_aids = []
                for anno in id2annotations[ea.libsbml_element.getId()]:
                    this_aid = database.database_instance.__getAbstractIdByIdDict__({anno.db:anno.id})
                    these_aids.append(this_aid)
                    these_vectors = [x[0] for x in ald.vector_and_annotation if x[1]==this_aid]
                    max_sim = 0.0
                    for other_vector in these_vectors:
                        this_sim = al.element_id2vector[ea.libsbml_element.getId()].angle(other_vector)
                        max_sim = max(max_sim, this_sim)
                    #print anno.db, anno.id, this_aid, max_sim
                #raw_input("?")
                this_position = 101
                for i in range(len(sorted_aid_and_value)):
                    [aid, value] = sorted_aid_and_value[i]
                    if aid in these_aids:
                        this_position = i
                        break
                positions.append(this_position)
    print(positions)
    print(sum(positions))
    print(float(sum(positions))/float(len(positions)))






    sbmlfilename = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/BIOMD0000000061.xml'
    doc = libsbml.readSBML(sbmlfilename)
    mea = semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True)
    id2urns = {}
    id2annotations = {}
    for ea in mea.getElementAnnotations():
        if not ea.libsbml_element.getId().startswith('v'):
            id2urns[ea.libsbml_element.getId()] = []
            id2annotations[ea.libsbml_element.getId()] = []
            for anno in ea.getAnnotations():
                id2urns[ea.libsbml_element.getId()].append(anno.db+':'+anno.id)
                id2annotations[ea.libsbml_element.getId()].append(anno)
            # unset the annotations
            ea.annotations = []
    al = AnnotationLookahead(doc, mea, ald, limited_propagation=lp)
    el_aid_val = []
    positions = []
    for ea in mea.getElementAnnotations():
        if ea.libsbml_element.getId() in id2urns:
            aid_and_value = al.predict_aid(ea.libsbml_element.getId(), ea.type)
            if aid_and_value is not None:
                #print ea.libsbml_element.getId(), id2urns[ea.libsbml_element.getId()]
                sorted_aid_and_value = sorted(aid_and_value, key=lambda l:l[1], reverse=True)[:100]
                for aid,value in sorted_aid_and_value:
                    #print aid, value ,database.getAnnotationObjectByAbstractId(aid), database.database_instance.__getNameByAbstractId__(aid)
                    el_aid_val.append([ea.libsbml_element.getId(), aid, value])
                #print
                these_aids = []
                for anno in id2annotations[ea.libsbml_element.getId()]:
                    this_aid = database.database_instance.__getAbstractIdByIdDict__({anno.db:anno.id})
                    these_aids.append(this_aid)
                    these_vectors = [x[0] for x in ald.vector_and_annotation if x[1]==this_aid]
                    max_sim = 0.0
                    for other_vector in these_vectors:
                        this_sim = al.element_id2vector[ea.libsbml_element.getId()].angle(other_vector)
                        max_sim = max(max_sim, this_sim)
                    #print anno.db, anno.id, this_aid, max_sim
                #raw_input("?")
                this_position = 101
                for i in range(len(sorted_aid_and_value)):
                    [aid, value] = sorted_aid_and_value[i]
                    if aid in these_aids:
                        this_position = i
                        break
                positions.append(this_position)
    print(positions)
    print(sum(positions))
    print(float(sum(positions))/float(len(positions)))


def main2():
    ald = AnnotationLookaheadDatabase()

if __name__=="__main__":
    main2()
