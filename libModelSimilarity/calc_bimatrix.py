#!/usr/bin/env python

import sys, os, semanticSBML.annotate, libsbml, libSBAnnotation.database_new, scipy, scipy.sparse
from . import cluster, fake_annotations

cluster.similarity_measure_version = 4
database_new_instance = libSBAnnotation.database_new.Meta_DB()
add_qsms = False
use_qsms = False
result_number = 3#sys.argv[2]#'3'
bimatrix_filename  = os.path.join( os.path.expanduser('~'), '.semanticSBML', 'bimatrix.csv' )
simmatrix_filename = os.path.join( os.path.expanduser('~'), '.semanticSBML', 'simmatrix.csv' )

use_libSBAnnotation = False
use_similarity_in_euklidian = False
use_semantic_density = False
if int(result_number) > 1:
    use_libSBAnnotation = True
if int(result_number) > 2:
    use_similarity_in_euklidian = True
if int(result_number) > 3:
    use_semantic_density = True

global annotation_count
annotation_count = {}

class mycsv:
    def __init__(self):
        self.matrix = []
    def load_from_dictionary(self, dictionary):
        self.primary_keys = list(dictionary.keys())
        self.secondary_keys = set([])
        [[self.secondary_keys.add(key2) for key2 in dictionary[key1]] for key1 in self.primary_keys]
        self.secondary_keys = list(self.secondary_keys)
        self.secondary_keys.sort()
        self.matrix = []
        self.matrix.append( [""]+self.secondary_keys )
        for key1 in self.primary_keys:
            line = [key1]
            for key2 in self.secondary_keys:
                val = '0.'
                if key2 in dictionary[key1]:
                    val = str(dictionary[key1][key2])
                line.append(val)
            self.matrix.append(line)
    def write_to_file(self,filename):
        print('Writing a matrix of size',len(self.matrix)-1,'x',len(self.matrix[1])-1)
        f = open(filename,'w')
        for line in self.matrix:
            f.write(','.join(line)+'\n')
        f.close()
    def load_from_file(self, filename):
        self.matrix = []
        f=open(filename,'r')
        for line in f:
            if not ',' in line:
                continue
            self.matrix.append( line.strip('\n').split(',') )
        self.primary_keys   = [line[0] for line in self.matrix[1:]]
        self.secondary_keys = self.matrix[0][1:]

class AnnotationInformationObject:
    def __init__(self, libsbml_document_or_model_annotation, model_name):
        global annotation_count
        if "SBMLDocument" in str(type(libsbml_document_or_model_annotation)):
            mea = semanticSBML.annotate.ModelElementsAnnotations(libsbml_document_or_model_annotation.getModel(), True)
        else:
            mea = libsbml_document_or_model_annotation
        self.model_name = model_name
        self.annotations = {}
        self.description_of_annotations = {}
        for ea in mea.getElementAnnotations():
            element_type = ea.libsbml_element.getElementName()
            if element_type not in 'compartment species reaction'.split():
                continue
            for annotation in ea.getAnnotations():
                abstract_id = database_new_instance.database_instance.__getAbstractIdByIdDict__({annotation.db:annotation.id})
                local_annotation_description = annotation.db+':'+annotation.id
                if abstract_id != None and use_libSBAnnotation:
                    local_annotation_id = str(abstract_id)
                else:
                    local_annotation_id = local_annotation_description
                if use_qsms:
                    qualifier_score = cluster.qualifier_similarity_factor["is"][annotation.qualifier]
                else:
                    qualifier_score = 1.
                if local_annotation_id in self.annotations:
                    annotation_score = self.annotations[local_annotation_id]
                else:
                    annotation_score = 0.
                if add_qsms:
                    annotation_score = annotation_score + qualifier_score
                else:
                    annotation_score = max(annotation_score, qualifier_score)
                self.description_of_annotations[local_annotation_id] = local_annotation_description
                self.annotations[local_annotation_id] = annotation_score
                if use_semantic_density:
                    if local_annotation_id not in annotation_count:
                        annotation_count[local_annotation_id] = 1
                    annotation_count[local_annotation_id] = annotation_count[local_annotation_id] + 1
                        

class AnnotationsInModels:
    def __init__(self, path):
        self.list_of_models = {}
        self.description_of_annotations = {}
        self.model_to_annotations = {}
        filenames = os.listdir(path)
        filenames.sort()
        for filename in filenames:
            if filename.startswith('BIOMD') and filename.endswith('.xml'):
                libsbml_document = libsbml.readSBML(path+os.sep+filename)
                model_name = str(int(filename[5:-4]))+' '+libsbml_document.getModel().getName()
                aio = AnnotationInformationObject(libsbml_document, model_name)
                self.model_to_annotations[aio.model_name] = aio.annotations
                for local_annotation_id in aio.description_of_annotations:
                    if local_annotation_id not in self.description_of_annotations:
                        self.description_of_annotations[local_annotation_id] = aio.description_of_annotations[local_annotation_id]
        _mycsv = mycsv()
        _mycsv.load_from_dictionary(self.model_to_annotations)
        _mycsv.write_to_file(bimatrix_filename)
                            

class CompareAnnotations:
    def __init__(self, max_distance, csv=None):
        global annotation_count
        log_total_annotation_count = scipy.log(sum(annotation_count.values())+1)
        self.max_distance = max_distance
        if csv == None:
            self.csv = mycsv()
            self.csv.load_from_file(bimatrix_filename)
        else:
            self.csv = csv
        self.annotation_keys = self.csv.secondary_keys
        self.known_annotation_keys = [int(x) for x in self.annotation_keys if not ":" in x]
        self.direct_neighbourhood_cache = {}
        self.neighbourhood_cache = {}
        i = 0
        j = len(self.known_annotation_keys)
        libSBAnnotation.stderr.downloadProgress(j,i,True,False,"local neighbourhoods")
        for known_annotation in self.known_annotation_keys:
            self.neighbourhood_cache[known_annotation] = self.generate_neighbourhood(known_annotation)
            i=i+1 
            libSBAnnotation.stderr.downloadProgress(j,i,False,False)
        libSBAnnotation.stderr.downloadProgress(j,i,False,True)
        self.similarity_matrix = []
        self.aid_to_SingleAnnotation = {}
        for index1 in range(len(self.annotation_keys)):
            anno1 = self.annotation_keys[index1]
            if not ':' in anno1:
                anno1 = int(anno1)
            similarity_line = []
            for index2 in range(len(self.annotation_keys)):
                anno2 = self.annotation_keys[index2]
                if not ':' in anno2:
                    anno2 = int(anno2)
                if index1 == index2:
                    this_similarity = 1.
                elif type(anno1) == type(1) and type(anno2) == type(1) and anno2 in self.neighbourhood_cache[anno1]:
                    sa1 = self.abstract_id_to_single_annotation(anno1)
                    sa2 = self.abstract_id_to_single_annotation(anno2)
                    this_similarity = sa1.similarity(sa2)
                else:
                    this_similarity = 0.
                if use_semantic_density:
                    local_semantic_density1 = 1. - (scipy.log(annotation_count[anno1])/log_total_annotation_count) if anno1 in annotation_count else 1. - (1./log_total_annotation_count)
                    local_semantic_density2 = 1. - (scipy.log(annotation_count[anno2])/log_total_annotation_count) if anno2 in annotation_count else 1. - (1./log_total_annotation_count)
                    this_similarity = this_similarity * local_semantic_density1 * local_semantic_density2
                similarity_line.append(this_similarity)
            self.similarity_matrix.append(similarity_line)
        self.similarity_matrix = scipy.sparse.csc_matrix(self.similarity_matrix)
    def abstract_id_to_single_annotation(self, abstract_id):
        if abstract_id not in self.aid_to_SingleAnnotation:
            db = "libSBAnnotation"
            id = abstract_id
            qualifier = 'is'
            self.aid_to_SingleAnnotation[abstract_id] = cluster.SingleAnnotation( fake_annotations.Annotation( db, id, qualifier ) )
            self.aid_to_SingleAnnotation[abstract_id].max_steps = self.max_distance
            self.aid_to_SingleAnnotation[abstract_id].construct_aid_distances()
        return self.aid_to_SingleAnnotation[abstract_id]
    def generate_neighbourhood(self, abstract_id):
        neighbourhood = [abstract_id]
        queue = []
        nextqueue = [abstract_id]
        for step in range(self.max_distance):
            queue = nextqueue
            nextqueue = []
            for aid in queue:
                neighbours = self.generate_direct_neighbourhood(aid)
                for neighbour in neighbours:
                    if neighbour not in neighbourhood:
                        nextqueue.append(neighbour)
                        neighbourhood.append(neighbour)
        return neighbourhood
    def generate_direct_neighbourhood(self, abstract_id):
        if abstract_id in self.direct_neighbourhood_cache:
            return self.direct_neighbourhood_cache[abstract_id]
        relation_object = database_new_instance.database_instance.__getRelationshipsByAbstractId__(abstract_id)
        related_abstract_ids = relation_object.getRelatedAbstractIds()
        self.direct_neighbourhood_cache[abstract_id] = related_abstract_ids
        return related_abstract_ids
    def dump_similarity_matrix(self, filename):
        import pickle
        f = open(filename, 'w')
        pickle.dump(self.similarity_matrix, f)
        f.close()

class FakeCompareAnnotations:
    def __init__(self):
        self.similarity_matrix = []
        self.csv = mycsv()
        self.csv.load_from_file(bimatrix_filename)
        self.annotation_keys = self.csv.secondary_keys
        self.known_annotation_keys = [int(x) for x in self.annotation_keys if not ":" in x]
    def load_similarity_matrix(self, filename):
        import pickle
        f = open(filename, 'r')
        self.similarity_matrix = pickle.load(f)
        f.close()

class EuklideanDistance:
    def __init__(self, compare_annotations):
        self.compare_annotations = compare_annotations
        self.model_to_vector = {}
        for line in self.compare_annotations.csv.matrix[1:]:
            model_name = line[0]
            vector = [float(x) for x in line[1:]]
            vector = scipy.matrix(vector).transpose()
            self.model_to_vector[model_name] = vector
        self.similarity_matrix = self.compare_annotations.similarity_matrix
    def similarity(self, vector1, vector2):
        len_vector1 = scipy.sqrt( (vector1.transpose()*vector1)[0,0] )
        len_vector2 = scipy.sqrt( (vector2.transpose()*vector2)[0,0] )
        if len_vector1 == 0 or len_vector2 == 0:
            similarity = 0.
        elif not use_similarity_in_euklidian:
            similarity = (vector1.transpose()*vector2)[0,0] / len_vector1 / len_vector2
        else:
            # the following lines should be commented in again in order to make the similarity measure the angle between the transformed vectors
            len_vector1 = scipy.sqrt( (vector1.transpose()*self.similarity_matrix*vector1)[0,0] )
            len_vector2 = scipy.sqrt( (vector2.transpose()*self.similarity_matrix*vector2)[0,0] )
            similarity = (vector1.transpose()*self.similarity_matrix*vector2)[0,0] / len_vector1 / len_vector2
        return similarity
    def compute_matrix(self,filename):
        # output in Falkos graph format
        f = open(filename, 'w')
        model_names = self.compare_annotations.csv.primary_keys
        for index1 in range(len(model_names)):
            print(index1)
            for index2 in range(index1+1, len(model_names)):
                sim = self.similarity(self.model_to_vector[model_names[index1]],
                                      self.model_to_vector[model_names[index2]])
                if sim > 0.:
                    bm1 = "BM"+model_names[index1].split()[0]
                    bm2 = "BM"+model_names[index2].split()[0]
                    f.write(bm1+'\t(sim)\t'+bm2+'\t=\t'+str(sim)+'\n')
        f.close()

class DatabaseAnnotation:
    def __init__(self):
        self.max_steps = 1
        self.empty = True
    def isEmpty(self):
        return self.empty
    def compute_matrix(self, filename):
        return self.ed.compute_matrix(filename)
    def load_from_path(self, path):
        self.aim = AnnotationsInModels(path)
        self.ca  = CompareAnnotations(self.max_steps)
        # putting the new neighbours into the bimatrix
        for known_aid in self.ca.neighbourhood_cache:
            neighbour_aids = self.ca.neighbourhood_cache[known_aid]
            for neighbour_aid in neighbour_aids:
                if neighbour_aid not in self.ca.known_annotation_keys:
                    # put this aid into the matrix, a similarity will be computed from it later
                    self.ca.csv.secondary_keys.append(str(neighbour_aid))
                    self.ca.known_annotation_keys.append(neighbour_aid)
                    self.ca.csv.matrix[0].append(str(neighbour_aid))
                    for line in self.ca.csv.matrix[1:]:
                        line.append('0.')
        # saving the new bimatrix
        self.ca.csv.write_to_file(bimatrix_filename)
        # reload the Compare Annotations instance
        self.ca = CompareAnnotations(self.max_steps)
        self.ca.dump_similarity_matrix(simmatrix_filename)
        self.load_from_cache()
    def load_from_cache(self):
        self.fca = FakeCompareAnnotations()
        self.fca.load_similarity_matrix(simmatrix_filename)
        self.ed = EuklideanDistance(self.fca)
        self.empty = False
    def compare_to_file(self, libsbml_document_or_filename, with_overlap=False):
        if type(libsbml_document_or_filename) == type(""):
            libsbml_document = libsbml.readSBML(libsbml_document_or_filename)
            filename = libsbml_document_or_filename
        else:
            libsbml_document = libsbml_document_or_filename
            filename = ""
        aio = AnnotationInformationObject(libsbml_document, filename)
        return self.compare_to_annotation_information_object(aio, with_overlap=with_overlap)
    def compare_model_annotation_object(self, model_annotation_object, with_overlap=False):
        aio = AnnotationInformationObject(model_annotation_object, "new file")
        return self.compare_to_annotation_information_object(aio, with_overlap=with_overlap)
    def compare_to_annotation_information_object(self, aio, with_overlap=False):
        vector = []
        for aid in self.fca.annotation_keys:
            value = 0.
            if aid in aio.annotations:
                value = aio.annotations[aid]
            vector.append(value)
        vector = scipy.matrix(vector).transpose()
        model_names = self.fca.csv.primary_keys
        model_names_to_sim = {}
        for index1 in range(len(model_names)):
            sim = self.ed.similarity(self.ed.model_to_vector[model_names[index1]],
                                     vector)
            overlap = (self.ed.model_to_vector[model_names[index1]].transpose()*vector)[0,0]
            returned_model_name = int(model_names[index1].split()[0])
            if with_overlap:
                model_names_to_sim[returned_model_name] = [sim, overlap]
            else:
                model_names_to_sim[returned_model_name] = sim
        return model_names_to_sim
    def write_dict_dict_file(self, filename):
        import pickle
        f=open(filename,'w')
        pickle.dump(self.biomodels_similarity(),f)
        f.close()
    def biomodels_similarity(self):
        dd = {}
        model_names = self.ed.compare_annotations.csv.primary_keys
        for model_name1 in model_names:
            model_number1 = int(model_name1.split()[0])
            dd[model_number1] = {}
            for model_name2 in model_names:
                model_number2 = int(model_name2.split()[0])
                sim = self.ed.similarity(self.ed.model_to_vector[model_name1],
                                         self.ed.model_to_vector[model_name2])
                dd[model_number1][model_number2] = sim
        return dd


if __name__ == "__main__":
    import datetime
    
    aim = AnnotationsInModels(sys.argv[1])
    ca = CompareAnnotations(1)
    ed = EuklideanDistance(ca)
    #ed.compute_matrix('results/ce'+result_number+'.attribute')
    
    da = DatabaseAnnotation()
    da.load_from_path(sys.argv[1])

    #da.load_from_cache()
    da.compute_matrix('results/ce'+result_number+'.attribute')
    print(da.compare_to_file("/home/schulzma/Desktop/cancer_relevant_receptors.xml"))
    #da.write_dict_dict_file("biomodels_dd.pickle")

