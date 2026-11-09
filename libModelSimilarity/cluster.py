#!/usr/bin/env python

import libsbml, semanticSBML.annotate, os, libSBAnnotation.database_new, libSBAnnotation.database, datetime, sys, math, scipy, libSBAnnotation.stderr
from . import html, qualifiers

similarity_measure_version = 6
path = './'

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="file to compare to the rest of the models", metavar="FILE", default=None)
    parser.add_option("-p", "--path", dest="path", help="path to the BioModels repository", metavar="PATH", default=None)
    parser.add_option("-s", "--smv", dest="similarity_measure_version", default="6", help="selects the similarity measure", metavar="1-6")
    parser.add_option("-d", "--dump", dest="dump", default=False, action="store_true", help="writes the cached annotated models to disk for faster comparison to single files")
    parser.add_option("-t", "--statistics", dest="show_statistics", default=False, action="store_true", help="prints the statistics of the dump")
    (options, args) = parser.parse_args()
    if options.path == None:
        print("You need to set a path!")
        parser.print_help()
        sys.exit(1)
    similarity_measure_version = int(options.similarity_measure_version)
    path = options.path


starting_time = datetime.datetime.now()
database_new_instance = libSBAnnotation.database_new.Meta_DB()
miriam2internal = dict([[libSBAnnotation.database.internal_to_miriam_names[x],x] for x in libSBAnnotation.database.internal_to_miriam_names])
qualifier_similarity_factor = qualifiers.Qualifiers().load_from_file()
relation_type_score = qualifiers.relation_type_score
height_cache = {}
all_tops = set([])
hash_count_dict = {}
pseudo_hash_count = 1
apply_caching = True

class SingleAnnotation:
    def __init__(self, annotation):
        self.db = str(annotation.db)
        self.id = str(annotation.id)
        self.qualifier = str(annotation.qualifier)
        self.aid_distance_cache = {}
        self.max_steps = 1
        [self.height, self.top] = [0, None]
        self.hash = self.db+":"+self.id
        if self.db in miriam2internal:
            self.db = miriam2internal[self.db]
            self.abstract_id = database_new_instance.database_instance.__getAbstractIdByIdDict__({self.db:self.id})
        elif self.db == 'libSBAnnotation':
            self.abstract_id = int(self.id)
        else:
            self.abstract_id = None
        if self.abstract_id != None:
            if similarity_measure_version >= 4:
                [self.height, self.top] = self.determine_height(self.abstract_id)
                all_tops.add(self.top)
            self.hash = self.abstract_id
            self.construct_aid_distances()
        global pseudo_hash_count
        if self.hash not in hash_count_dict:
            hash_count_dict[self.hash] = 1
            pseudo_hash_count = pseudo_hash_count + 1
        hash_count_dict[self.hash] = hash_count_dict[self.hash] + 1
        pseudo_hash_count = pseudo_hash_count + 1
    def __cmp__(self, other):
        if self.abstract_id != None and other.abstract_id != None:
            return self.abstract_id == other.abstract_id
        else:
            return self.db == other.db and self.id == other.id and self.qualifier == other.qualifier
    def __repr__(self):
        return self.qualifier+": "+self.db+": "+self.id+" identified as "+str(self.abstract_id)
    def information_factor(self):
        # if similarity_measure_version >=6:
        #     max_hc = max(hash_count_dict.values())
        #     return 1. - (math.log(hash_count_dict[self.hash])/math.log(max_hc))**3.
        if similarity_measure_version >=5:
            return 1. - (math.log(hash_count_dict[self.hash])/math.log(pseudo_hash_count))
        else:
            return 1.
    def determine_height(self, aid, max_depth=20):
        if aid in height_cache:
            return height_cache[aid]
        parents = list(database_new_instance.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, "is_a", 1, "to"))
        if len(parents) == 0 or max_depth == 0:
            return_val = [0, aid]
        else:
            max_dist = -1
            top = None
            for parent_aid in parents:
                [this_dist,this_top] = self.determine_height(parent_aid, max_depth-1)
                if this_dist > max_dist:
                    max_dist = this_dist
                    top = this_top
            return_val = [max_dist+1, top]
        height_cache[aid] = return_val
        return return_val
    def construct_aid_distances(self):
        other_aids = database_new_instance.database_instance.__getRelatedAbstractIdsByAbstractId__(self.abstract_id)
        for other_aid in other_aids:
            self.__judge_abstract_id_distance__(self.abstract_id, other_aid)
    def judge_abstract_id_distance(self, aid1, aid2):
        aid_hash = str(aid1)+" "+str(aid2)
        if aid_hash in self.aid_distance_cache:
            return self.aid_distance_cache[aid_hash]
        else:
            return [0., []]
    def __judge_abstract_id_distance__(self, aid1, aid2, max_steps=None):
        """Function to judge the distance between two abstract ids"""
        if max_steps == None:
            max_steps = self.max_steps
        aid_hash = str(aid1)+" "+str(aid2)
        if aid_hash in self.aid_distance_cache:
            return self.aid_distance_cache[aid_hash]
        if max_steps == 0:
            if aid1 == aid2:
                return [1., []]
            else:
                return [0., []]
        else:
            maximum_sim = 0.
            max_trace = None
            relation_object = database_new_instance.database_instance.__getRelationshipsByAbstractId__(aid1)
            for [new_aid, d2id, relation_type, reverse] in relation_object.relations:
                [child_score, trace] = self.__judge_abstract_id_distance__(new_aid, aid2, max_steps-1)
                this_score = relation_type_score[relation_type] * child_score
                if this_score > maximum_sim:
                    maximum_sim = this_score
                    max_trace = [[aid1, new_aid, relation_type, reverse]] + trace
            self.aid_distance_cache[aid_hash] = [maximum_sim, max_trace]
            return [maximum_sim, max_trace]

    def similarity(self,other):
        information_factor = self.information_factor() * other.information_factor()
        height_correction_factor = 1.
        tangens_height_correction_factor = 1.
        if similarity_measure_version >=4 and self.height != 0 and other.height != 0:
            height_correction_factor = 2./float(self.height+other.height)
            tangens_height_correction_factor = scipy.tanh(height_correction_factor*.1)
        if similarity_measure_version >=3 and self.abstract_id != None and other.abstract_id != None:
            # compare the two elements based on their abstract_id
            qsfs = qualifier_similarity_factor[self.qualifier][other.qualifier]
            if self.abstract_id == other.abstract_id:
                return qsfs
            [aid_similarity, trace] = self.judge_abstract_id_distance(self.abstract_id, other.abstract_id)
            return aid_similarity**height_correction_factor * qsfs * information_factor
            #return aid_similarity * tangens_height_correction_factor * qsfs * information_factor
        elif similarity_measure_version >= 2:
            if self.db == other.db and self.id == other.id:
                # use the qualifiers to compute the similarity
                return qualifier_similarity_factor[self.qualifier][other.qualifier]
            else:
                return 0.
        elif similarity_measure_version >= 1:
            # similarity without regarding the qualifiers
            if self.db == other.db and self.id == other.id and self.qualifier == other.qualifier:
                return 1.
            else:
                return 0.
            
            

class ElementAnnotation:
    def __init__(self,annotation_list):
        self.libsbml_element_id = annotation_list.libsbml_element.getId()
        self.libsbml_element_name = annotation_list.libsbml_element.getName()
        self.libsbml_element_type = str(type(annotation_list.libsbml_element))
        self.single_annotation_list = []
        for anno in annotation_list.getAnnotations():
            self.single_annotation_list.append(SingleAnnotation(anno))
        self.abstract_id_list = list(set( [sa.abstract_id for sa in self.single_annotation_list if sa.abstract_id != None and sa.qualifier == "is"] ))
    def __cmp__(self,other):
        if self.abstract_id != None and other.abstract_id != None:
            return self.abstract_id == other.abstract_id
        else:
            for sa in self.single_annotation_list:
                if sa not in other.single_annotation_list:
                    return False
            for sa in other.single_annotation_list:
                if sa not in self.single_annotation_list:
                    return False
            return True
    def __repr__(self):
        return str(self.libsbml_element_id)+" ("+str(self.libsbml_element_name)+") "+", ".join([repr(sa) for sa in self.single_annotation_list])
    def similarity(self,other):
        # first of all check, whether the libsbml elements which are described
        # are of the same type. If not, then the similarity should be 0
        if self.libsbml_element_type != other.libsbml_element_type and not "instance" in self.libsbml_element_type and not "instance" in other.libsbml_element_type:
            return 0.
        # check if annotations are available in both elements
        relation_divisor = float(len(self.single_annotation_list) + len(other.single_annotation_list))
        if relation_divisor == 0.:
            return 0.
        # go through the annotations
        total_sum = 0.
        for [sum_list, max_list] in [[self.single_annotation_list, other.single_annotation_list],
                                     [other.single_annotation_list, self.single_annotation_list]]:
            for first_elem in sum_list:
                maximum = 0.
                for second_elem in max_list:
                    dist = first_elem.similarity(second_elem)
                    if dist > maximum:
                        maximum = dist
                total_sum = total_sum + maximum
        return total_sum / relation_divisor

class ModelAnnotation:
    def __init__(self, model_elements_annotations, model_name):
        self.model_name = model_name
        self.element_annotations = []
        i = 0
        j = len(model_elements_annotations.getElementAnnotations())
        libSBAnnotation.stderr.downloadProgress(j,i,True,False,"model "+model_name+" into cache")
        for anel in model_elements_annotations.getElementAnnotations():
            element_type = str(type(anel.libsbml_element)).split("\'")[1].split('.')[1]
            if element_type in 'Model Parameter Event RateRule AssignmentRule'.split():
                continue
            self.element_annotations.append(ElementAnnotation(anel))
            i=i+1 
            libSBAnnotation.stderr.downloadProgress(j,i,False,False)
        libSBAnnotation.stderr.downloadProgress(j,i,False,True)
        self.hash = str(self)
        self.similarity_cache = {}
    def similarity(self, other, absolute=False):
        max_sim = max(self._similarity(self, absolute), other._similarity(other, absolute))
        if max_sim != 0.:
            return self._similarity(other, absolute) / max_sim
        else:
            return 0.
    def _similarity(self,other, absolute=False):
        if apply_caching and other.hash in self.similarity_cache:
            return self.similarity_cache[other.hash]
        relation_divisor = float(len(self.element_annotations) + len(other.element_annotations))
        if relation_divisor == 0.:
            sim = 0.
        else:
            if absolute:
                relation_divisor = 1.
            total_sum = 0.
            for [sum_list, max_list] in [[self.element_annotations, other.element_annotations],
                                         [other.element_annotations, self.element_annotations]]:
                for first_elem in sum_list:
                    maximum = 0.
                    for second_elem in max_list:
                        dist = first_elem.similarity(second_elem)
                        if dist > maximum:
                            maximum = dist
                    total_sum = total_sum + maximum
            sim = total_sum / relation_divisor
        # store the similarity on both sides
        self.similarity_cache[other.hash] = sim
        other.similarity_cache[self.hash] = sim
        return sim

class Cluster:
    def __init__(self, elements, subclusters):
        self.elements    = elements
        self.subclusters = subclusters
        self.similarity_cache = {}
    def similarity(self,other):
        other_hash = str(other)
        if apply_caching and other_hash in self.similarity_cache:
            return self.similarity_cache[other_hash]
        distance_list = []
        for e1 in self.elements:
            for e2 in other.elements:
                distance_list.append(e1.similarity(e2))
        dist = min(distance_list) # single_linkage
        #dist = max(distance_list) # complete_linkage
        #dist = sum(distance_list)/float(len(distance_list)) # average_linkage
        self.similarity_cache[other_hash] = dist
        return dist

        

class DatabaseAnnotation(object):
    _shared_state = {}
    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls._shared_state
        if 'model_annotations' not in self.__dict__:
            self.model_annotations = []
        return self
    def isEmpty(self):
        return len(self.model_annotations) == 0
    def cluster(self):
        # create single clusters
        cluster_list = [Cluster([ma],[]) for ma in self.model_annotations]
        while len(cluster_list) > 1:
            print(len(cluster_list),"clusters left")
            # determine the distances
            distances = {}
            for index1 in range(len(cluster_list)):
                cluster1 = cluster_list[index1]
                for index2 in range(index1+1, len(cluster_list)):
                    cluster2 = cluster_list[index2]
                    distances[cluster1.similarity(cluster2)] = [index1, index2]
                    print("\t",(index1,index2))
            # choose clusters to be merged
            minimal_distance = max(distances.keys())
            [index1, index2] = distances[minimal_distance]
            cluster1 = cluster_list[index1]
            cluster2 = cluster_list[index2]
            # create new cluster
            new_cluster = Cluster(cluster1.elements+cluster2.elements, [cluster1, cluster2])
            # delete old from the list
            cluster_list.remove(cluster1)
            cluster_list.remove(cluster2)
            cluster_list.append(new_cluster)
        return cluster_list[0]
    def matrix(self,sort=True):
        """Creates a matrix from the distances and displays it as an html file."""
        if sort:
            elements = self.cluster().elements
        else:
            elements = self.model_annotations
        dim = len(elements)+1
        # create the matrix
        out_matrix = []
        for i in range(dim):
            out_matrix.append([""]*dim)
        # fill the border
        for index in range(len(elements)):
            out_matrix[index][0] = elements[index].model_name
            out_matrix[-1][index+1] = elements[index].model_name
        # fill in the elements
        for index1 in range(len(elements)):
            for index2 in range(len(elements)):
                out_matrix[index1][index2+1] = elements[index1].similarity(elements[index2])
        # creating the html
        h=html.html()
        h.set_matrix(out_matrix)
        h.tofile("cluster"+str(similarity_measure_version)+".html")
    def compare_to_file(self, filename, absolute=False):
        doc = libsbml.readSBML(filename)
        model_annotation = ModelAnnotation( semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True), "New Model" )
        return self.compare_model_annotation_object(model_annotation, absolute)
    def compare_model_annotation_object(self, model_annotation, absolute=False):
        model2sim = {}
        i = 0
        j = len(self.model_annotations)
        libSBAnnotation.stderr.downloadProgress(j, i, True, False, "model similarities")
        for ma in self.model_annotations:
            if absolute:
                model2sim[ma.model_name] = model_annotation._similarity( ma, absolute=True )
            else:
                model2sim[ma.model_name] = model_annotation.similarity( ma, absolute=False )
            i = i + 1
            libSBAnnotation.stderr.downloadProgress(j, i, False, False)
        libSBAnnotation.stderr.downloadProgress(j, i, False, True)
        return model2sim
        #values = list(set(model2sim.values()))
        #values.sort()
        #values.reverse()
        #for val in values:
        #    for name in [name for name,val2 in model2sim.iteritems() if val == val2]:
        #        print val,":",name
    def load_from_path(self, path):
        if not path.endswith(os.sep):
            path = path+os.sep
        modelfiles = [path+f for f in os.listdir(path) if f.startswith('BIOMD') and f.endswith('.xml')]
        modelfiles.sort()
        #self.doc_list = []
        # testing models
        glycolysis = [42, 61, 64, 70, 71, 172, 176, 177, 206, 211]
        glyc_sides = [13, 17, 23, 51]
        circ_clock = [16, 21, 22, 24, 25, 34, 36, 55, 73, 74, 78]
        cell_cycle = [5, 6, 7, 8, 56, 109, 110, 111, 144, 150, 168, 169, 181, 195, 196, 207, 208]
        folate = [18, 213]
        mapk = [9, 10, 11, 14, 19, 26, 27, 28, 29, 30, 31, 49, 146]
        calcium = [39,43,44,45,57,58,59,98,100,113,114,115,122,123,124,128,138,145,161,162,166,184]
        four_groups = [70,71,211,16,21,22,5,7,8,111,26,27,28,29]
        for filename in modelfiles:
            modelnumber = str(int(filename.split(os.sep)[-1][5:].split('.')[0]))
            #if int(modelnumber) not in four_groups:
            #    continue
            doc = libsbml.readSBML(filename)
            #self.doc_list.append(doc)
            self.model_annotations.append( ModelAnnotation( semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True), "BM"+modelnumber ) )
        # updating the height cache for ways up
        if similarity_measure_version >=6:
            global done_aid
            done_aid = []
            for aid in range(len(database_new_instance.database_instance.abstract_item.rows)):
                parents = list(database_new_instance.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, "is_a", 1, "from"))
                if len(parents)==0:
                    self._update_step_down(aid,0)
    def _update_step_down(self, aid, step):
        global hash_count_dict, done_aid
        this_count = 0
        if aid in hash_count_dict:
            this_count = hash_count_dict[aid]
            if aid in done_aid:
                return this_count
        children = list(database_new_instance.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, "is_a", 1, "to"))
        for child in children:
            child_count = self._update_step_down(child,step+1)
            this_count = this_count + child_count
        hash_count_dict[aid] = this_count
        done_aid.append(aid)
        return this_count
    def countAnnotations(self):
        count_dict = {}
        for ma in self.model_annotations:
            for ea in ma.element_annotations:
                if len(ea.single_annotation_list) == 0:
                    continue
                length = len(ea.abstract_id_list)
                if length not in count_dict:
                    count_dict[length] = 0
                count_dict[length] = count_dict[length]+1
        keys = list(count_dict.keys())
        keys.sort()
        for key in keys:
            print(key,':',count_dict[key])
    def annotationStatistics(self):
        count_dict = {}
        for ma in self.model_annotations:
            for ea in ma.element_annotations:
                for sa in ea.single_annotation_list:
                    if sa.abstract_id != None:
                        if sa.abstract_id not in count_dict:
                            count_dict[sa.abstract_id] = 0
                        count_dict[sa.abstract_id] = count_dict[sa.abstract_id] + 1
        count_dict2 = {}
        for aid,num in count_dict.items():
            if num not in count_dict2:
                count_dict2[num] = []
            count_dict2[num].append(aid)
        keys = list(count_dict2.keys())
        keys.sort()
        for key in keys[:40]:
            aids = count_dict2[key]
            print(key,"is the number of models", end=' ')
            print(", ".join([database_new_instance.database_instance.__getNameByAbstractId__(aid) for aid in aids]), end=' ')
            print("appear in")
    def falkos_files(self):
        sim_file = file("results/c"+str(similarity_measure_version)+".sim","w")
        att_file = file("results/c"+str(similarity_measure_version)+".attribute","w")
        for i in range(len(self.model_annotations)):
            model_i = self.model_annotations[i]
            for j in range(i+1, len(self.model_annotations)):
                model_j = self.model_annotations[j]
                this_sim = model_i.similarity(model_j)
                if this_sim > 0.0:
                    sim_file.write(model_i.model_name+"\tsim\t"+model_j.model_name+"\n")
                    att_file.write(model_i.model_name+"\t(sim)\t"+model_j.model_name+"\t=\t"+str(this_sim)+"\n")
        sim_file.close()
        att_file.close()
    def falkos_statistics(self):
        import popen2
        sim_list = []
        num_bins = 10
        popen2.popen2("octave --persist --eval \"hist("+str(sim_list)+","+num_bins+")\"")
        
        

if __name__ == "__main__":
    if options.filename == None:
        if not options.show_statistics:
            # cluster all biomodels
            da=DatabaseAnnotation()
            da.load_from_path(path)
            da.falkos_files()
            if options.dump:
                print("Dumping the biomodels to disk!")
                print("This may take a while.")
                f=open("dump","w")
                import pickle
                pickle.dump([da, hash_count_dict, pseudo_hash_count],f)
                f.close()
        else:
            try:
                f=open("dump","r")
                print("Reading dumped biomodels")
                import pickle
                [da, hash_count_dict, pseudo_hash_count] = pickle.load(f)
                f.close()
                da.annotationStatistics()
            except IOError:
                # file not present
                print("dump not present")
    else:
        # compare one model to the rest of the biomodels
        # try to load a dumped biomodels cache
        try:
            f=open("dump","r")
            print("Reading dumped biomodels")
            import pickle
            [da, hash_count_dict, pseudo_hash_count] = pickle.load(f)
            f.close()
        except IOError:
            # file not present
            print("Reconstructing the database on the fly")
            da=DatabaseAnnotation()
            da.load_from_path(path)
        # doing the comparison
        for filename in options.filename.split(','):
            print("Comparing file",filename, "to the database")
            da.compare_to_file(filename)
        


