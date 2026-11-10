#!/usr/bin/env python
import numpy
import copy, datetime, math
import libsbml
import scipy
import scipy.sparse
import scipy.sparse.linalg
import scipy.linalg
import semanticSBML.annotate
import libSBAnnotation.database_new
import libSBAnnotation.config
import libSBAnnotation.cofactors
import libModelSimilarity.cluster
import libModelSimilarity.qualifiers

database_new_instance = libSBAnnotation.database_new.Meta_DB()
qualifier_similarity_factor = libModelSimilarity.qualifiers.Qualifiers().load_from_file(libSBAnnotation.config.Config().getpath('resources','qualifiers'))

def annotation2aid(annotation):
    return database_new_instance.database_instance.__getAbstractIdByIdDict__({annotation.db:annotation.id})

def annotation2key(annotation):
    """Converts an annotation into a string representation"""
    abstract_id = annotation2aid(annotation)
    if abstract_id is not None:
        return str(abstract_id)
    else:
        return "%s:%s"%(annotation.db, annotation.id)

def is_int_str(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def logit(x):
    return numpy.log(x / (1.00001 - x))

def logistic(x):
    return 1. / (1. + numpy.exp(-x))

# True reproduces the original program (Python 2 integer division in SigmaMatrix), see CHANGES.md
HEIGHT_INTEGER_DIVISION = False

class SigmaMatrix:
    def __init__(self, SBMLDocuments, model_ids, max_distance_neighbours = 3, use_full_database=False):
        start = datetime.datetime.now()
        self.SBMLDocuments = SBMLDocuments
        self.model_ids = model_ids
        self.max_distance_neighbours = max_distance_neighbours
        meas = []
        for doc in self.SBMLDocuments:
            meas.append( semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True) )
        self.use_full_database = use_full_database
        if self.use_full_database:
            anno_count = len(database_new_instance.database_instance.abstract_item)
            i = 0
            self.keys = []
            while i<anno_count:
                self.keys.append(str(i))
                i=i+1
        else:
            keys = set()
            [[[keys.add(annotation2key(annotation)) for annotation in ea.getAnnotations()] for ea in mea.getElementAnnotations()] for mea in meas]
            self.keys = list(keys)
            self.generate_neighbours()
        self.keys.sort()
        self.compute_heights()
        self.generate_matrix()
    def _compute_single_height(self, aid, max_height=30, blacklist=[]):
        if max_height <=0:
            return 1
        if aid in blacklist:
            return 30 # max_height
        parents = list(database_new_instance.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, "is_a", 1, "to"))
        if parents == []:
            return 1
        else:
            new_blacklist = [aid]+blacklist
            min_height = min([self._compute_single_height(p, max_height-1, new_blacklist) for p in parents])
            return min_height + 1
    def _parents(self, aid):
        if aid not in self._parent_cache:
            self._parent_cache[aid] = list(database_new_instance.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, "is_a", 1, "to"))
        return self._parent_cache[aid]
    def _height(self, aid, max_height=30):
        """
        Height of an entry in its ontology: 1 + number of is_a steps on the shortest way down
        to an entry that has no more specific entries below it.
        Gives the same values as L{_compute_single_height} (the original, recursive
        function, which visits every path and needs minutes for large models), but
        searches breadth first and remembers the parents of each entry.
        """
        seen = set([aid])
        level = [aid]
        steps = 0
        while level and steps < max_height:
            next_level = []
            for a in level:
                parents = self._parents(a)
                if not parents:
                    return steps + 1
                for p in parents:
                    if p not in seen:
                        seen.add(p)
                        next_level.append(p)
            level = next_level
            steps += 1
        return max_height + 1
    def compute_heights(self):
        self.heights = {}
        self._parent_cache = {}
        for key in self.keys:
            if is_int_str(key):
                self.heights[key] = self._height(int(key))
            else:
                self.heights[key] = 1
    def expand_single_aid(self, expand):
        relation_object = database_new_instance.database_instance.__getRelationshipsByAbstractId__(expand)
        return relation_object.getRelatedAbstractIds()
    def generate_neighbours(self):
        all_aids = [int(aid) for aid in self.keys if is_int_str(aid)]
        queue = copy.deepcopy(all_aids)
        for i in range(self.max_distance_neighbours):
            new_queue = []
            for expand in queue:
                related_abstract_ids = self.expand_single_aid(expand)
                for aid in related_abstract_ids:
                    if aid not in all_aids:
                        all_aids.append(aid)
                        new_queue.append(aid)
            queue = new_queue
        for aid in all_aids:
            said = str(aid)
            if said not in self.keys:
                self.keys.append(said)
    def _expand_relation(self, start_aid, current_aid, sim, depth):
        if depth == 0:
            return
        aids2value = self.aid2related[current_aid]
        start_index = self.keys.index(str(start_aid))
        for next_aid in aids2value:
            value = sim * aids2value[next_aid]
            next_index = self.keys.index(str(next_aid))
            if self.sigma[start_index, next_index] < value:
                self.sigma[start_index, next_index] = value
                self.sigma[next_index, start_index] = value
                self._expand_relation(start_aid, next_aid, value, depth-1)
    def generate_matrix(self):
        n = len(self.keys)
        self.sigma = scipy.sparse.lil_matrix((n, n))
        for i in range(n):
            self.sigma[i, i] = 1.
        self.aid2related = {}
        for (index, aid) in [(index, int(aid)) for (index, aid) in enumerate(self.keys) if is_int_str(aid)]:
            self.aid2related[aid] = {}
            relation_object = database_new_instance.database_instance.__getRelationshipsByAbstractId__(aid)
            for [other_aid, foo, relation_type, bar] in relation_object.relations:
                if str(other_aid) in self.keys:
                    value = libModelSimilarity.cluster.relation_type_score[relation_type]
                    self.aid2related[aid][other_aid] = value
        for aid in self.aid2related:
            self._expand_relation(aid, aid, 1., self.max_distance_neighbours)
        for i in range(len(self.sigma.rows)):
            for k in range(len(self.sigma.rows[i])):
                j = self.sigma.rows[i][k]
                height_sum = self.heights[self.keys[i]] + self.heights[self.keys[j]]
                height = 2 // height_sum if HEIGHT_INTEGER_DIVISION else 2 / height_sum
                self.sigma[i, j] = self.sigma[i, j] ** height
        self.sigma = self.sigma.tocsr()
    def _element_annotation_to_vector(self, ea):
        vec = numpy.zeros(len(self.keys))
        for a in ea.getAnnotations():
            str_anno = annotation2key(a)
            index_anno = self.keys.index(str_anno)
            if index_anno == -1:
                raise Exception('Annotation not in Sigma matrix! Maybe SigmaMatrix has not been initialized with the correct models.')
            qsm = qualifier_similarity_factor["is"][a.qualifier]
            vec[index_anno] = qsm
        return numpy.matrix([vec])
    def _sigma_multiply(self, v1, v2):
        return (v1 * self.sigma * v2.transpose())[0,0]
    def compute_similarity_sbml_elements(self, element1, element2):
        eas = [semanticSBML.annotate.ElementAnnotations(e) for e in [element1, element2]]
        return self.compute_similarity_eas(*eas)
    def compute_similarity_eas(self, ea1, ea2):
        if ea1.type != ea2.type:
            return 0.
        [vec1, vec2] = [self._element_annotation_to_vector(ea) for ea in [ea1, ea2]]
        return self._sigma_multiply(vec1, vec2) / numpy.sqrt(self._sigma_multiply(vec1, vec1)) / numpy.sqrt(self._sigma_multiply(vec2, vec2))
    def compute_pairwise_similarities(self):
        self.esm = ElementSimilarityMatrix()
        n = len(self.SBMLDocuments)
        for i in range(n):
            doc1 = self.SBMLDocuments[i]
            model1 = self.model_ids[i]
            mea1 = semanticSBML.annotate.ModelElementsAnnotations(doc1.getModel(), True)
            for j in range(i+1, n):
                doc2 = self.SBMLDocuments[j]
                model2 = self.model_ids[j]
                mea2 = semanticSBML.annotate.ModelElementsAnnotations(doc2.getModel(), True)
                for ea1 in mea1.getElementAnnotations():
                    for ea2 in mea2.getElementAnnotations():
                        sim = self.compute_similarity_eas(ea1, ea2)
                        if sim > 0.:
                            self.esm.add_similarity(model1, model2, ea1.id, ea2.id, sim)
        return self.esm
            
        
class ElementSimilarityMatrix:
    def __init__(self):
        self.matrix = {}
        self.tuples = []
        self.tuples_copied = None
    def rebuild_stack(self):
        self.tuples_copied = sorted(self.tuples, key=lambda x:x[0], reverse=True)
    def pop_similarity(self):
        """Returns a list [sim, model1, model2, id1, id2] containing the element with the biggest similarity."""
        if self.tuples_copied is None:
            self.rebuild_stack()
        if len(self.tuples_copied) == 0:
            return None
        return self.tuples_copied.pop(0)
    def add_similarity(self, model1, model2, id1, id2, sim):
        self._add_similarity([model1, model2, id1, id2], self.matrix, sim)
        self._add_similarity([model2, model1, id2, id1], self.matrix, sim)
        self.tuples.append([sim, model1, model2, id1, id2])
    def _add_similarity(self, keys, _dict, sim):
        if len(keys) == 1:
            _dict[keys[0]] = sim
        else:
            if keys[0] not in _dict:
                _dict[keys[0]] = {}
            self._add_similarity(keys[1:], _dict[keys[0]], sim)
    def get_similarity(self, model1, model2, id1, id2):
        try:
            return self.matrix[model1][model2][id1][id2]
        except KeyError:
            return 0.
    def _get_model_and_id(self, sbml_element):
        ea = semanticSBML.annotate.ElementAnnotations(sbml_element)
        model = sbml_element.getModel().getId()
        return [model, ea.id]
    #def get_similarity_for_sbml_elements(self, element1, element2):
    #    [[model1, id1], [model2, id2]] = [self._get_model_and_id(e) for e in [element1, element2]]
    #    return self.get_similarity(model1, model2, id1, id2)
    def __str__(self):
        s = ""
        if len(self.matrix) == 2:
            #comparing two models
            mee2s = self.matrix[list(self.matrix.keys())[0]]
            ee2s = mee2s[list(mee2s.keys())[0]]
            for elem1 in ee2s:
                line = []
                for elem2 in ee2s[elem1]:
                    line.append(["", elem2, str(ee2s[elem1][elem2])])
                line[0][0] = elem1
                line = "\n".join(["\t".join(l) for l in line])
                s += line + "\n"
        return s

class GlobalSigma:
    dependencies = {"Compartment":[[["Outside"], "Compartments", "Id", "=="],
                                   [["CompartmentType"], "CompartmentTypes", "Id", "=="]],
                    "Species":[[["Compartment"], "Compartments", "Id", "=="],
                               [["SpeciesType"], "SpeciesTypes", "Id", "=="]],
                    "Reaction":[[["ListOfReactants","Species"], "Species", "Id", "in"],
                                [["ListOfModifiers","Species"], "Species", "Id", "in"],
                                [["ListOfProducts" ,"Species"], "Species", "Id", "in"]],
                    "Event":[],
                    "Constraint":[],
                    "InitialAssignment":[],
                    "AlgebraicRule":[],
                    "AssignmentRule":[],
                    "RateRule":[],
                    }
    half_dependencies = {"Compartment":[],
                    "Species":[[["SpeciesType"], "SpeciesTypes", "Id", "=="]],
                    "Reaction":[[["ListOfReactants","Species"], "Species", "Id", "in"],
                                #[["ListOfModifiers","Species"], "Species", "Id", "in"],
                                [["ListOfProducts" ,"Species"], "Species", "Id", "in"]],
                    "Event":[],
                    "Constraint":[],
                    "InitialAssignment":[],
                    "AlgebraicRule":[],
                    "AssignmentRule":[],
                    "RateRule":[],
                    }
    def __init__(self, SBMLDocuments, model_ids, meas, sigma_matrix, limited_propagation=None):
        self.SBMLDocuments = SBMLDocuments
        self.model_ids = model_ids
        self.meas = meas
        self.sigma_matrix = sigma_matrix
        self.limited_propagation = limited_propagation
        self.cofactors = libSBAnnotation.cofactors.Cofactors()
        self.local_esm = self.sigma_matrix.compute_pairwise_similarities()
        self.global_esm = ElementSimilarityMatrix()
        self.get_global_esm()
    def _get_type(self, element):
        try:
            a = semanticSBML.annotate.ElementAnnotations(element, True)
            return a.type.capitalize()
        except:
            return None
    def _get_model_elements(self, document):
        elements = []
        model = document.getModel()
        list_functions = [x for x in dir(model) if x.startswith("getListOf")]
        for list_function in list_functions:
            #for element in eval("model.%s()"%list_function):
            #    if "getId" in dir(element) and self._get_type(element) is not None:
            for element in getattr(model, list_function)():
                if hasattr(element, "getId") and self._get_type(element) is not None:
                    elements.append(element)
        return elements
    def _make_global_vectors(self, document):
        elements = self._get_model_elements(document)
        element_ids = [x.getId() for x in elements]
        R = numpy.zeros((len(elements), len(elements)))
        r = 0.5
        vecs = []
        for i in range(len(elements)):
            element = elements[i]
            i_is_cofactor = self.cofactors.is_cofactor_libsbml_element(element)
            for other_element in self.get_dependent_elements(document, element):
                j_is_cofactor = self.cofactors.is_cofactor_libsbml_element(other_element)
                j = elements.index(other_element)
                if not j_is_cofactor and not i_is_cofactor:
                    R[i, j] = r
                    R[j, i] = r
            ea = semanticSBML.annotate.ElementAnnotations(element, True)
            vecs.append( self.sigma_matrix._element_annotation_to_vector(ea).tolist()[0] )
        I = numpy.matrix(numpy.eye(len(elements)))
        largest_eigenvalue = float(max(scipy.linalg.eigvals(R).real))
        if largest_eigenvalue != 0.:
            _lambda = 0.5 / largest_eigenvalue
        else:
            _lambda = 0.5
        #_lambda = 0.1
        v = numpy.matrix(vecs)
        if self.limited_propagation is None:
            pi = (I - _lambda * R).getI()
        elif self.limited_propagation == 0:
            pi = I
        elif self.limited_propagation == 1:
            pi = I + _lambda * R
        elif self.limited_propagation == 2:
            pi = I + _lambda * R + (_lambda**2 * R * R)
        p = pi * v
        global_vecs = []
        for global_vec in p.tolist():
            global_vecs.append( numpy.matrix([global_vec]) )
        return global_vecs
    def get_global_esm(self):
        self.global_vector = [self._make_global_vectors(doc) for doc in self.SBMLDocuments]
        self.global_esm = ElementSimilarityMatrix()
        n = len(self.SBMLDocuments)
        for i in range(n):
            for j in range(i+1, n):
                for k,e1 in enumerate(self._get_model_elements(self.SBMLDocuments[i])):
                    v1 = self.global_vector[i][k]
                    for l,e2 in enumerate(self._get_model_elements(self.SBMLDocuments[j])):
                        v2 = self.global_vector[j][l]
                        if type(e1) == type(e2):
                            sim = self.sigma_matrix._sigma_multiply(v1, v2) / numpy.sqrt(self.sigma_matrix._sigma_multiply(v1, v1)) / numpy.sqrt(self.sigma_matrix._sigma_multiply(v2, v2))
                            if sim > 0.:
                                self.global_esm.add_similarity(self.model_ids[i], self.model_ids[j], e1.getId(), e2.getId(), sim)
        return self.global_esm
    def get_dependent_elements(self, doc, element):
        _type = self._get_type(element)
        deps = []
        if _type not in self.dependencies:
            return deps
        for [command, _list, attr, operator] in self.dependencies[_type]:
            try:
                command_available = eval("element.isSet%s()"%command[0])
                if not command_available:
                    continue
            except AttributeError:
                # the isSet method is not defined because the command attribute if obligatory
                pass
            val = eval("element.get%s()"%command[0])
            if len(command) > 1:
                for com in command[1:]:
                    val = [eval("x.get%s()"%com) for x in val]
            for other in eval("doc.getModel().getListOf%s()"%_list):
                if eval("other.get%s() %s val"%(attr, operator)):
                    deps.append(other)
        return deps

class GlobalSigmaSimProp(GlobalSigma):
    def __init__(self, SBMLDocuments, model_ids, meas, sigma_matrix):
        GlobalSigma.__init__(self, SBMLDocuments, model_ids, meas, sigma_matrix)
    def get_global_esm(self):
        import scipy.sparse
        import scipy.sparse.linalg
        import copy
        self.global_esm = ElementSimilarityMatrix()
        n = len(self.SBMLDocuments)
        _lambda = 0.2
        r = 0.5
        for i in range(n):
            ilen = len(self._get_model_elements(self.SBMLDocuments[i]))
            model1 = self.model_ids[i]
            for j in range(i+1, n):
                jlen = len(self._get_model_elements(self.SBMLDocuments[j]))
                model2 = self.model_ids[j]
                # compute similarities for each model pair
                #R = scipy.sparse.lil_matrix((ilen*jlen, ilen*jlen))
                R = numpy.zeros((ilen*jlen, ilen*jlen))
                sigma = numpy.matrix([[0.]]*(ilen*jlen))
                for k,e1 in enumerate(self._get_model_elements(self.SBMLDocuments[i])):
                    kdep = [semanticSBML.annotate.ElementAnnotations(x).id  for x in self.get_dependent_elements(self.SBMLDocuments[i], e1)]
                    for l,e2 in enumerate(self._get_model_elements(self.SBMLDocuments[j])):
                        # compute entries for the \delta - R matrix
                        ldep = [semanticSBML.annotate.ElementAnnotations(x).id for x in self.get_dependent_elements(self.SBMLDocuments[j], e2)]
                        for o,e3 in enumerate(self._get_model_elements(self.SBMLDocuments[i])):
                            oid = semanticSBML.annotate.ElementAnnotations(e3).id
                            if oid not in kdep:
                                continue
                            for p,e4 in enumerate(self._get_model_elements(self.SBMLDocuments[j])):
                                pid = semanticSBML.annotate.ElementAnnotations(e4).id
                                if pid not in ldep:
                                    continue
                                R[(k*jlen+l), (o*jlen+p)] = r*r
                                R[(o*jlen+p), (k*jlen+l)] = r*r
                        # compute entries for the sigma
                        local_sim = self.local_esm.get_similarity(model1, model2, semanticSBML.annotate.ElementAnnotations(e1).id, semanticSBML.annotate.ElementAnnotations(e2).id)
                        if local_sim != 0.:
                            sigma[(k*jlen+l), 0] = local_sim
                # transform R for faster computations
                #R = R.tocsc()
                R = scipy.sparse.csc_matrix(R)
                # rescale R
                #largest_eigenvalue = scipy.sparse.linalg.eigen(R, k=1, which='LR').tolist()[0]
                largest_eigenvalue = float(max(scipy.linalg.eigvals(R.toarray()).real))
                # R = I - _lambda*R
                R = (-0.5 / largest_eigenvalue) * R
                R.setdiag(numpy.ones(ilen*jlen))
                # solve rho = R^{-1} * sigma
                (rho, info) = scipy.sparse.linalg.cgs(R, sigma, x0=copy.deepcopy(sigma))
                if info > 0:
                    (rho, info) = scipy.sparse.linalg.bicg(R, sigma, x0=copy.deepcopy(sigma))
                    if info > 0:
                        # does not work and will crash
                        rho = R.getI() * sigma
                # set entries in the global_esm
                rho = numpy.matrix(rho).transpose()
                i_elements = self._get_model_elements(self.SBMLDocuments[i])
                j_elements = self._get_model_elements(self.SBMLDocuments[j])
                # find row and column maxima
                ks_max = [0.] * len(i_elements)
                ls_max = [0.] * len(j_elements)
                for k in range(len(i_elements)):
                    for l in range(len(j_elements)):
                        rkl = rho[(k*jlen+l), 0]
                        for (_list, li) in [(ks_max, k), (ls_max, l)]:
                            if rkl > _list[li]:
                                _list[li] = rkl
                # normalize by them
                #for k in range(len(i_elements)):
                #    for l in range(len(j_elements)):
                #        for (_list, li) in [(ks_max, k), (ls_max, l)]:
                #            if _list[li] > 0.0000000001:
                #                rho[(k*jlen+l), 0] = rho[(k*jlen+l), 0] / math.sqrt(_list[li])
                # put them into the global element similarity matrix
                for k,e1 in enumerate(i_elements):
                    for l,e2 in enumerate(j_elements):
                        if rho[(k*jlen+l), 0] > 0.0000000001 or rho[(k*jlen+l), 0] < -0.0000000001:
                            #print semanticSBML.annotate.ElementAnnotations(e1).id, semanticSBML.annotate.ElementAnnotations(e2).id, rho[(k*jlen+l), 0]
                            self.global_esm.add_similarity(model1, model2, semanticSBML.annotate.ElementAnnotations(e1).id, semanticSBML.annotate.ElementAnnotations(e2).id, rho[(k*jlen+l), 0])
                
class Match:
    def __init__(self, models, ids, esm_instance):
        self.models = models
        self.ids = ids
        self.esm = esm_instance
    def to_Falko(self):
        l = []
        for i in range(len(self.models)):
            l.append({"mid":self.models[i],
                      "element_id":self.ids[i]})
        return l
    def overlap(self, other):
        for i in range(len(self.models)):
            for j in range(len(other.models)):
                if self.models[i] == other.models[j] and self.ids[i] == other.ids[j]:
                    return True
        return False
    def exclude(self, other):
        for i in range(len(self.models)):
            for j in range(len(other.models)):
                if self.models[i] == other.models[j] and self.ids[i] != other.ids[j]:
                    return True
        return False
    def logit_score(self):
        logit_score = 0.
        if len(self.models) > 1:
            for i in range(len(self.models)):
                for j in range(i+1, len(self.models)):
                    this_sim = self.esm.get_similarity(self.models[i], self.models[j], self.ids[i], self.ids[j])
                    logit_score += logit(this_sim)
        return logit_score
    def score(self):
        offset = -3. # random choice for general similarity
        return logistic(self.logit_score() + offset)
    def merge(self, other):
        new = copy.deepcopy(self)
        for i in range(len(other.models)):
            if [True for j in range(len(new.models)) if new.models[j] == other.models[i] and new.ids[j] == other.ids[j]].count(True) == 0:
                # model id combination is not in new object
                new.models.append(other.models[i])
                new.ids.append(other.ids[i])
        return new
            
        
    
                    
class ModelMatcher:
    def __init__(self, SBMLDocuments, model_ids, max_distance_neighbours=3, use_feature_propagation=True):
        self.SBMLDocuments = SBMLDocuments
        self.model_ids = model_ids
        self.meas = [semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True) for doc in self.SBMLDocuments]
        self.sigma = SigmaMatrix(self.SBMLDocuments, model_ids, max_distance_neighbours)
        self.global_sigma_anno_prop = None
        self.global_sigma_sim_prop = None
        if use_feature_propagation:
            self.global_sigma_anno_prop = GlobalSigma(self.SBMLDocuments, self.model_ids, self.meas, self.sigma)
            self.global_sigma = self.global_sigma_anno_prop
        else:
            self.global_sigma_sim_prop = GlobalSigmaSimProp(self.SBMLDocuments, self.model_ids, self.meas, self.sigma)
            self.global_sigma = self.global_sigma_sim_prop
        self.global_esm = self.global_sigma.global_esm
    def pickle(self):
        for i in range(len(self.SBMLDocuments)):
            self.SBMLDocuments[i] = self.SBMLDocuments[i].toSBML()
        del self.meas
        for o in [self.sigma, self.global_sigma_anno_prop, self.global_sigma_sim_prop]:
            if o is not None:
                try:
                    del o.meas
                except:
                    pass
    def unpickle(self):
        for i in range(len(self.SBMLDocuments)):
            self.SBMLDocuments[i] = libsbml.readSBMLFromString("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"+self.SBMLDocuments[i])
        self.meas = [semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), True) for doc in self.SBMLDocuments]
        self.global_sigma.meas = self.meas
    def get_similarity(self, models, ids):
        match = Match(models, ids, self.global_esm)
        return match.score()
    def pick_propagation_type(self, propagation_type=None):
        if propagation_type is None:
            return
        elif propagation_type == 'local':
            self.global_sigma = self.sigma
            self.global_esm = self.global_sigma.esm
            self.global_sigma.meas = self.meas
        elif propagation_type == 'annotation_propagation':
            if self.global_sigma_anno_prop is None:
                self.global_sigma_anno_prop = GlobalSigma(self.SBMLDocuments, self.model_ids, self.meas, self.sigma)
            self.global_sigma = self.global_sigma_anno_prop
            self.global_esm = self.global_sigma.global_esm
        elif propagation_type == 'similarity_propagation':
            if self.global_sigma_sim_prop is None:
                self.global_sigma_sim_prop = GlobalSigmaSimProp(self.SBMLDocuments, self.model_ids, self.meas, self.sigma)
            self.global_sigma = self.global_sigma_sim_prop
            self.global_esm = self.global_sigma.global_esm
    def match(self, minimum_similarity=0.):
        matches = [x.to_Falko() for x in self._match(minimum_similarity)]
        # setting the element type for all the elements in the matches
        for match in matches:
            for d in match:
                mid = d["mid"]
                eid = d["element_id"]
                mea = self.meas[self.model_ids.index(mid)]
                _type = [ea.type for ea in mea.getElementAnnotations() if ea.id == eid][0]
                d["element_type"] = _type
        self.global_esm.rebuild_stack()
        return matches
    def _match(self, minimum_similarity=0.):
        self.all_matches = []
        while True:
            pop = self.global_esm.pop_similarity()
            if pop is None:
                break
            [sim, model1, model2, id1, id2] = pop
            if sim < minimum_similarity:
                break
            new_match = Match([model1, model2], [id1, id2], self.global_esm)
            conflicts = [x for x in self.all_matches if new_match.overlap(x)]
            if conflicts == []:
                # no conflicting matches
                self.all_matches.append(new_match)
            else:
                # look whether the conflict is resolvable
                exclusions = [x for x in conflicts if new_match.exclude(x)]
                if exclusions == []:
                    # we can handle the conflict
                    merged = copy.deepcopy(new_match)
                    for conflict in conflicts:
                        merged = merged.merge(conflict)
                    # look for use of two ids from the same model in the merged ids
                    if len(list(set(merged.models))) < len(merged.models):
                        continue
                    # look if the merged contributes to the overall similarity
                    single_logit_scores = sum([x.logit_score() for x in conflicts])
                    new_logit_score = merged.logit_score()
                    if new_logit_score <= single_logit_scores:
                        continue
                    # otherwise take the new, merged Match
                    self.all_matches = [merged] + [x for x in self.all_matches if not new_match.overlap(x)]
        return self.all_matches
            
            
        
        
            
def main2():
    docs = [libsbml.readSBML(f) for f in ['/home/schulzma/Desktop/shortchain.xml']*2]
    mm = ModelMatcher(docs,[1, 2], 1, False)
    for match in mm.match(.0001):
        print(match)
    
def main():
    import os, datetime
    path = 'libModelSimilarity/curated/'
    files = [path+x for x in os.listdir(path) if x.endswith('.xml')]
    files.sort()
    sbml_docs = [libsbml.readSBML(f) for f in files]
    docs = [sbml_docs[60],sbml_docs[63]]
    mm = ModelMatcher(docs,[1, 2], 1, True)
    for match in mm.match(.8):
        print(match)
    print(mm.get_similarity([1, 2], ['AMP', 'AMP']))

    import pickle
    mm.pickle()
    ms = pickle.dumps(mm)
    mm = pickle.loads(ms)
    mm.unpickle()

if __name__ == "__main__":
    main()
