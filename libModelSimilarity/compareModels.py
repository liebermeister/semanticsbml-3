import semanticSBML.annotate, libSBAnnotation, scipy
from . import calc_bimatrix

class ComparableModel:
    def __init__(self, compare_models_instance, filename = None, sbmldoc = None, model_elements_annotations = None):
        self.compare_models_instance = compare_models_instance
        self.libsbml_element_ids = []
        self.libsbml_element_names = []
        self.sbmldoc = None
        if model_elements_annotations != None:
            self.mea = model_elements_annotations
            if self.mea.libsbml_document != None:
                self.sbmldoc = self.mea.libsbml_document
        elif filename != None or sbmldoc != None:
            if filename != None:
                self.sbmldoc = self._filename_to_libsbml_document_instance(filename)
            else:
                self.sbmldoc = sbmldoc
            self.mea = semanticSBML.annotate.ModelElementsAnnotations(self.sbmldoc.getModel(), True)
        else:
            raise Hell()
        self.matrix = self._model_elements_annotation_to_element_matrix(self.mea)
    def _filename_to_libsbml_document_instance(self, filename):
        if type(filename) != type("string"):
            # assuming we already have a libsbml document
            return filename
        import libsbml
        return libsbml.readSBML(filename)
    def _model_elements_annotation_to_element_vector(self, mea):
        matrix = self._model_elements_annotation_to_element_matrix(mea)
        vector = [0.]* len(self.compare_models_instance.annotation_keys)
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > .5:
                    vector[j] = 1.
        return vector
    def _model_elements_annotation_to_element_matrix(self, mea):
        matrix = []
        self.description_of_annotations = {}
        for ea in mea.getElementAnnotations():
            annotation_strings = []
            element_type = ea.type#libsbml_element.getElementName()
            if element_type not in 'compartment species reaction'.split():
                continue
            self.libsbml_element_ids.append(ea.id)#libsbml_element.getId())
            self.libsbml_element_names.append(ea.name)#libsbml_element.getName())
            for annotation in ea.getAnnotations():
                abstract_id = self.compare_models_instance.database_new_instance.database_instance.__getAbstractIdByIdDict__({annotation.db:annotation.id})
                local_annotation_description = annotation.db+':'+annotation.id
                if abstract_id != None:
                    local_annotation_id = str(abstract_id)
                else:
                    local_annotation_id = local_annotation_description
                self.description_of_annotations[local_annotation_id] = local_annotation_description
                if local_annotation_id not in annotation_strings:
                    annotation_strings.append(local_annotation_id)
            vector = []
            for s in self.compare_models_instance.annotation_keys:
                if s in annotation_strings:
                    vector.append(1.)
                else:
                    vector.append(0.)
            matrix.append(vector)
        return matrix
    def __to_annotation_element_json_graph__(self, postfix = ''):
        json = []
        for element_annotation in self.mea.getElementAnnotations():
            this_json_dict = {}
            this_json_dict["id"] = element_annotation.id+postfix
            this_json_dict["name"] = element_annotation.name
            for annotation in element_annotation.getAnnotations():
                suggestion = annotation.getName() + " (" + this_json_dict["id"] + ")"
                if suggestion != None and suggestion != "" and annotation.qualifier == "is":
                    this_json_dict["name"] = suggestion
                    break
            this_json_dict["data"] = {"$dim":10.,
                                      "$type":"circle",
                                      "$color":"#888888"
                                      }
            this_json_dict["adjacencies"] = []
            json.append(this_json_dict)
        return json
    def to_json_graph(self,postfix = ''):
        if self.sbmldoc == None:
            return self.__to_annotation_element_json_graph__(postfix)
        m = self.sbmldoc.getModel()
        json = []
        for species in m.getListOfSpecies():
            this_json_dict = {}
            this_json_dict["id"] = species.getId()+postfix
            this_json_dict["name"] = species.getId()
            if species.isSetName():
                this_json_dict["name"] = species.getName() + " (" + this_json_dict["id"] + ")"
            element_annotation = semanticSBML.annotate.ElementAnnotations( species )
            for annotation in element_annotation.getAnnotations():
                suggestion = annotation.getName() + " (" + this_json_dict["id"] + ")"
                if suggestion != None and suggestion != "" and annotation.qualifier == "is":
                    this_json_dict["name"] = suggestion
                    break
            this_json_dict["data"] = {"$dim":10.,
                                      "$type":"circle",
                                      "$color":"#888888"
                                      }
            this_json_dict["adjacencies"] = []
            for reaction in m.getListOfReactions():
                reactants = [x.getSpecies() for x in reaction.getListOfReactants()]
                products = [x.getSpecies() for x in reaction.getListOfProducts()]
                modifiers = [x.getSpecies() for x in reaction.getListOfModifiers()]
                if species.getId() in reactants+products+modifiers:
                    if species.getId() in reactants:
                        data = {"$color":"#00ff00",
                                "$direction":[species.getId()+postfix,reaction.getId()+postfix],
                                "$type":"arrow",
                                "role":"reactant",
                                "distance":20,
                                "draw":1
                                }
                    elif species.getId() in products:
                        data = {"$color":"#ff0000",
                                "$direction":[reaction.getId()+postfix,species.getId()+postfix],
                                "$type":"arrow",
                                "role":"product",
                                "distance":20,
                                "draw":1
                                }
                    elif species.getId():
                        data = {"$color":"#0000ff",
                                "role":"modifier",
                                "distance":100,
                                "draw":1
                                }
                    data["weight"] = 1
                    this_json_dict["adjacencies"].append( {"nodeTo":reaction.getId()+postfix,
                                                           "data":data
                                                           } )
            json.append(this_json_dict)
        for reaction in m.getListOfReactions():
            this_json_dict = {}
            this_json_dict["id"] = reaction.getId()+postfix
            this_json_dict["name"] = reaction.getId()#species.getName()
            this_json_dict["data"] = {"$dim":10.,
                                      "$type":"square",
                                      "$color":"#888888"
                                      }
            this_json_dict["adjacencies"] = []
            participants = [x.getSpecies() for x in reaction.getListOfReactants()] + [x.getSpecies() for x in reaction.getListOfProducts()]
            for participant in participants:
                this_json_dict["adjacencies"].append( {"nodeTo":participant+postfix,
                                                       "data":{"weight":1,
                                                               "draw":0}
                                                       } )
            json.append(this_json_dict)
        return json


class CompareModels:
    def __init__(self):
        self.database_new_instance = libSBAnnotation.database_new.Meta_DB()
        self.database_annotation = calc_bimatrix.DatabaseAnnotation()
        self.database_annotation.load_from_cache()
        self.csv = calc_bimatrix.mycsv()
        self.csv.load_from_file(calc_bimatrix.bimatrix_filename)
        self.annotation_keys = self.csv.secondary_keys
        fca = calc_bimatrix.FakeCompareAnnotations()
        fca.load_similarity_matrix(calc_bimatrix.simmatrix_filename)
        self.similarity_matrix = fca.similarity_matrix
        self.current_element_matrix = None
    def similarity(self, vector1, vector2):
        len_vector1 = scipy.sqrt( (vector1.transpose()*vector1)[0,0] )
        len_vector2 = scipy.sqrt( (vector2.transpose()*vector2)[0,0] )
        if len_vector1 == 0 or len_vector2 == 0:
            similarity = 0.
        else:
            len_vector1 = scipy.sqrt( (vector1.transpose()*self.similarity_matrix*vector1)[0,0] )
            len_vector2 = scipy.sqrt( (vector2.transpose()*self.similarity_matrix*vector2)[0,0] )
            similarity = (vector1.transpose()*self.similarity_matrix*vector2)[0,0] / len_vector1 / len_vector2
        return similarity

    def compare_models(self, model1, model2):
        self.comparable_model1 = ComparableModel(self, model_elements_annotations=model1)
        self.comparable_model2 = ComparableModel(self, model_elements_annotations=model2)
        matrix = []
        for i in range(len(self.comparable_model1.matrix)):
            line = []
            vector1 = scipy.matrix(self.comparable_model1.matrix[i]).transpose()
            for j in range(len(self.comparable_model2.matrix)):
                vector2 = scipy.matrix(self.comparable_model2.matrix[j]).transpose()
                this_sim = self.similarity(vector1, vector2)
                line.append(this_sim)
            matrix.append(line)
        self.current_element_matrix = matrix
        return matrix

    def matrix_to_greedy_sim_list(self,threshold=.25):
        best_sims = []
        deadi = []
        deadj = []
        for n in range(min([len(self.current_element_matrix),len(self.current_element_matrix[0])])):
            # compute the biggest entry in the matrix
            max_val = -1
            max_ij  = [-1,-1]
            for i in range(len(self.current_element_matrix)):
                if i in deadi:
                    continue
                for j in range(len(self.current_element_matrix[i])):
                    if j in deadj:
                        continue
                    if self.current_element_matrix[i][j] > max_val:
                        max_val = self.current_element_matrix[i][j]
                        max_ij = [i,j]
            if max_val < threshold:
                break
            best_sims.append( [self.comparable_model1.libsbml_element_ids[max_ij[0]],
                               self.comparable_model2.libsbml_element_ids[max_ij[1]],
                               max_val] )
            deadi.append(max_ij[0])
            deadj.append(max_ij[1])
        return best_sims

    def matrix_to_sim_list(self,threshold=.25):
        sims = []
        for i in range(len(self.current_element_matrix)):
            for j in range(len(self.current_element_matrix[i])):
                
                sims.append( [self.comparable_model1.libsbml_element_ids[i],
                              self.comparable_model2.libsbml_element_ids[j],
                              self.comparable_model1.libsbml_element_names[i],
                              self.comparable_model2.libsbml_element_names[j],
                              self.current_element_matrix[i][j]] )
        return sims

    def element_to_max_sim(self):
        sim1 = {}
        sim2 = {}
        for i in range(len(self.current_element_matrix)):
            sim1[self.comparable_model1.libsbml_element_ids[i]] = max(self.current_element_matrix[i])
        for j in range(len(self.current_element_matrix[0])):
            sim2[self.comparable_model2.libsbml_element_ids[j]] = max([x[j] for x in self.current_element_matrix])
        return [sim1,sim2]

    def __split_nodes__(self,nodes, max_degree=3):
        import copy
        split_indices = []
        # finding the nodes which have to be split
        for node,index in zip(nodes,list(range(len(nodes)))):
            if len(node['adjacencies'])>max_degree and node['data']['$type'] == 'circle':
                split_indices.append(index)
        # removing them from the list
        split_indices.reverse()
        split_nodes = []
        for index in split_indices:
            split_nodes.append(nodes[index])
            del nodes[index]
        # split them and put them into the other list
        for node in split_nodes:
            for adjacency in node['adjacencies']:
                new_node = copy.deepcopy(node)
                new_id = node['id']+"_adjacent_split_to_"+adjacency['nodeTo']
                new_node['id'] = new_id
                new_node['adjacencies'] = [adjacency]
                nodeTo = adjacency['nodeTo']
                for other_node in nodes+split_nodes:
                    if other_node['id'] == nodeTo:
                        for other_adjacency in other_node['adjacencies']:
                            if other_adjacency['nodeTo'] == node['id']:
                                other_adjacency['nodeTo'] = new_id
                nodes.append(new_node)
                                                                                                                                                                                                            
    def __sim2hex__(self,sim,rgb='rgb'):
        hex_num = str(hex(int((1.-sim)*100+154)))[2:].rjust(2,"0")
        return_string = '#'
        for letter in 'rgb':
            if letter in rgb:
                return_string = return_string + hex_num
            else:
                return_string = return_string + '00'
        return return_string
                                            
    def to_global_graph(self,threshold=.25,max_degree=3):
        json1 = self.comparable_model1.to_json_graph("_1")
        json2 = self.comparable_model2.to_json_graph("_2")
        if json1 == [] or json2 == []:
            # no sbml document available
            return []
        sims = self.matrix_to_greedy_sim_list(threshold)
        [element_to_sim1, element_to_sim2] = self.element_to_max_sim()
        for [id1,id2,sim] in sims:
            try:
                node1 = [j for j in json1 if j['id']==id1+"_1"][0]
                node2 = [j for j in json2 if j['id']==id2+"_2"][0]
            except IndexError:
                # In the sims list compartments might be present which will not appear in the json
                # objects. Therefore, we will just neglect these similarities
                continue
            # check for the number of adjacencies in both nodes if they exceed 3, then dont put a link there
            if len(node1['adjacencies'])>max_degree or len(node2['adjacencies'])>max_degree:
                continue
            node1['adjacencies'].append({'nodeTo':id2+"_2",'data':{'$color': '#ff7f00',
                                                                   'distance': 1,
                                                                   'draw': 1,
                                                                   'weight': 1}})
            node2['adjacencies'].append({'nodeTo':id1+"_1",'data':{'$color': '#ff7f00',
                                                                   'distance': 1,
                                                                   'draw': 1,
                                                                   'weight': 1}})
        graph = json1+json2
        for [e2s,letter,postfix] in [[element_to_sim1,'r','_1'],[element_to_sim2,'b','_2']]:
            for node in graph:
                this_id = node['id']
                if this_id.endswith(postfix):
                    node['data']['$color'] = self.__sim2hex__(e2s[this_id[:-(len(postfix))]],letter)
        self.__split_nodes__(graph,max_degree=max_degree)
        return graph

def main():
    cm = CompareModels()
    cm.compare_models('curated/BIOMD0000000070.xml','curated/BIOMD0000000071.xml')
    print(cm.matrix_to_sim_list())
    print(cm.matrix_to_greedy_sim_list())
    print(cm.to_global_graph())
    print(cm.element_to_max_sim())

def bugtest():
    import libsbml, os.path
    cm = CompareModels()
    doc1 = libsbml.readSBML(os.path.expanduser(os.path.join('~','Desktop','KEGG_threonine.xml')))
    model1 = semanticSBML.annotate.ModelElementsAnnotations(None, libsbml_document=doc1, suppress_errors='')
    doc2 = libsbml.readSBML(os.path.expanduser(os.path.join('~','semanticsbml','trunk','libModelSimilarity','curated','BIOMD0000000066.xml')))
    model2 = semanticSBML.annotate.ModelElementsAnnotations(None, libsbml_document=doc2, suppress_errors='')
    cm.compare_models(model1, model2)
    cm.matrix_to_greedy_sim_list()
    cm.to_global_graph(max_degree=9)

if __name__ == "__main__":
    main()
    #bugtest()
