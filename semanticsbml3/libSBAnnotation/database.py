import copy
from . import indexedTables, config, stderr, file
from functools import cmp_to_key

table_list = "abstract_item string_repr abstract_item_relation takes_part localisation".split()
ext_db_names = [
    ["KEGG Compound","kegg_compound"],["KEGG Reaction","kegg_reaction"],
    ["KEGG Formula","kegg_formula"],["Formula","formula"],   # formula: from ChEBI (added in the port)
    ["ChEBI","chebi"],["3DMET","3dmet"],["CAS","cas"],
    ["PDB-CCD","pdb_ccd"],["PubChem","pubchem_substance"],
    ["PubChem-substance","pubchem_substance"],
    ["PubChem-compound","pubchem_compound"],
    ["Knapsack","knapsack"],["LipidBank","lipidbank"],
    ["LipidMaps","lipidmaps"],["Gmelin","gmelin"],["MSDchem","msdchem"],
    ["Beilstein","beilstein"],["UM-BBD","um_bbd"],["MolBase","molbase"],
    ["DrugBank","drugbank"],["Protein Data Bank","pdb"],
    ["Enzyme Nomenclature","ec"],["GeneOntology","go"],
    ["Gene Ontology","go"],["MetaCyc","metacyc"],
    ["Reactome","reactome"],["AraCyc","aracyc"],["BIND","bind"],
    ["Ensembl","ensembl"],["UniProt","uniprot"],["Taxonomy","taxonomy"],
    ["BioModels Database","biomodels"],["ICD","icd"],
    ["IntAct","intact"],["InterPro","interpro"],
    ["KEGG Pathway","kegg_pathway"],["PubMed","pubmed"],["OMIM","omim"],
    ["PIRSF","pirsf"],["DOI","doi"],["SGD","sgd"],
    ["Systems Biology Ontology","sbo"],["KEGG Drug","kegg_drug"],
    ["LigandBox","ligandbox"], ["KEGG Organism","kegg_organism"],
    ["HMDB","hmdb"],["Nikkaji","nikkaji"],["Entrez Gene","entrez_gene"],
    ["RefSeq","refseq"],["MIM","mim"],["HPRD","hprd"],["HGNC","hgnc"],
    ["IMGT/GENE-DB","imgt_gene_db"],["RGD","rgd"],["TAIR","tair"],
    ["MaizeGDB","maizegdb"],["GeneDB","genedb"],["Wormbase","wormbase"],
    ["ZFIN","zfin"],["Xenbase","xenbase"],["MGI","mgi"],
    ["RatMap","ratmap"],["FlyBase","flybase"],
    ["Nucleotide Sequence Database","nsd"],["Pfam","pfam"],
    ["International Protein Index","ipi"],["PROSITE","prosite"],
    ["BioCyc","biocyc"],["VectorBase","vectorbase"],["CGNC","cgnc"],
    ["PBR","pbr"],["VBRC","vbrc"],["Pathema","pathema"],
    ["Pseudocap","pseudocap"],["ApiDB-CryptoDB","apidb_cryptodb"],
    ["EcoCyc","ecocyc"],["EcoGene","ecogene"],["NIF-Subcellular","nif_subcellular"],
    ["CORUM","corum"],["PMID","pmid"],["TCDB","tcdb"],["DIP","dip"],
    ["NCBI protein GI","ncbi_protein_gi"],["GermOnline","germonline"],
    ["BioGRID","biogrid"],["NCBI Gene ID","ncbi_gene_id"],["GenBank","genbank"],
    ["locus","locus"],["STITCH","stitch"],["BindingDB","bindingdb"],
    ["allyson","allyson"],["druggable","druggable"]
   ]
internal_to_miriam_names = {"kegg_compound":"KEGG Compound",
                            "kegg_reaction":"KEGG Reaction",
                            "kegg_drug":"KEGG Drug",
                            "kegg_enzyme":"KEGG Enzyme",
                            "chebi":"ChEBI",
                            "3dmet":"3DMET",
                            "cas":"CAS",
                            "taxonomy":"Taxonomy",
                            "pubchem_substance":"PubChem-substance",
                            "go":"Gene Ontology",
                            "reactome":"Reactome",
                            "uniprot":"UniProt",
                            "ec":"Enzyme Nomenclature",
                            "entrez_gene":"Entrez Gene",
                            "refseq":"RefSeq",
                            "nsd":"Nucleotide Sequence Database",
                            "interpro":"InterPro",
                            "icd":"ICD",
                            "omim":"OMIM",
                            "intact":"IntAct",
                            "bind":"BIND",
                            "ensemble":"Ensembl",
                            "pirsf":"PIRSF",
                            "sgd":"SGD",
                            "druggable":"druggable",
                            "":"",
                            "":"",
                            "kegg_organism":"KEGG Organism"}

class AnnotationObject:
    """Fake annotation object for current programming purposes."""
    def __init__(self, dictionary):
        self.dictionary = dictionary
    def toDatabase2id(self):
        return self.dictionary
    def toInternalCompliantDatabase2id(self):
        mcdict = {}
        for key,val in self.dictionary.items():
            if key in list(internal_to_miriam_names.values()):
                new_key = [x for x in internal_to_miriam_names if internal_to_miriam_names[x] == key][0]
            else:
                new_key = key
            mcdict[new_key] = val
        return mcdict
    def toMiriamCompliantDatabase2id(self):
        mcdict = {}
        for key,val in self.dictionary.items():
            if key in internal_to_miriam_names:
                mcdict[internal_to_miriam_names[key]] = val
            #else:
            #    print "KEY",key,"unknown to internal_to_miriam_names"
        return mcdict
    def __repr__(self):
        return "Annotation Object "+repr(self.dictionary)
    def iteritems(self):
        return iter(self.dictionary.items())
    def items(self):
        return self.dictionary.items()
    
class RelationObject:
    """Fake relation object for current programming purposes."""
    def __init__(self, abstract_id, database2id=None):
        self.abstract_id = abstract_id
        self.database2id = database2id
        self.relations = []
    def addRelation(self,related_abstract_id,related_database2id,
                    relation_type,reverse=False):
        self.relations.append( [related_abstract_id,related_database2id,
                                relation_type,reverse] )
    def getRelatedAbstractIds(self):
        return [x[0] for x in self.relations]
    def __repr__(self):
        string  = "RelationObject:\n"
        string += "\n".join( [",".join([str(x) for x in row])
                              for row in self.relations] )
        return string
    

class InconsistentDatabaseError(Exception):
    def __init__(self,value):
        Exception.__init__(self,value)
        self.value = value
    def __str__(self):
        return repr(self.value)

class Database(object):
    def __init__(self):
        tablesModule = indexedTables
        self.registered_tables = []
        self.abstract_item = tablesModule.AbstractItem()
        self.registered_tables.append("abstract_item")
        self.string_repr = tablesModule.StringRepr()
        self.registered_tables.append("string_repr")
        self.string_repr_from_database = tablesModule.StringReprFromDatabase()
        self.registered_tables.append("string_repr_from_database")
        self.abstract_item_relation = tablesModule.AbstractItemRelation()
        self.registered_tables.append("abstract_item_relation")
        self.takes_part = tablesModule.TakesPart()
        self.registered_tables.append("takes_part")
        self.localisation = tablesModule.Localisation()
        self.registered_tables.append("localisation")
        # automatically create external item tables
        for display_name, internal_name in ext_db_names:
            table_name = "external_item_"+internal_name
            self.__dict__[table_name] = tablesModule.ExternalItem(display_name)
            if table_name not in self.registered_tables:
                self.registered_tables.append(table_name)
        self.config = config.Config()
        self.filename = self.config.getpath("database","database_dir")+"dump"
        self.generalized_items = []
        self.list_of_imported_databases = []
        self.list_of_internal_names_for_external_resources = [x[1] for x in ext_db_names]
        self.list_of_all_names = None
        self.abstract_id_to_annotation_object = {}
        self._relation_cache = None

        try:
            # get the information on which databases are already loaded before the database content is loaded
            import pickle
            f=open(self.config.getpath("database","database_dir")+"internal_db.main","rb")
            [self.list_of_imported_databases] = pickle.load(f)
            f.close()
        except IOError:
            pass

    @property
    def relation_cache(self):
        """abstract id -> relation type -> direction -> {abstract id: None}; built on first use"""
        if self._relation_cache is None:
            cache = {}
            for row in self.abstract_item_relation.rows:
                for (from_aid,to_aid,direction) in [(row[0],row[1],'from'),(row[1],row[0],'to')]:
                    if from_aid not in cache:
                        cache[from_aid] = {}
                    rcf = cache[from_aid]
                    if row[2] not in rcf:
                        rcf[row[2]] = {}
                    rcfr = rcf[row[2]]
                    if direction not in rcfr:
                        rcfr[direction] = {}
                    rcfr[direction][to_aid] = None
            self._relation_cache = cache
        return self._relation_cache

    def getListOfImportedDatabases(self):
        return [internal_to_miriam_names[x] if x in internal_to_miriam_names else x for x in self.list_of_imported_databases]
    def addToListOfImportedDatabases(self,database):
        self.list_of_imported_databases.append(database)
        
    def resetGeneralizedItems(self):
        self.generalized_items = []
        
    def getListOfAllRelationTypes(self):
        return list(self.abstract_item_relation.global_index["relation_type"].keys())
        
    def dump(self):
        for table_name in self.registered_tables:
            self.__dict__[table_name].dump()
        import pickle
        f=open(self.config.getpath("database","database_dir")+"internal_db.main","wb")
        pickle.dump([self.list_of_imported_databases],f)
        f.close()
        
    def load(self):
        number_of_tables = len(self.registered_tables)
        import os
        if not os.path.exists(self.config.getpath("database","database_dir")+"internal_db.main"):
            return None
        stderr.downloadProgress(number_of_tables,0,True,False,"internal database")
        index = 0
        try:
            for table_name in self.registered_tables:
                self.__dict__[table_name].load()
                index = index + 1
                stderr.downloadProgress(number_of_tables,index,False,False)
        except IOError:
            stderr.downloadProgress(number_of_tables,number_of_tables,False,True)
            return None
        stderr.downloadProgress(number_of_tables,number_of_tables,False,True)

        # testing a cache for abstract item relations
        #self.air_cache = {}
        #for row in self.abstract_item_relation.rows:
        #    for (from_aid,to_aid) in [(row[0],row[1]),(row[1],row[0])]:
        #        if from_aid not in self.air_cache:
        #            self.air_cache[from_aid] = []
        #        self.air_cache[from_aid].append(to_aid)
        # end of explicit caching og abstract item relations

        self._relation_cache = None     # built when first needed, see relation_cache
        # end of second cache for relations
        return self
    
    def __str__(self):
        for variable in self.__dict__:
            if variable.startswith("external_item_"):
                table_list.append(variable)
        return "\n\n".join(str[self.__dict__[x]] for x in table_list)

    def addReaction(self,database2id,substrates,products,modifiers,enzymes, source_database):
        """Adds a reaction. Substrates and products is a list of lists
        with the entries [database2id,name_list,stoich]."""
        # adding the abstract reaction
        abstract_reaction_id = len(self.abstract_item)
        self.abstract_item.insertRow(["reaction"])
        # adding the external database reaction
        for database_name,database_id in database2id.items():
            external_db_table = self.__dict__["external_item_"+database_name]
            external_db_table.insertRow([database_id,
                                         abstract_reaction_id,
                                         False])
        # adding substrates, products, modifiers, enzymes
        for role,list_of_participants in [["substrate", substrates],
                                          ["product", products],
                                          ["modifier", modifiers],
                                          ["enzyme", enzymes]]:
            for participant in list_of_participants:
                database2id = participant[0]
                name_list = participant[1]
                stoich = participant[2] if len(participant)==3 else 1
                # get the abstract id of the participant
                abstract_participant_id = self.__getAbstractIdByIdDict__(database2id)
                if abstract_participant_id == None and len(name_list)>0:
                    abstract_participant_id = self.__getAbstractIdByNameList__(name_list)
                if abstract_participant_id == None:
                    # create the species
                    abstract_participant_id = self.addSpecies(database2id,name_list, source_database)
                    #abstract_participant_id = self.__getAbstractIdByIdDict__(database2id)
                # make an entry in the takes part table
                self.takes_part.insertRow([abstract_participant_id,
                                           abstract_reaction_id, role,
                                           stoich])


    def addSpecies(self, database2id, name_list, source_database):
        return self.addItem(database2id ,name_list, "species", source_database)
    
    def addOrganism(self, database2id, name_list, source_database):
        return self.addItem(database2id ,name_list, "organism", source_database)
    
    def addSBOTerm(self, database2id, name_list, source_database):
        return self.addItem(database2id ,name_list, "sbo_term", source_database)
    
    def addGene(self, database2id, name_list, source_database):
        return self.addItem(database2id ,name_list, "gene", source_database)
    
    def addItem(self, database2id, name_list, item_type, source_database):
        """Given a dictionary of databases to IDs in those databases
        and a list of names in which the first name in the favourite
        description."""

        # problems with tabs in names
        for i in range(len(name_list)):
            if '\t' in name_list[i]:
                name_list[i] = name_list[i].replace('\t',' ')

        # problems with the automated resolution of entries mapped to multiple entries in a different database
        deldb=[]
        for db in database2id:
            if type(database2id[db]) == type([]):
                deldb.append(db)
        for db in deldb:
            del database2id[db]

        more_generalized_item_links = []
        more_generalized_abstract_items = []
        more_specialized_abstract_items = []
        # take care of the fact that the values in database2id can be lists
        remove_databases_from_database2id = []
        for database_name,list_of_ids in database2id.items():
            if type(list_of_ids) == type([]):
                remove_databases_from_database2id.append(database_name)
                for database_id in list_of_ids:
                    # create abstract items
                    more_specialized_abstract_items.append(
                        self.addItem({database_name:database_id},[],item_type, source_database)
                        )
        for database_name in remove_databases_from_database2id:
            del database2id[database_name]

        # first try to find whether one of the IDs is already in the database
        list_of_abstract_item__database_name_and_database_id = []
        for database_name, database_id in database2id.items():
            this_abstract_item_id = self.__getAbstractIdByIdDict__(
                {database_name:database_id}
                )
            if this_abstract_item_id != None:
                if [database_name, database_id] in self.generalized_items:
                    more_generalized_item_links.append([database_name, database_id])
                    more_generalized_abstract_items.append(this_abstract_item_id)
                else:
                    list_of_abstract_item__database_name_and_database_id.append(
                        [this_abstract_item_id, database_name, database_id]
                        )
        # maybe there are some generalized items referenced
        for database_name,database_id in more_generalized_item_links:
            # remove them from the directly linked items
            del database2id[database_name]
        if len(list_of_abstract_item__database_name_and_database_id) == 0:
            abstract_item_id = None
        else:
            # check the list_of_abstract_item__database_name_and_database_id
            # do they contain a link to an external database which is
            # described differently in the database2id dictionary?
            generalized_abstract_items = []
            for this_abstract_item_id,this_database_name,this_database_id in \
                    list_of_abstract_item__database_name_and_database_id:
                for new_ai_db_name,new_ai_db_id in \
                        self.getAnnotationObjectByAbstractId(this_abstract_item_id).items():
                    if new_ai_db_name in database2id and \
                       database2id[new_ai_db_name] != new_ai_db_id:
                        if this_abstract_item_id not in generalized_abstract_items:
                            generalized_abstract_items.append(this_abstract_item_id)
            # check if there was a generalized_abstract_item found
            for generalized_abstract_item in generalized_abstract_items:
                # create a new abstract item
                new_abstract_item = len(self.abstract_item)
                self.abstract_item.insertRow([item_type])

                # reroute the external database items to this new abstract
                for this_abstract_item_id, this_database_name,this_database_id in \
                        list_of_abstract_item__database_name_and_database_id:
                    if this_abstract_item_id == generalized_abstract_item:
                        # find the entry this_database_name to this_database_id
                        database_table = self.__dict__["external_item_"+this_database_name]
                        for rowIndex in database_table.findAllRowIndices(
                            {"data_id":this_database_id}
                            ):
                            # reroute the abstract item to the new abstract item
                            new_row = copy.deepcopy(database_table.rows[rowIndex])
                            new_row[1] = new_abstract_item
                            database_table.updateRow(rowIndex,new_row)

                # create an is_a relation
                self.abstract_item_relation.insertRow(
                    [generalized_abstract_item, new_abstract_item,
                     "is_a", "automatically generated"]
                    )
                for ai,database_name,database_id in \
                        list_of_abstract_item__database_name_and_database_id:
                    if ai == generalized_abstract_item:
                        # add a new abstract item to the generalized items list
                        self.generalized_items.append([database_name,
                                                       database_id])
                        # remove the database link from the database2id
                        # dictionary which will be inserted later
                        del database2id[database_name]
                # go on with the rest
                list_of_abstract_item__database_name_and_database_id = \
                    [x for x in \
                     list_of_abstract_item__database_name_and_database_id \
                     if x[0] != generalized_abstract_item]
            set_of_abstract_ids_found_from_database_links = \
                set( [x[0] for x in \
                      list_of_abstract_item__database_name_and_database_id]\
                   )
            clean_abstract_ids = \
                list(set_of_abstract_ids_found_from_database_links.difference( \
                         set( generalized_abstract_items ) ))
            if len(clean_abstract_ids) == 1:
                # found unique id
                abstract_item_id = clean_abstract_ids[0]
            elif len(clean_abstract_ids) > 1:
                # here we have found abstract items which are supposed
                # to be together according to this database
                abstract_item_id = min(clean_abstract_ids)
                # set the abstract id of the other external resources
                # to this abstract_item_id
                for database_name, database_id in database2id.items():
                    this_abstract_item_id = self.__getAbstractIdByIdDict__( \
                        {database_name:database_id} \
                        )
                    if this_abstract_item_id != abstract_item_id:
                        external_table = \
                            self.__dict__["external_item_"+database_name]
                        rowIndices =  external_table.findAllRowIndices( \
                                          {"data_id":database_id})
                        if len(rowIndices) == 0:
                            # entry not in database
                            continue
                        rowIndex = rowIndices[0]
                        row = copy.deepcopy(external_table.rows[rowIndex])
                        row[1] = abstract_item_id
                        external_table.updateRow(rowIndex, row)
            else:
                # could not find associated abstract item or an error occured
                abstract_item_id = None



        # do I need to create an abstract item
        if abstract_item_id == None:
            abstract_item_id = len(self.abstract_item.rows)
            self.abstract_item.insertRow([item_type])
        # add external items
        for database_name in database2id:
            id_in_database = database2id[database_name]
            database_table = self.__dict__["external_item_"+database_name]
            # is the entry already in the database
            if not database_table.entryInDatabase({"data_id":id_in_database}):
                database_table.insertRow([id_in_database, abstract_item_id, False])
        # adding the names
        if len(name_list)>0:
            # before adding the name, check if the abstract id already has
            # a favourite description
            this_rows = [x for x in self.string_repr.findAllRows({"abstract_item":abstract_item_id}) if x[2] == True]
            needs_favourite_description = len(this_rows) == 0
            self.string_repr_from_database.insertRow([name_list[0], abstract_item_id, source_database, name_list[0].lower()])

            if not self.string_repr.entryInDatabaseFaster({"repr":name_list[0],
                                                           "abstract_item":abstract_item_id,
                                                           "is_favourite_description":needs_favourite_description,
                                                           "lower_case":name_list[0].lower()}):
                self.string_repr.insertRow([name_list[0], abstract_item_id,
                                            needs_favourite_description,
                                            name_list[0].lower()])

            for name in name_list[1:]:
                #if not self.string_repr.entryInDatabase({"repr":name}):
                if not self.string_repr.entryInDatabaseFaster({"repr":name,
                                                               "abstract_item":abstract_item_id,
                                                               "is_favourite_description":False,
                                                               "lower_case":name.lower()}):
                    self.string_repr.insertRow([name, abstract_item_id, False,
                                                name.lower()])

        # adding the extra relations:
        for related_abstract_item in more_generalized_abstract_items:
            if not self.abstract_item_relation.entryInDatabaseFaster({"from_item":abstract_item_id,
                                                                      "to_item":related_abstract_item,
                                                                      "relation_type":"is_a",
                                                                      "evidence":"automatically generated"}):
                self.abstract_item_relation.insertRow([abstract_item_id,
                                                       related_abstract_item,
                                                       "is_a",
                                                       "automatically generated"])
        for related_abstract_item in more_specialized_abstract_items:
            # the entries in the more_specialized_abstract_items could already have been included into the current item
            if related_abstract_item == abstract_item_id:
                continue
            if not self.abstract_item_relation.entryInDatabaseFaster({"from_item":related_abstract_item,
                                                                      "to_item":abstract_item_id,
                                                                      "relation_type":"is_a",
                                                                      "evidence":"automatically generated"}):
                self.abstract_item_relation.insertRow([related_abstract_item,
                                                       abstract_item_id,
                                                       "is_a",
                                                       "automatically generated"])
        return abstract_item_id

    def addName(self, abstract_item_id, name_list):
        for name in name_list:
            if not self.string_repr.entryInDatabaseFaster({"repr":name,
                                                           "abstract_item":abstract_item_id,
                                                           "is_favourite_description":False,
                                                           "lower_case":name.lower()}):
                self.string_repr.insertRow([name, abstract_item_id, False,
                                            name.lower()])

    def addExternalItem(self, abstract_item_id, database2id):
        # add external items
        for database_name in database2id:
            id_in_database = database2id[database_name]
            database_table = self.__dict__["external_item_"+database_name]
            #
            if database_table.entryInDatabase({"abstract_item":abstract_item_id}):
                r = database_table.findAllRows({"abstract_item":abstract_item_id})
                raise Exception("Duplicate abstract id inserting\n%s\ninto\n%s"%(str(database2id),str(r)))
            # is the entry already in the database
            if not database_table.entryInDatabase({"data_id":id_in_database}):
                database_table.insertRow([id_in_database, abstract_item_id, False])
    
    def addRelationship(self, source, from_id, relation_list):
        """Adds a list of relations for a certain database.
        The relation_list is a list of two elements lists giving
        the relationship and the id."""
        from_abstract_item = self.__getAbstractIdByIdDict__(
            {source:from_id}
            )
        for [relation,to_id] in relation_list:
            to_abstract_item = self.__getAbstractIdByIdDict__(
                {source:to_id}
                )
            if to_abstract_item == None:
                # the to item has to be created in the database
                to_abstract_item = self.addSpecies({source:to_id},[], source)
            self.__addRelationshipByAbstractId__(source, relation, from_abstract_item, to_abstract_item)
            
    def addRelationshipByNames(self, source, relation_type,
                               from_name, to_name):
        from_abstract_id = self.__getAbstractIdByNameList__([from_name])
        to_abstract_id   = self.__getAbstractIdByNameList__([to_name])
        if from_abstract_id == None or to_abstract_id == None:
            raise InconsistentDatabaseError("Could not add relationship between names "+\
                                            from_name+" and "+to_name)
        self.__addRelationshipByAbstractId__(source, relation_type, from_abstract_id, to_abstract_id)
        
    def __addRelationshipByAbstractId__(self, source, relation_type,
                                      from_abstract_id, to_abstract_id):
        if not self.abstract_item_relation.entryInDatabaseFaster({"from_item":from_abstract_id,
                                                            "to_item":to_abstract_id,
                                                            "relation_type":relation_type,
                                                            "evidence":source}):
            self.abstract_item_relation.insertRow([from_abstract_id,
                                                   to_abstract_id,
                                                   relation_type,
                                                   source])
            
    def __getNameByAbstractReactionId__(self, abstract_id):
        def part_of_line(self,rows):
            """Display the part of the line that corresponds to
            all substrates / products / whatever"""
            string_list = []
            for [species_id,foo,role,stoichiometry] in rows:
                species_name = self.__getNameByAbstractId__(species_id)
                if species_name == None:
                    species_name = "unnamed substance"
                string_list.append(
                    (str(stoichiometry)+" ")*(stoichiometry != 1)\
                    +species_name )
            return " + ".join(string_list)
        
        participation_rows = self.takes_part.findAllRows(
            {"reaction_id":abstract_id})
        substrates = part_of_line(self,
                                  [x for x in participation_rows \
                                   if x[2] == "substrate"])
        products = part_of_line(self,
                                [x for x in participation_rows \
                                 if x[2] == "product"])
        modifiers = part_of_line(self,
                                 [x for x in participation_rows \
                                  if x[2] != "substrate" and \
                                  x[2] != "product"])
        return substrates + " <=> " + products + \
               " [ " + modifiers + " ] "
    
    def __getAllNamesByAbstractId__(self,abstract_id):
        rows = self.string_repr.findAllRows({"abstract_item":abstract_id})
        return [x[0] for x in rows]
    
    def __getNameByAbstractId__(self,abstract_id,try_parents=1):
        if abstract_id == None:
            return None
        # get the type of the abstract item
        type_of_abstract_item = self.abstract_item.rows[abstract_id][0]
        if type_of_abstract_item == "reaction":
            return self.__getNameByAbstractReactionId__(abstract_id)
        else:
            rows = self.string_repr.findAllRows({"abstract_item":abstract_id})
            rows = [x for x in rows if x[2] == True]
            
            if len(rows) == 0:
                # no entry found
                if try_parents > 0:
                    relationships = self.__getRelationshipsByAbstractId__(abstract_id,"is_a")
                    # direct search
                    # avoid depth-first-search by setting the recursion depth to 0
                    for related_abstract_id in relationships.getRelatedAbstractIds():
                        name = self.__getNameByAbstractId__(related_abstract_id,0)
                        if name != None:
                            return name
                    # recursive step
                    for related_abstract_id in relationships.getRelatedAbstractIds():
                        name = self.__getNameByAbstractId__(related_abstract_id,try_parents-1)
                        if name != None:
                            return name
                    return None
                else:
                    return None
            return rows[0][0]
        
    def __getDirectNameByAbstractId(self, abstract_id, database):
        # returns the standard name from a certain database
        rows = self.string_repr_from_database.findAllRows({"abstract_item":abstract_id,
                                                           "database":database})
        if len(rows) > 0:
            direct_name = rows[0][0]
            if direct_name != None and direct_name != "":
                return direct_name
        # otherwise get a name the normal way
        return self.__getNameByAbstractId__(abstract_id)
        
    def __getNameByIdDict__(self,database2id, try_direct_name=False):
        abstract_id = self.__getAbstractIdByIdDict__(database2id)
        if abstract_id == None:
            return None
        if try_direct_name and len(list(database2id.keys())) == 1:
            return self.__getDirectNameByAbstractId(abstract_id, list(database2id.keys())[0])
        return self.__getNameByAbstractId__(abstract_id)
    
    def __getAbstractIdByNameList__(self,name_list):
        possible_abstract_ids = set([y for y in \
                                     [self.__getAbstractIdByName__(x) for x in name_list] \
                                     if y != None])
        if len(possible_abstract_ids) > 1:
            return None
            #raise InconsistentDatabaseError("More than one semanticSBML ID for a set "+\
            #                                "of names was found!\n"+"\n".join(name_list))
        elif len(possible_abstract_ids) == 1:
            return possible_abstract_ids.pop()
        else:
            return None
        
    def __getAbstractIdListByName__(self,name):
        rows = self.string_repr.findAllRows({"lower_case":name.lower()})
        if len(rows) == 0:
            return None
        # check if we have a favourite description
        favourite_rows = [x for x in rows if x[2] == True]
        if len(favourite_rows)>0:
            rows = favourite_rows
        return list(set([row[1] for row in rows]))
    
    def __getAbstractIdByName__(self,name):
        aid_list = self.__getAbstractIdListByName__(name)
        if aid_list == None:
            return None
        return aid_list[0]
    
    def __getAbstractIdByIdDict__(self,database2id):
        for database_name in database2id:
            id_in_database = database2id[database_name]
            if database_name not in self.list_of_internal_names_for_external_resources:
                if database_name not in dict(ext_db_names):
                    continue
                    #raise InconsistentDatabaseError("Database "+database_name+\
                    #                                " unknown to external_database_list.")
                else:
                    database_name = dict(ext_db_names)[database_name]
            database_table = self.__dict__["external_item_"+database_name]
            rows = database_table.findAllRows({"data_id":id_in_database})
            if len(rows) == 1:
                return rows[0][1]
        return None
    
    def getAnnotationObjectByIdDict(self,database2id):
        abstract_id = self.__getAbstractIdByIdDict__(database2id)
        if abstract_id == None:
            return None
        return self.getAnnotationObjectByAbstractId(abstract_id)

    def getAnnotationObjectByAbstractId(self,abstract_id):
        if abstract_id in self.abstract_id_to_annotation_object:
            return self.abstract_id_to_annotation_object[abstract_id]
        description_dict = {}
        for display_name, internal_name in ext_db_names:
            database_table = self.__dict__["external_item_"+internal_name]
            rows = database_table.findAllRows({"abstract_item":abstract_id})
            if len(rows)>0:
                description_dict[internal_name] = rows[0][0]
        if len(description_dict) == 0:
            returnme = None
        else:
            returnme = AnnotationObject(description_dict)
        self.abstract_id_to_annotation_object[abstract_id] = returnme
        return returnme

    def getAnnotationObjectListByName(self,name):
        abstract_id_list = self.__getAbstractIdListByName__(name)
        if abstract_id_list == None:
            return None
        return [self.getAnnotationObjectByAbstractId(aid) for aid in abstract_id_list]

    def getAnnotationObjectByName(self,name):
        abstract_id = self.__getAbstractIdByName__(name)
        if abstract_id == None:
            return None
        return self.getAnnotationObjectByAbstractId(abstract_id)

    def getListOfAllNames(self):
        """Returns all lower case names in the database."""
        if self.list_of_all_names == None:
            self.list_of_all_names = [x[3] for x in self.string_repr.rows]
        return self.list_of_all_names

    def getCloseNames(self,name,fuzziness):
        from . import mydifflib
        return mydifflib.get_close_matches(name.lower(),self.getListOfAllNames(),10,fuzziness)

    def getSuperstringNames(self, name, fuzziness):
        max_len = float(len(name))/fuzziness
        results = []
        for other_name in self.getListOfAllNames():
            if other_name.find(name) != -1:
                results.append(other_name)
        results.sort(key=cmp_to_key(lambda x,y: len(x)-len(y)))
        return results

    def isNameInDatabase(self,name):
        return self.string_repr.entryInDatabase({"lower_case":name.lower()})

    def __getRelationshipsByAbstractIdWithoutCaching__(self,abstract_id,
                                         relation_type="all",direction="both", source=None):
        """This function will give back the entries in the
        abstractItemRelation table.

        @return: Returns a Relation Object
        @rtype: RelationObject"""
        to_rows = []
        from_rows = []
        if relation_type == "all":
            if direction != "from":
                to_rows = self.abstract_item_relation.findAllRows({"from_item":abstract_id})
            if direction != "to":
                from_rows = self.abstract_item_relation.findAllRows({"to_item":abstract_id})
        else:
            if direction != "from":
                to_rows = self.abstract_item_relation.findAllRows({"relation_type":relation_type,
                                                                   "from_item":abstract_id})
            if direction != "to":
                from_rows = self.abstract_item_relation.findAllRows({"relation_type":relation_type,
                                                                     "to_item":abstract_id})

        
        relation_object = RelationObject(abstract_id)
        for rows,direction_is_reverse,index_new_abstract_id in [[from_rows,True,0],
                                                                [to_rows,False,1]]:
            for row in rows:
                # check whether we restrict our search to a certain source
                # and then skip other rows
                if source != None and source != row[3]:
                    continue
                related_abstract_id = row[index_new_abstract_id]
                relation_type = row[2]
                related_database2id = self.getAnnotationObjectByAbstractId(related_abstract_id)
                if related_database2id != None:
                    related_database2id = related_database2id.toDatabase2id()
                relation_object.addRelation(related_abstract_id, related_database2id,
                                            relation_type, direction_is_reverse)
        return relation_object
    
    def __getRelationshipsByAbstractId__(self,abstract_id,
                                         relation_type="all",direction="both"):
        """This function will give back the entries in the
        abstractItemRelation table.

        @return: Returns a Relation Object
        @rtype: RelationObject"""
        rows = []
        relation_object = RelationObject(abstract_id)

        if not abstract_id in self.relation_cache:
            return relation_object
        this_relation_cache = self.relation_cache[abstract_id]
        if relation_type == 'all':
            relation_keys = list(this_relation_cache.keys())
        else:
            relation_keys = [relation_type] if relation_type in this_relation_cache else []
        for this_relation_type in relation_keys:
            direction_to_aid_dict = this_relation_cache[this_relation_type]
            if direction == 'both':
                direction_keys = list(direction_to_aid_dict.keys())
            else:
                direction_keys = [direction] if direction in direction_to_aid_dict else []
            for this_direction in direction_keys:
                aid_to_none = direction_to_aid_dict[this_direction]
                for aid in aid_to_none:
                    rows.append( [abstract_id, aid, this_relation_type, this_direction] )
        
        for row in rows:
            related_abstract_id = row[1]
            relation_type = row[2]
            related_database2id = self.getAnnotationObjectByAbstractId(related_abstract_id)
            if related_database2id != None:
                related_database2id = related_database2id.toDatabase2id()
            relation_object.addRelation(related_abstract_id, related_database2id,
                                        relation_type, row[3]=="to")
        return relation_object

    def __getDepth__(self, abstract_id):
        raids = self.__getRelatedAbstractIdsByAbstractId__(abstract_id, relation_type="is_a", direction="from")
        if len(raids) == 0:
            return 0
        else:
            return max([self.__getDepth__(raid) for raid in raids]) + 1

    def __getRelatedAbstractIdsByAbstractId__(self, abstract_id, relation_type="all",
                                              depth=1, direction="both"):
        """Returns a set of related abstract ids."""
        if depth == 1:
            return set( self.__getRelationshipsByAbstractId__(abstract_id,
                                                              relation_type,
                                                              direction)\
                        .getRelatedAbstractIds() )
        else:
            # find all children
            children_set = set( self.__getRelationshipsByAbstractId__(abstract_id,
                                                                      relation_type,
                                                                      direction)\
                                .getRelatedAbstractIds() )
            for child in list(children_set):
                # recursion
                children_set = children_set.union(set( \
                    self.__getRelatedAbstractIdsByAbstractId__(child,
                                                               relation_type,
                                                               depth-1,
                                                               direction) ))
            # add yourself
            return children_set.union([abstract_id])

    def getListOfRelatedIdDictsByIdDict(self, database2id, relation_type="all",
                                        depth=1, direction="both"):
        abstract_id = self.__getAbstractIdByIdDict__(database2id)
        if abstract_id == None:
            return []
        list_of_related_abstract_ids = list(\
            self.__getRelatedAbstractIdsByAbstractId__(abstract_id,
                                                       relation_type,
                                                       depth,
                                                       direction))
        return [self.getAnnotationObjectByAbstractId(x).toDatabase2id() \
                for x in list_of_related_abstract_ids]

    def getRelationshipsByName(self, name, relation_type="all", direction="both"):
        abstract_id = self.__getAbstractIdByName__(name)
        if abstract_id == None:
            return None
        return self.__getRelationshipsByAbstractId__(abstract_id,
                                                     relation_type,
                                                     direction)

    def getRelationshipsByIdDict(self, database2id,
                                 relation_type="all", direction="both"):
        abstract_id = self.__getAbstractIdByIdDict__(database2id)
        if abstract_id == None:
            return None
        return self.__getRelationshipsByAbstractId__(abstract_id,
                                                     relation_type,
                                                     direction)

    def __isRelatedByAbstractId__(self, abstract_id1,
                                  abstract_id2, directed=False):
        direction_right_rows = self.abstract_item_relation\
                               .findAllRows({"from_item":abstract_id1,
                                             "to_item":abstract_id2})
        direction_left_rows = self.abstract_item_relation\
                              .findAllRows({"from_item":abstract_id2,
                                            "to_item":abstract_id1}) \
                              if not directed else []
        if len(direction_right_rows)+len(direction_left_rows) > 0:
            return True
        else:
            return False

    def isRelatedByIdDict(self, database2id1,
                          database2id2, directed=False):
        abstract_id1 = self.__getAbstractIdByIdDict__(database2id1)
        abstract_id2 = self.__getAbstractIdByIdDict__(database2id2)
        if abstract_id1 == None or abstract_id2 == None:
            return False
        return self.__isRelatedByAbstractId__(abstract_id1,
                                              abstract_id2,
                                              directed)



    def __getReactionAISetByAbstractId__(self,abstract_id):
        return set([x[1] for x in \
                    self.takes_part.findAllRows({"species_id":abstract_id})])

    def __getReactionAISetByListOfAbstractIds__(self,list_of_abstract_ids):
        reaction_ai_set = self.__getReactionAISetByAbstractId__(list_of_abstract_ids[0])
        for abstract_id in list_of_abstract_ids[1:]:
            reaction_ai_set = reaction_ai_set.intersection( \
                self.__getReactionAISetByAbstractId__(abstract_id) )
        return reaction_ai_set

    def getReactionAISetByListOfIdDicts(self,list_of_database2ids):
        abstract_ids = []
        for database2id in list_of_database2ids:
            abstract_id = self.__getAbstractIdByIdDict__(database2id)
            if abstract_id == None:
                return None
            abstract_ids.append(abstract_id)
        return self.__getReactionAISetByListOfAbstractIds__(abstract_ids)
        

if __name__ == "__main__":
    d = Database()
    d.load()
    print(d.__getNameByIdDict__({"chebi":"CHEBI:15343"}))
    #database2id1 = {"kegg_compound":"C00005","chebi":"CHEBI:16474"}
    #database2id2 = {"kegg_compound":"C00009","chebi":"CHEBI:18367"}
    #database2id3 = {"cas":"53-57-6"}
    #database2id4 = {"cas":"14265-44-2"}
    #print bigd.__getNameByIdDict__(database2id1)
    #print bigd.__getNameByIdDict__(database2id2)
    #print bigd.getRelationshipsByIdDict(database2id1)
    #print bigd.getRelationshipsByIdDict(database2id2)
    #print bigd.__getNameByIdDict__(database2id3)
    #print bigd.__getNameByIdDict__(database2id4)

