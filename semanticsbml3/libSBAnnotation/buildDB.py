#!/usr/bin/env python

# Allysons telomere data

# SABIO-RK Reaction # webservice
# Reactome
# Pathway Commons
# SGD
# Saccharomyces genome database pathways
# WikiPathways
# BiGG Database
# Biocyc

import re, xml.sax
from . import stderr

tempvar = []

class2size = {"NcbiTaxonomyNodesRead":1532,
              "NcbiTaxonomyNamesRead":1532,
              "ReactomeRead":11967,
              "GeneOntologyRead":29027,
              "ChebiRead":22274,
              "SboRead":422,
              "KeggTaxonomyRead":False,
              "KeggGenomeRead":1075,
              "KeggCompoundRead":15441,
              "KeggDrugRead":8904,
              "KeggEnzymeRead":5061,
              "KeggReactionRead":7910,
              "SGDCrossRefRead":71085,
              "InterproRead":20329,
              "AllysonNodesRead":False,
              "AllysonNamesRead":False,
              "IntEnzRead":5172,
              "BioGridRead":17424000,
              "DruggableRead":130,
              "MatadorRead":15844,
              "DrugbankRead":4774,
              "StitchActionsRead":24422984,
              "StitchChemicalSourcesRead":59189164,
              "StringProteinAliasesRead":20624144,
              "BindingDBTXTRead":251447,
              "BindingDBSDFRead":314527}

class Read:
    """Basic class for the processing of new databases. This abstract class has to be derived from twice.
    The first derived class is the implementation of the processing of a particular type of input format.
    This class has to implement the read function which is processing the file_contents and producing an entry
    dictionary for every entry in the database. This entry is then processed in the processEntry function
    which has to be implemented in the second derived class."""
    def __init__(self, database, file_content, size_known=False):
        """Initialize the flat file processing.

        @param database: database to be filled
        @type database: database
        @param file_content: flat file to be processed
        @type: string"""
        self.database = database
        self.current_species = 0
        self.class_name = str(self.__class__).split(".")[-1]
        self.size_known = size_known
        if self.size_known == False and self.class_name in class2size:
            self.size_known = class2size[self.class_name]
        stderr.downloadProgress(self.size_known, 0, True, False, url=str(self.__class__).split(".")[-1][:-4])
        self.read(file_content)
        self.database.resetGeneralizedItems()
        stderr.downloadProgress(False,self.current_species,False,True)
    def getDatabaseName(self, name):
        """Converts the names occuring the the flat files to internal names for databases.

        @param name: name in the flat file
        @type name: string
        @return: name in the database
        @rtype: string"""
        newname = name.lower().replace("-","_")
        if newname == "pubchem":
            newname = "pubchem_substance"
        return newname
    def reportReadEntry(self):
        """Function being called when an entry in the flat file has been read."""
        self.current_species = self.current_species + 1
        stderr.downloadProgress(self.size_known, self.current_species, url=str(self.__class__))
        return
    def read(self, file_content):
        """Function converting the file content and adding new entries to the database with the processEntry function.
        For this purpose it creates an \"entry\" dictionary for each database entry and feeds this into the processEntry function."""
        raise NotImplementedError
    def safeProcessEntry(self, entry):
        try:
            self.processEntry(entry)
        except Exception as e:
            print(repr(e))
            #raise
    def processEntry(self, entry):
        """Add an antry dictionary to the database."""
        raise NotImplementedError


class TabLineRead(Read):
    def __init__(self,database,file_content, size_known=False):
        Read.__init__(self,database,file_content, size_known=size_known)
    def read(self,file_content):
        if 'newline' in dir(file_content):
            func = file_content.newline
            args = []
        else:
            func = file_content.split
            args = ['\n']
        for line in func(*args):
            if len(line) == 0 or line.startswith("#"):
                continue
            self.safeProcessEntry(line.split("\t"))
            self.reportReadEntry()
class SGDCrossRefRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        self.ignore = ["NCBINCBI protein GI", "NCBIGene ID", "NCBIDNA accession ID", "NCBIRefSeq Accession",
                       "NCBIRefSeq protein version ID", "GenBank/EMBL/DDBJDNA accession ID",
                       "GenBank/EMBL/DDBJDNA version ID", "GenBank/EMBL/DDBJProtein version ID",
                       "EBIUniProt/TrEMBL ID", "EBIUniProt/Swiss-Prot ID", "EBIUniParc ID",
                       "EUROSCARFGene ID", "CGDGene ID", "TCDBTC number", "DIPGene ID", "GermOnlineGene ID"]
        self.source2name = {"IUBMBEC number":"ec",
                            "TCDBTC number":"tcdb",
                            "DIPGene ID":"dip",
                            "GermOnlineGene ID":"germonline",
                            "BioGRIDGene ID":"biogrid",
                            "MetaCycPathway ID":"metacyc"}
        self.source2relation = {"BioGRIDGene ID":None,
                                "IUBMBEC number":"hasFunction"} # hasFunction hasProcess isLocated
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        [xrefid, xrefsource, xreftype, systematic_name, sgdid, standard_name] = entry
        aid = self.database.__getAbstractIdByIdDict__({'sgd':sgdid})
        if aid is None:
            aid = self.database.addSpecies({'sgd':sgdid}, [standard_name, systematic_name], 'sgd')
        if xrefsource+xreftype in self.ignore:
            return
        # get the reference ID
        sourcename = self.source2name[xrefsource+xreftype]
        relation = self.source2relation[xrefsource+xreftype]
        raid = self.database.__getAbstractIdByIdDict__({sourcename:xrefid})
        if raid is None:
            if relation is None:
                # add the reference directly to the aid
                self.database.addExternalItem(aid, {sourcename:xrefid})
            else:
                # create a new species
                raid = self.database.addSpecies({sourcename:xrefid}, [], 'sgd')
        if relation is not None:
            # add a relation link between aid and raid
            self.database.__addRelationshipByAbstractId__('sgd', relation, aid, raid)
        
        
        
class AllysonNodesRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        [productId, organism, canonicalName, aliases] = entry
        aid = self.database.addSpecies({"allyson":productId},[canonicalName],"allyson")
        # add the organism
        tax_aids = self.database.__getAbstractIdListByName__(organism)
        if tax_aids != None and len(tax_aids) == 1:
            tax_aid = tax_aids[0]
            self.database.__addRelationshipByAbstractId__("allyson", "hasTaxonomy", aid, tax_aid)
        # process the aliases
class AllysonEdgesRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        [source, type, target, weight, foo1, foo2, foo3] = entry
        dbname = "allyson"
        source_abstract_id = self.database.__getAbstractIdByIdDict__({dbname:source})
        if source_abstract_id == None:
            source_abstract_id = self.database.addSpecies({"allyson":source},[],"allyson")
        target_abstract_id = self.database.__getAbstractIdByIdDict__({dbname:target})
        if target_abstract_id == None:
            target_abstract_id = self.database.addSpecies({"allyson":target},[],"allyson")
        self.database.__addRelationshipByAbstractId__(weight, type, source_abstract_id, target_abstract_id)
class ReactomeRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        [reactome_id, long_name, foo, cross_id_line] = entry
        if "[" in long_name:
            short_name = long_name.split("[")[0].strip()
            compartment = long_name.split("[")[1].strip().strip("]")
        else:
            compartment = None
        kegg_compound_id = None
        for cross_id in cross_id_line.split("|"):
            if cross_id.startswith("C") and len(cross_id) == 6:
                # is KEGG Compound ID
                kegg_compound_id = cross_id
        # adding the species
        self.database.addSpecies({"reactome":reactome_id},[long_name],"reactome")
        

class KeggRead(Read):
    def __init__(self,database,file_content, size_known=False):
        Read.__init__(self,database,file_content, size_known=size_known)
    def read(self,file_content):
        entry = {}
        current_header = None
        for line in file_content.split("\n"):
            if len(line)==0:
                continue
            if line.startswith("///"):
                self.safeProcessEntry(entry)
                self.reportReadEntry()
                entry = {}
            elif line.startswith(" "):
                entry[current_header] = entry[current_header]+"\n"+line.strip()
            else:
                try:
                    current_header = line.split()[0]
                except IndexError:
                    print(">"+line+"<")
                    raise
                entry[current_header] = (" ".join(line.split()[1:])).strip()


class KeggTaxonomyTreeRead(Read):
    def __init__(self, database, file_content, size_known=False):
        Read.__init__(self, database, file_content, size_known=size_known)
    def read(self, file_content):
        current_hierarchie = []
        for line in file_content.split("\n"):
            if len(line) == 0:
                continue
            if line.startswith("#"):
                # hierarchical entry
                depth = line.count("#")-1
                name = line.strip("#").strip()
                entry = {}
                entry["NAME"] = [name]
                current_hierarchie = current_hierarchie[:depth]
                # put the name and its parent into the database
                if len(current_hierarchie)>0:
                    entry["PARENT"] = current_hierarchie[-1]
                self.safeProcessEntry(entry)
                current_hierarchie = current_hierarchie+[name]
            else:
                split_line = line.split("\t")
                t_number    = split_line[0]
                kegg_tax_id = split_line[1]
                short_name  = split_line[2]
                long_name   = split_line[3]
                entry = {}
                entry["ENTRY"] = kegg_tax_id
                name_list = []
                # process long name
                is_est = "(EST)" in long_name
                if is_est:
                    long_name = " ".join(long_name.split()[:-1])
                # check for a trivial name
                if "(" in long_name:
                    trivial_name = long_name.split("(")[1].strip().strip(")")
                    name_list.append(trivial_name)
                    long_name = long_name.split("(")[0].strip()
                name_list.append(long_name)
                name_list.append(short_name)
                entry["NAME"] = name_list
                entry["PARENT"] = current_hierarchie[-1]
                self.safeProcessEntry(entry)
            self.reportReadEntry()

class OboRead(Read):
    def __init__(self,database,file_content, size_known=False):
        self.database_name_translation = {"CAS Registry Number":"cas","KEGG COMPOUND":"kegg_compound","KEGG DRUG":"kegg_drug","KEGG REACTION":"kegg_reaction","Gmelin Registry Number":"gmelin","MSDchem":"msdchem","Beilstein Registry Number":"beilstein","UM-BBD compID":"um_bbd","MolBase":"molbase","DrugBank":"drugbank","PDB":"pdb","EC":"ec","MetaCyc":"metacyc","Reactome":"reactome","UM-BBD_enzymeID":"um_bbd","UM-BBD_pathwayID":"um_bbd","UM-BBD_reactionID":"um_bbd","AraCyc":"aracyc", "PubMed citation":"pubmed", "NIF_Subcellular":"nif_subcellular","CORUM":"corum","PMID":"pmid"}
        self.ignore_databases = ["Patent","COMe","LIPID MAPS instance","LIPID MAPS class","RESID","WebElements","Wikipedia","GOC","TC","IMG","KEGG GLYCAN","UniPathway","PDBeChem","PO"]
        Read.__init__(self,database,file_content, size_known=size_known)
    def read(self,file_content):
        entry = {}
        current_header = None
        for line in file_content.split("\n"):
            if line.startswith("[T"):
                if entry != {}:
                    self.safeProcessEntry(entry)
                    self.reportReadEntry()
                entry = {}
            elif line.startswith("\n"):
                continue
            else:
                current_header = line.split(":")[0]
                rest_of_line = (":".join(line.split(":")[1:])).strip()
                if current_header not in entry:
                    entry[current_header] = []
                entry[current_header].append(rest_of_line)
        self.safeProcessEntry(entry)









class KeggCompoundRead(KeggRead):
    def __init__(self, database, file_content, size_known=False):
        KeggRead.__init__(self, database, file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("kegg_compound")
    def processEntry(self,entry):
        kegg_id = entry["ENTRY"].split()[0]
        dblinks = {"kegg_compound":kegg_id}
        if "DBLINKS" in entry:
            database_name = None
            for line in entry["DBLINKS"].split("\n"):
                # strip off the name of the database
                if ":" in line:
                    database_name = self.getDatabaseName(line.split(":")[0].strip())
                    database_id = (" ".join(line.split(":")[1:])).strip()
                # or keep the database name of the last row
                else:
                    database_id = line.strip()
                # search the database_id for spaces, they denote a list of IDs
                if " " in database_id:
                    database_id = database_id.split()
                # add the database_id(s) to the dblinks
                if database_name not in dblinks:
                    dblinks[database_name] = database_id
                else:
                    # we are appending new values to the dblinks[database_name] list
                    if type(dblinks[database_name]) != type([]):
                        dblinks[database_name] = [dblinks[database_name]]
                    if type(database_id) != type([]):
                        database_id = [database_id]
                    dblinks[database_name] = dblinks[database_name]+database_id
        # extracting the formulas for the kegg entry
        if "FORMULA" in entry:
            dblinks["kegg_formula"] = entry["FORMULA"]
        # put a CHEBI: in front of every chebi ID as this is missing in KEGG Compound
        if "chebi" in dblinks:
            if type(dblinks["chebi"]) == type([]):
                for i in range(len(dblinks["chebi"])):
                    if not dblinks["chebi"][i].startswith("CHEBI:"):
                        dblinks["chebi"][i] = "CHEBI:"+dblinks["chebi"][i]
            else:
                if not dblinks["chebi"].startswith("CHEBI:"):
                    dblinks["chebi"] = "CHEBI:"+dblinks["chebi"]
        if "NAME" in entry:
            names = "".join(entry["NAME"].split("\n"))
            names = names.split(";")
        else:
            names = []
        self.database.addSpecies(dblinks,names,"kegg_compound")


class KeggEnzymeRead(KeggRead):
    def __init__(self,database,file_content, size_known=False):
        KeggRead.__init__(self,database,file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("kegg_enzyme")
    def processEntry(self,entry):
        if "ENTRY" not in entry:
            return
        ecnumber = entry["ENTRY"].split()[1]
        # allow only precise EC numbers
        if "-" in ecnumber:
            return
        namelist = entry["NAME"].replace("\n","").split(";")
        # prepend the sysname if available
        if "SYSNAME" in entry:
            namelist.insert(0," ".join([x.strip() for x in entry["SYSNAME"].split("\n")]))
        # TODO: GENES STRUCTURES(PDB)
        dblinks = {"ec":ecnumber}
        if "DBLINKS" in entry:
            for line in entry["DBLINKS"].split("\n"):
                database_name = self.getDatabaseName(line.split(":")[0].strip())
                database_id = line.split(":")[1].strip()
                if database_id != ecnumber:
                    dblinks[database_name] = database_id
        self.database.addSpecies(dblinks,namelist, "kegg_enzyme")

class KeggDrugRead(KeggRead):
    def __init__(self,database,file_content, size_known=False):
        KeggRead.__init__(self,database,file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("kegg_drug")
    def processEntry(self,entry):
        dblinks_orig = {"kegg_drug":entry["ENTRY"].split()[0]}
        dblinks = {}
        namelist = []
        for name in entry["NAME"].replace("\n","").split(";"):
            if "(" in namelist:
                namelist.append(name.split("(")[0].strip())
            else:
                namelist.append(name)
        # TODO: TARGET
        if "DBLINKS" in entry:
            database_name = None
            for line in entry["DBLINKS"].split("\n"):
                # strip off the name of the database
                if ":" in line:
                    database_name = self.getDatabaseName(line.split(":")[0].strip())
                    database_id = (" ".join(line.split(":")[1:])).strip()
                # or keep the database name of the last row
                else:
                    database_id = line.strip()
                # search the database_id for spaces, they denote a list of IDs
                if " " in database_id:
                    database_id = database_id.split()
                # add the database_id(s) to the dblinks
                if database_name not in dblinks:
                    dblinks[database_name] = database_id
                else:
                    # we are appending new values to the dblinks[database_name] list
                    if type(dblinks[database_name]) != type([]):
                        dblinks[database_name] = [dblinks[database_name]]
                    if type(database_id) != type([]):
                        database_id = [database_id]
                    dblinks[database_name] = dblinks[database_name]+database_id
        # put a CHEBI: in front of every chebi ID as this is missing in KEGG Compound
        if "chebi" in dblinks:
            if type(dblinks["chebi"]) == type([]):
                for i in range(len(dblinks["chebi"])):
                    if not dblinks["chebi"][i].startswith("CHEBI:"):
                        dblinks["chebi"][i] = "CHEBI:"+dblinks["chebi"][i]
            else:
                if not dblinks["chebi"].startswith("CHEBI:"):
                    dblinks["chebi"] = "CHEBI:"+dblinks["chebi"]
        abstract_id = self.database.addSpecies(dblinks_orig,namelist, "kegg_drug")
        # add relations for the dblinks dict
        for dbname,id in dblinks.items():
            if type(id) == type([]):
                continue
            parent_abstract_id = self.database.__getAbstractIdByIdDict__({dbname:id})
            if parent_abstract_id != None and parent_abstract_id != abstract_id:
                # add the relation
                self.database.__addRelationshipByAbstractId__("kegg_drug", "is_a", abstract_id, parent_abstract_id)
        

class KeggTaxonomyRead(KeggTaxonomyTreeRead):
    def __init__(self, database, file_content, size_known=False):
        KeggTaxonomyTreeRead.__init__(self, database, file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("kegg_taxonomy")
    def processEntry(self, entry):
        dblinks = {}
        if "ENTRY" in entry:
            dblinks["kegg_organism"] = entry["ENTRY"]
        self.database.addOrganism(dblinks,entry["NAME"],"taxonomy")
        if "PARENT" in entry:
            # add is_a relationship by names
            current_name = entry["NAME"][0]
            parent_name = entry["PARENT"]
            self.database.addRelationshipByNames("kegg_organism",
                                                 "is_a",
                                                 current_name,
                                                 parent_name)

class KeggGenomeRead(KeggRead):
    def __init__(self, database, file_content, size_known=False):
        KeggRead.__init__(self, database, file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("kegg_genome")
    def processEntry(self, entry):
        dblinks = {"kegg_organism":entry["ENTRY"]}
        name_list = []
        definition = entry["DEFINITION"]
        # remove newlines and white spaces from definition
        definition = " ".join([x.strip() for x in definition.split("\n")])
        if "(" in definition:
            name_list.append(definition.split("(")[0].strip())
            name_list.append(definition.split("(")[1].strip(")"))
        else:
            name_list.append(definition)
        if "NAME" in entry:
            for name in entry["NAME"].split(","):
                name = name.strip()
                # check if the name is a taxonomy number
                try:
                    int(name)
                    dblinks["taxonomy"] = name
                    # check if the link is in the database
                    rows = self.database.external_item_taxonomy.findAllRows({"data_id":name})
                    if len(rows)>0:
                        # entry in taxonomy found
                        # make a link to Kegg Organism
                        self.database.addOrganism(dblinks,name_list, "kegg_organism")
                except ValueError:
                    # name is not a taxonomy number
                    if name.lower() not in name_list:
                        # the file also contains additional upper case names, WTF
                        name_list.append(name)
        entry_abstract_id = self.database.addOrganism(dblinks,name_list, "kegg_organism")
        # BREAK FROM HERE
        return
        # walk through the taxonomy
        if "TAXONOMY" in entry:
            taxonomy_line = entry["TAXONOMY"].split("LINEAGE")[1].strip().replace("\n"," ")
            taxonomy_list = [x.strip() for x in taxonomy_line.split(";")]
            for index in range(len(taxonomy_list)):
                tax_name = taxonomy_list[index]
                # test whether this name exists
                if self.database.__getAbstractIdByNameList__([tax_name]) == None:
                    # name does not exist
                    # create organism
                    tax_abstract_id = self.database.addOrganism({},[tax_name], "kegg_organism")
                    # add a connection to the parent
                    if index != 0:
                        parent_name = taxonomy_list[index-1]
                        self.database.addRelationshipByNames("kegg_organism",
                                                             "is_a",
                                                             tax_name,
                                                             parent_name)
            # add a relationship from the bottom of the taxonomy_list to our current entry
            self.database.addRelationshipByNames("kegg_organism", "is_a",
                                                 name_list[0], taxonomy_list[-1])
        

class KeggReactionRead(KeggRead):
    def __init__(self,database,file_content, size_known=False):
        KeggRead.__init__(self,database,file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("kegg_reaction")
    def __processEquationSide__(self,side):
        substrates = []
        for part in side.split(" + "):
            part = part.strip()
            if " " in part:
                [stoich,kegg_id] = part.split()
                try:
                    stoich = int(stoich)
                except ValueError:
                    # convert something with n into a number,
                    # n will be treated as 1
                    stoich = eval(stoich.replace("n","1"))
            else:
                stoich = 1
                kegg_id = part
            if "(" in kegg_id:
                kegg_id = kegg_id.split("(")[0]
            substrates.append([{"kegg_compound":kegg_id},[],stoich])
        return substrates
    def __processEquation__(self,equation):
        return [self.__processEquationSide__(x) for x in equation.split("<=>")]
    def processEntry(self,entry):
        kegg_id = entry["ENTRY"].split()[0]
        enzyme_name_list = []
        if "NAME" in entry:
            for name in entry["NAME"].split("\n"):
                enzyme_name_list.append(name)
        equation_line = " ".join([x.strip() for x in entry["EQUATION"].split("\n")])
        [substrates,products] = self.__processEquation__(equation_line)
        enzyme_dblinks = {}
        if "ENZYME" in entry and "-" not in entry["ENZYME"] and entry["ENZYME"].count(".") < 4:
            enzyme_dblinks["ec"] = entry["ENZYME"]
        enzyme_list = []
        if len(enzyme_name_list)>0 or len(enzyme_dblinks)>0:
            enzyme_list.append([enzyme_dblinks,enzyme_name_list])
        self.database.addReaction({"kegg_reaction":kegg_id},substrates,products,[],enzyme_list, "kegg_reaction")


class ChebiRead(OboRead):
    def __init__(self,database,file_content, size_known=False):
        OboRead.__init__(self,database,file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("chebi")
    def processEntry(self,entry):
        if "is_obsolete" in entry or "is_transitive" in entry or "is_cyclic" in entry:
            return
        if "id" not in entry:
            # some internal definitions
            return
        chebi_id = entry["id"][0]
        names = entry["name"]
        if "synonym" in entry:
            for synonym_line in entry["synonym"]:
                # check whether this line is a formula or SMILE or InChI
                is_formula_or_smile = False
                #for part_of_line in "FORMULA".split():
                #    if part_of_line in synonym_line:
                #        is_formula_or_smile = True
                #        break
                # the line looks like
                # "name" type_of_relation source
                if not is_formula_or_smile:
                    name = synonym_line.split("\"")[1]
                    if name not in names:
                        names.append(name)
        relationships = []
        if "is_a" in entry:
            for line in entry["is_a"]:
                relationships.append(["is_a",line])
        if "relationship" in entry:
            for line in entry["relationship"]:
                relationships.append(line.split())
        database2id = {"chebi":chebi_id}
        if "xref" in entry:
            for line in entry["xref"]:
                source = line.split(":")[0]
                foreign_id = line.split(":")[1].split()[0]
                database_description = line.split(":")[1].split("\"")[1]
                if database_description not in self.ignore_databases:
                    if self.database_name_translation[database_description] in database2id:
                        # check if the entry is already in the database2id dictionary
                        if type(database2id[self.database_name_translation[database_description]]) != type([]):
                            if database2id[self.database_name_translation[database_description]] == foreign_id:
                                continue
                            database2id[self.database_name_translation[database_description]] = [database2id[self.database_name_translation[database_description]]]
                        if foreign_id not in database2id[self.database_name_translation[database_description]]:
                            database2id[self.database_name_translation[database_description]].append(foreign_id)
                    else:
                        database2id[self.database_name_translation[database_description]] = foreign_id
        abstract_id = self.database.addSpecies(database2id,names, "chebi")
        self.database.addRelationship("chebi",chebi_id,relationships)
        if "alt_id" in entry:
            for alt_id in entry["alt_id"]:
                #self.database.external_item_chebi.insertRow([alt_id,abstract_id,False])
                # add a new species with this alt_id
                alt_aid = self.database.addSpecies({"chebi":alt_id},[], "chebi")
                # link this alt_aid with an is_a relationship
                self.database.addRelationship("chebi",chebi_id,[["is_a",alt_id]])


class GeneOntologyRead(OboRead):
    def __init__(self,database, file_content, size_known=False):
        OboRead.__init__(self,database,file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("go")
    def processEntry(self,entry):
        if "is_obsolete" in entry or "is_transitive" in entry or "is_cyclic" in entry:
            return
        if "id" not in entry or "namespace" not in entry:
            return
        go_id = entry["id"][0]
        if "namespace" not in entry:
            print(entry)
        namespace = entry["namespace"][0]
        names = entry["name"]
        if "synonym" in entry:
            for line in entry["synonym"]:
                names.append(line.split("\"")[1])
        relationships = []
        if "is_a" in entry:
            for line in entry["is_a"]:
                relationships.append(["is_a",line.split("!")[0].strip()])
        if "relationship" in entry:
            for line in entry["relationship"]:
                relation_type = line.split()[0]
                relation_id = " ".join(line.split()[1:]).split("!")[0].strip()
                relationships.append([relation_type,relation_id])
        database2id = {"go":go_id}
        if "xref" in entry:
            for line in entry["xref"]:
                foreign_id = line.split(":")[1]
                database_description = line.split(":")[0]
                if database_description not in self.ignore_databases:
                    # skip strange ec numbers 
                    if database_description == "KEGG":
                        if foreign_id.startswith("R"):
                            database_description = "KEGG REACTION"
                        else:
                            print(foreign_id)
                    if self.database_name_translation[database_description] == "ec" and (foreign_id.count(".")!=3 or "-" in foreign_id):
                        continue
                    
                    if self.database_name_translation[database_description] in database2id:
                        if type(database2id[self.database_name_translation[database_description]]) != type([]):
                            database2id[self.database_name_translation[database_description]] = [database2id[self.database_name_translation[database_description]]]
                        database2id[self.database_name_translation[database_description]].append(foreign_id)
                    else:
                        database2id[self.database_name_translation[database_description]] = foreign_id
        # remove reactome IDs as they are incorrect in GO
        if "reactome" in database2id:
            del database2id["reactome"]
        self.database.addItem(database2id,names,namespace, "go")
        self.database.addRelationship("go",go_id,relationships)
        
class NcbiTaxonomyNodesRead(Read):
    def __init__(self, database, file_content, size_known=False):
        Read.__init__(self, database, file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("taxonomy")
    def read(self, file_content):
        for line in file_content.split("\t|\n"):
            columns = line.split("\t|\t")
            if len(columns)<2:
                continue
            entry = {}
            entry["ENTRY"] = columns[0]
            if columns[1] != columns[0] and columns[1] != None:
                entry["PARENT"] = columns[1]
            self.safeProcessEntry(entry)
            self.reportReadEntry()
    def processEntry(self, entry):
        this_abstract_id = self.database.addOrganism({"taxonomy":entry["ENTRY"]},[], "taxonomy")
        if "PARENT" in entry:
            # add a relation
            # get parent abstract id
            parent_abstract_id = self.database.__getAbstractIdByIdDict__({"taxonomy":entry["PARENT"]})
            if parent_abstract_id == None:
                # parent is not in the database yet
                # create it
                parent_abstract_id = self.database.addOrganism({"taxonomy":entry["PARENT"]},[],"taxonomy")
            self.database.__addRelationshipByAbstractId__("taxonomy", "is_a",
                                                          this_abstract_id, parent_abstract_id)
class NcbiTaxonomyNamesRead(Read):
    def __init__(self, database, file_content, size_known=False):
        Read.__init__(self, database, file_content, size_known=size_known)
    def read(self, file_content):
        last_tax_id = "1"
        name_list = []
        for line in file_content.split("\t|\n"):
            if len(line) == 0:
                continue
            [tax_id,name,foo,meaning] = line.split("\t|\t")
            if tax_id == last_tax_id:
                if name not in name_list:
                    if meaning.startswith("scientific"):
                        name_list.insert(0,name)
                    else:
                        name_list.append(name)
            else:
                self.safeProcessEntry({"NAMES":name_list,"DBLINKS":{"taxonomy":last_tax_id}})
                self.reportReadEntry()
                name_list = []
                last_tax_id = tax_id
                if name not in name_list:
                    if meaning.startswith("scientific"):
                        name_list.insert(0,name)
                    else:
                        name_list.append(name)
        self.safeProcessEntry({"NAMES":name_list,"DBLINKS":{"taxonomy":last_tax_id}})
        self.reportReadEntry()
    def processEntry(self, entry):
        self.database.addOrganism(entry["DBLINKS"],entry["NAMES"], "taxonomy")


class SboRead(OboRead):
    def __init__(self, database, file_content, size_known=False):
        OboRead.__init__(self, database, file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("sbo")
    def processEntry(self,entry):
        if "id" not in entry:
            return
        sbo_id = entry["id"][0]
        dblinks = {"sbo":sbo_id}
        names = entry["name"]
        self.database.addSBOTerm(dblinks,names, "sbo")
        if "is_a" in entry:
            for line in entry["is_a"]:
                parent_sbo = line.split("!")[0].strip()
                self.database.addRelationship("sbo",sbo_id,[["is_a",parent_sbo]])


class EntrezGene2RefSeqRead(TabLineRead):
    def __init__(self, database, file_content, size_known=False):
        TabLineRead.__init__(self, database, file_content, size_known)
        self.database.addToListOfImportedDatabases("entrez_gene")
    def processEntry(self, entry):
        [tax_id, GeneID, status, RNA_nucleotide_accession, RNA_nucleotide_gi, protein_accession, protein_gi, genomic_nucleotide_accession, genomic_nucleotide_gi, start, end, orientation, assembly] = entry
        entrez_gene_aid = self.database.addSpecies({"entrez_gene":GeneID},[], "entrez_gene")
        protein_aid = self.database.addSpecies({"refseq":protein_accession.split(".")[0]},[], "entrez_gene")
        genomic_nucleotide_aid = self.database.addSpecies({"refseq":genomic_nucleotide_accession},[], "entrez_gene")
        tax_aid = self.database.addOrganism({"taxonomy":tax_id},[], "entrez")
        self.database.__addRelationshipByAbstractId__("entrez_gene", "encodes",    entrez_gene_aid, protein_aid)
        self.database.__addRelationshipByAbstractId__("entrez_gene", "isPartOf",   entrez_gene_aid, genomic_nucleotide_aid)
        self.database.__addRelationshipByAbstractId__("entrez_gene", "inOrganism", entrez_gene_aid, tax_aid)

class EntrezGene2GoRead(TabLineRead):
    def __init__(self, database, file_content, size_known=False):
        TabLineRead.__init__(self, database, file_content, size_known)
    def processEntry(self, entry):
        [tax_id, GeneID, go_term, foo1, foo2, foo3, foo4, go_type] = entry
        entrez_gene_aid = self.database.addSpecies({"entrez_gene":GeneID},[], "entrez_gene")
        go_aid = self.database.addItem({"taxonomy":tax_id},[], go_type, "entrez_gene")
        go_type2relation = {"Process":"inProcess","Function":"hasFunction","Component":"isLocated"}
        self.database.__addRelationshipByAbstractId__("entrez_gene", go_type2relation[go_type], entrez_gene_aid, go_aid)

class EntrezGeneRefseqUniprotRead(TabLineRead):
    def __init__(self, database, file_content, size_known=False):
        TabLineRead.__init__(self, database, file_content, size_known)
    def processEntry(self, entry):
        [refseq_id, uniprot_id] = entry
        aid = self.database.__getAbstractIdByIdDict__({"refseq":refseq_id})
        if aid == None:
            # refseq id not known
            return
            #self.database.addSpecies({"refseq":refseq_id,"uniprot":uniprot_id}, [], "entrez_gene")
        else:
            self.database.external_item_uniprot.insertRow([uniprot_id, aid, False])
        
class EntrezGeneInfoRead(TabLineRead):
    def __init__(self, database, file_content, size_known=False):
        TabLineRead.__init__(self, database, file_content, size_known)
    def processEntry(self, entry):
        [tax_id, GeneID, Symbol, LocusTag, Synonyms, dbXrefs, chromosome, map_location, description, type_of_gene, Symbol_from_nomenclature_authority, Full_name_from_nomenclature_authority, Nomenclature_status, Other_designations, Modification_date] = entry
        aid = self.database.__getAbstractIdByIdDict__({"entrez_gene":GeneID})
        if aid == None:
            aid = self.database.addSpecies({"entrez_gene":GeneID}, [], "entrez_gene")
        name_list = []
        if Symbol != "-" and Symbol != "NEWENTRY":
            name_list.append(Symbol)
        if Full_name_from_nomenclature_authority != "-":
            name_list.append(Full_name_from_nomenclature_authority)
        if Synonyms != "-":
            name_list = name_list+Synonyms.split("|")
        if description != "-" and len(name_list) == 0:
            name_list.append(description)
        dblinks = {"entrez_gene":GeneID}
        if dbXrefs != "-":
            for ref in dbXrefs.split("|"):
                [db, id] = ref.split(":")
                db = db.lower()
                for char in "/-":
                    db = db.replace(char,"_")
                dblinks[db] = id
        first_name=True
        for name in name_list:
            self.database.string_repr.insertRow([name, aid, first_name, name.lower()])
            first_name = False
        for db,id in dblinks.items():
            if db == "uniprotkb_swiss_prot":
                db = uniprot
            database_table = self.database.__dict__["external_item_"+db]
            if not database_table.entryInDatabaseFaster({"data_id":id,
                                                         "abstract_item":aid,
                                                         "is_reverse":False}):
                database_table.insertRow([id, aid, False])
        
class SwissProtFileRead(Read):
    def __init__(self, database, file_content, size_known=False):
        self.db2name = {"EMBL":"nsd",
                        "HGNC":"hgnc",
                        "GO":"go",
                        "InterPro":"interpro",
                        "Pfam":"pfam",
                        "PROSITE":"prosite",
                        "RefSeq":"refseq",
                        "IPI":"ipi"}
        Read.__init__(self, database, file_content, size_known=size_known)
    def read(self, file_content):
        entry = {}
        for line in file_content.split("\n"):
            if line.startswith("//"):
                self.safeProcessEntry(entry)
                self.reportReadEntry()
                entry = {}
            else:
                key = line[:2].strip()
                if len(key) != 2:
                    continue
                rest_line = line[2:].strip()
                if key not in entry:
                    entry[key] = []
                entry[key].append(rest_line)
        self.safeProcessEntry(entry)
        self.reportReadEntry()

class UniprotRead(SwissProtFileRead):
    def __init__(self, database, file_content, size_known=False):
        self.forgotten_names = []
        SwissProtFileRead.__init__(self, database, file_content, size_known=size_known)
        self.database.addToListOfImportedDatabases("uniprot")
    def processEntry(self, entry):
        if "ID" not in entry:
            print(entry)
            return
        uniprot_id = entry["ID"][0].split()[0]
        ac = entry["AC"][0].split(";")[0]
        tax_id = entry["OX"][0].split("=")[1].split(";")[0]

        names = []
        ec_number = None
        if "DE" in entry:
            for line in entry["DE"]:
                if "Contains:" in line:
                    break
                if "=" in line:
                    name = line.split("=")[1].strip().strip(";")
                    if line.split("=")[0].endswith("EC"):
                        ec_number = name
                    else:
                        names.append(name)
        names.append(uniprot_id)
            
        dblinks = {}
        protein_links = {"uniprot":ac}
        gene_links = {}
        if "DR" in entry:
            for link in entry["DR"]:
                dblinkline = link.strip(".").split("; ")
                database = dblinkline[0]
                if database in self.db2name:
                    database = self.db2name[database]
                else:
                    if database not in self.forgotten_names:
                        self.forgotten_names.append(database)
                    continue
                dbid = dblinkline[1]
                if database in dblinks:
                    if type(dblinks[database]) != type([]):
                        dblinks[database] = [dblinks[database]]
                    dblinks[database].append(dbid)
                else:
                    dblinks[database] = dbid
            
                if database == "nsd":
                    embl_protein_id = dblinkline[2]
                    protein_links[database] = embl_protein_id
                    del dblinks[database]
                    
        # move the other protein links out
        for protein_db in ["ipi"]:
            if protein_db in dblinks:
                protein_links[protein_db] = dblinks[protein_db]
                del dblinks[protein_db]
        # move the go link out
        go_id = None
        if "go" in dblinks:
            go_id = dblinks["go"]
            del dblinks["go"]
            if type(go_id) != type([]):
                go_id = [go_id]

        # move dblinks to protein or gene links
        for db,id in dblinks.items():
            if db in 'hgnc refseq'.split():
                gene_links[db] = id
            elif db in 'pfam interpro prosite'.split():
                protein_links[db] = id
            else:
                input("db class unknown "+db)
                
        protein_abstract_id = self.database.addSpecies(protein_links, names, "uniprot")
        if ec_number is not None:
            ec_abstract_id = self.database.__getAbstractIdByIdDict__({"ec":ec_number})
            if ec_abstract_id is None:
                ec_abstract_id = self.database.addSpecies({"ec":ec_number}, [], "uniprot")
            self.database.__addRelationshipByAbstractId__("uniprot",
                                                          "hasProcess",
                                                          protein_abstract_id,
                                                          ec_abstract_id)
        gene_abstract_id = None
        if len(list(gene_links.keys()))>0:
            gene_abstract_id = self.database.addGene(gene_links,[], "uniprot")
            self.database.__addRelationshipByAbstractId__("uniprot",
                                                          "encodes",
                                                          gene_abstract_id,
                                                          protein_abstract_id)
            
        if go_id:
            for gid in go_id:
                go_abstract_id = self.database.__getAbstractIdByIdDict__({"go":gid})
                if not go_abstract_id:
                    continue
                self.database.__addRelationshipByAbstractId__("uniprot",
                                                              "hasProcess",
                                                              protein_abstract_id,
                                                              go_abstract_id)
                if gene_abstract_id != None:
                    self.database.__addRelationshipByAbstractId__("uniprot",
                                                                  "hasProcess",
                                                                  gene_abstract_id,
                                                                  go_abstract_id)
        # adding the taxonomy id
        tax_abstract_id = self.database.__getAbstractIdByIdDict__({"taxonomy":tax_id})
        if not tax_abstract_id:
            tax_abstract_id = self.database.addOrganism({"taxonomy":tax_id},[], "uniprot")
        self.database.__addRelationshipByAbstractId__("uniprot",
                                                      "inOrganism",
                                                      protein_abstract_id,
                                                      tax_abstract_id)
        if gene_abstract_id != None:
            self.database.__addRelationshipByAbstractId__("uniprot",
                                                          "inOrganism",
                                                          protein_abstract_id,
                                                          tax_abstract_id)
            
        


class SdfRead(Read):
    def __init__(self, database, file_content, size_known=False):
        Read.__init__(self, database, file_content, size_known=size_known)
    def read(self, file_content):
        entry = {}
        key = None
        for line in file_content.split("\n"):
            if line.startswith("$$$$"):
                self.safeProcessEntry(entry)
                self.reportReadEntry()
                entry = {}
                key = None
            elif line.startswith(">"):
                key = line.split(">")[1].split("<")[1]
            else:
                if not key:
                    continue
                if key not in entry:
                    entry[key] = []
                entry[key].append(line)
        self.safeProcessEntry(entry)
        self.reportReadEntry()

class LargeSdfRead(Read):
    def __init__(self, database, file_content, size_known=False):
        Read.__init__(self, database, file_content, size_known=size_known)
    def read(self, file_content):
        entry = {}
        key = None
        for line in file_content.newline():
            if line.startswith("$$$$"):
                self.safeProcessEntry(entry)
                self.reportReadEntry()
                entry = {}
                key = None
            elif line.startswith(">") and "<" in line:
                key = line.split(">")[1].split("<")[1]
            else:
                if not key:
                    continue
                if key not in entry:
                    entry[key] = []
                entry[key].append(line)
        self.safeProcessEntry(entry)
        self.reportReadEntry()

class PubchemRead(Read):
    def __init__(self, database, file_content, size_known=False):
        Read.__init__(self, database, file_content, size_known=size_known)
    def read(self, file_content):
        entry = {}
        oldid = 0
        for line in file_content.split("\n"):
            if "\t" not in line:
                continue
            newid = line.split("\t")[0]
            name = "\t".join(line.split("\t")[1:])
            if newid != oldid:
                self.safeProcessEntry(entry)
                self.reportReadEntry()
                oldid = newid
                entry = {}
                entry["id"] = newid
                entry["names"] = []
            if not name.startswith("CID"):
                if name.endswith("]"):
                    name = name.split("[")[0].strip()
                for namepart in name.split("\t"):
                    entry["names"].append(namepart)
        self.safeProcessEntry(entry)
        self.reportReadEntry()


class PubchemCompoundRead(PubchemRead):
    def __init__(self, database, file_content, size_known=False):
        ignore_regex = ["LS-\\d+", # leadscope
                        "[A-Z]?\\d+\\_ALDRICH", # sigma aldrich
                        "I\\d{2}-\\d+", # IS chemical technology
                        "CCRIS \\d+",
                        "HSDB \\d+",
                        "WLN: [A-Z0-9]+",
                        "[A-Z]{2}\\d-\\d+",
                        "BRN \\d+",
                        "EINECS \\d+-\\d+-\\d+",
                        "NSC.\\d+",
                        "Caswell No.*",
                        "WLN: [A-Z]+ [A-Z]+ [A-Z]+",
                        "MLS\\d{9}",
                        "AIDS[-]?\\d{6}",
                        "ALD-[A-Z]?\\d{6}",
                        "DRG-\\d{4}",
                        "STK\\d{6}",
                        "ZINC\\d{8}",
                        "NCGC\\d+-\\d+",
                        "SMR\\d{9}",
                        "EPA Pesticide Chemical Code .*",
                        "[A-Z]?\\d{2,8}\\_(SIAL|RIEDEL|FLUKA|SUPELCO)",
                        "ST\\d{7}",
                        "\\d+-\\d+-\\d+-\\d+.*", # beilstein
                        "C\\d{6}", # comparative toxigenomics database
                        "InChI=.*"]
        self.database2regex = {"kegg_compound":re.compile("^C\\d{5}$"),
                               "chebi":re.compile("^CHEBI:\\d+$"),
                               "cas":re.compile("^\\d+-\\d+-\\d+$"),
                               "biocyc":re.compile("^CPD-\\d+$"),
                               "ignore":re.compile("^(" + "|".join(ignore_regex) + ")$"),
                               "um_bbd":re.compile("^c\\d{4,}$")}
        PubchemRead.__init__(self, database, file_content, size_known=size_known)
    def processEntry(self, entry):
        if len(list(entry.keys())) == 0:
            return
        dblinks = {"pubchem_compound":entry["id"]}
        namelist = []
        for name in entry["names"]:
            matches_db = False
            for db,regex in self.database2regex.items():
                if regex.match(name):
                    matches_db = True
                    if db == "ignore":
                        continue
                    if db not in dblinks:
                        dblinks[db] = name
                    else:
                        if type(dblinks[db]) != type([]):
                            dblinks[db] = [dblinks[db]]
                        dblinks[db].append(name)
            if not matches_db:
                namelist.append(name)
        self.database.addSpecies(dblinks,namelist, "pubchem_compound")


class PubchemCompound2SubstanceRead(PubchemRead):
    def __init__(self, database, file_content, size_known=False):
        PubchemRead.__init__(self, database, file_content, size_known=size_known)
    def processEntry(self, entry):
        if len(list(entry.keys())) == 0:
            return
        compound_id = entry["id"]
        caid = self.database.__getAbstractIdByIdDict__({"pubchem_compound":compound_id})
        if caid == None:
            raise Hell()
        substance_links = [entry["names"][index] for index in range(len(entry["names"])) if index % 2 == 0]
        for substance_id in substance_links:
            said = self.database.addSpecies({"pubchem_substance":substance_id},[], "pubchem_compound")
            self.database.__addRelationshipByAbstractId__("pubchem_compound", "is_a", said, caid)

class PubchemCompoundParentRead(PubchemRead):
    def __init__(self, database, file_content, size_known=False):
        PubchemRead.__init__(self, database, file_content, size_known=size_known)
    def processEntry(self, entry):
        if len(list(entry.keys())) == 0:
            return
        child = entry["id"]
        parents = [p for p in entry["names"] if p != child]
        if len(parents) == 0:
            return
        caid = self.database.__getAbstractIdByIdDict__({"pubchem_compound":child})
        if caid == None:
            raise Hell()
        for parent in parents:
            paid = self.database.__getAbstractIdByIdDict__({"pubchem_compound":parent})
            if paid == None:
                raise Hell()
            self.database.__addRelationshipByAbstractId__("pubchem_compound", "is_a", caid, paid)

class InterproRead(Read):

    db2internal = {'swissprot':'uniprot',
                   'trembl':'uniprot',
                   '':''}
    ignoredb = 'profile smart gene3d panther prints ssf prositedoc prodom tigrfams hamap cazy'.split()

    class InterproHandler(xml.sax.ContentHandler):
        entry = {}
        element_stack = []
        relation_blacklist = ['structure_db_links','publication','external_doc_list','li','p','pre','b']
        def startElement(self, name, attrs):
            self.element_stack.append(str(name))
            if name == 'interpro':
                self.entry = {}
                self.entry['id'] = str(attrs['id'])
                self.entry['short_name'] = str(attrs['short_name'])
                self.entry['relations'] = []
                self.entry['names'] = []
            elif name == 'db_xref':
                db = str(attrs['db'])
                id = str(attrs['dbkey'])
                rel_type = self.element_stack[-2]
                if rel_type not in self.relation_blacklist:
                    self.entry['relations'].append( [db,id,rel_type] )
            elif name == 'rel_ref':
                db = str('interpro')
                id = str(attrs['ipr_ref'])
                rel_type = self.element_stack[-2]
                if rel_type not in self.relation_blacklist:
                    self.entry['relations'].append( [db,id,rel_type] )

        def endElement(self, name):
            if name == "interpro":
                self.processEntry(self.entry)
            del self.element_stack[-1]

        def characters(self, data):
            if self.element_stack[-1] == "name":
                text = data.strip()
                if text != "":
                    self.entry['names'].append( str(text) )
                    
    def processEntry(self, entry):
        self.reportReadEntry()
        id = entry['id']
        names = entry['names']+[entry['short_name']]
        ipaid = self.database.addSpecies({"interpro":id}, names, "interpro")
        for [db, dbid, relation] in entry['relations']:
            db = db.lower()
            if db in self.ignoredb: continue
            if db in self.db2internal:
                db = self.db2internal[db]
            other_aid = self.database.__getAbstractIdByIdDict__({db:dbid})
            if other_aid == None:
                other_aid = self.database.addSpecies({db:dbid}, [], "interpro")
            if relation == "parent_list":
                self.database.__addRelationshipByAbstractId__("interpro", "is_a", ipaid, other_aid)
            elif relation == "contains":
                self.database.__addRelationshipByAbstractId__("interpro", "has_part", other_aid, ipaid)
            elif relation == "abstract":
                pass
            elif relation == "member_list":
                self.database.__addRelationshipByAbstractId__("interpro", "is_a", ipaid, other_aid)
            elif relation == "child_list":
                self.database.__addRelationshipByAbstractId__("interpro", "is_a", other_aid, ipaid)
            elif relation == "example":
                self.database.__addRelationshipByAbstractId__("interpro", "is_a", other_aid, ipaid)
            elif relation == "found_in":
                self.database.__addRelationshipByAbstractId__("interpro", "has_part", ipaid, other_aid)
            else:
                print("Unknown relation",relation)
                input("?")

    def read(self, file_content):
        # this code is getting quite long, because we need to shut down the dtd validation
        # creating our handler
        _handler = self.InterproHandler()
        _handler.processEntry = self.processEntry
        # fast import of String IO
        try:
            from io import StringIO
        except ImportError:
            from io import StringIO
        # creating the parser
        xml_parser = xml.sax.make_parser()
        # setting our own handler class
        xml_parser.setContentHandler(_handler)
        # shut off dtd validation !!!
        xml_parser.setFeature(xml.sax.handler.feature_external_ges, 0)
        # preparing the input
        input_source = xml.sax.InputSource()
        input_source.setByteStream(StringIO(file_content))
        # go!
        xml_parser.parse(input_source)
            

class IntEnzRead(SwissProtFileRead):
    def __init__(self,database,file_content, size_known=False):
        SwissProtFileRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        if not 'ID' in entry:
            return
        ec_number = entry['ID'][0]
        names = []
        for key in 'DE AN CA'.split():
            if key in entry:
                for name in entry[key]:
                    names.append(name)
        prosite_links = []
        if 'PR' in entry:
            prosite_links = [x.split(';')[1].strip() for x in entry['PR']]
        uniprot_links = []
        if 'DR' in entry:
            for line in entry['DR']:
                for ue in line.strip(';').split(';'):
                    [uid, ugene] = [x.strip() for x in ue.split(',')]
                    if ugene.endswith('_HUMAN'):
                        uniprot_links.append(uid)
        ec_aid = self.database.__getAbstractIdByIdDict__({'ec':ec_number})
        if not ec_aid:
            ec_aid = self.database.addSpecies({'ec':ec_number},names,'intenz')
        # process UniProt and prosite IDs
        for [db_name, idlist, relation_type] in [['uniprot', uniprot_links, 'is_a'], ['prosite', prosite_links, 'is_a']]:
            for _id in idlist:
                db2id = {db_name:_id}
                uaid = self.database.__getAbstractIdByIdDict__(db2id)
                if not uaid:
                    uaid = self.database.addSpecies(db2id, [], 'intenz')
                self.database.__addRelationshipByAbstractId__('intenz', relation_type, uaid, ec_aid)

class BioGridRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        self.ready = False
        self.source2name = {'ENTREZ_GENE': 'entrez_gene',
                            'ENTREZ_GENE_ETG': None,
                            'DICTYBASE': None,
                            'XENBASE': 'xenbase',
                            'GENBANK_DNA_ACCESSION': None,
                            'GENBANK_PROTEIN_ACCESSION': None,
                            'GENBANK_GENOMIC_DNA_GI': None,
                            'GENBANK_GENOMIC_DNA_ACCESSION': None,
                            'GENBANK_PROTEIN_GI': None,
                            'GENBANK_DNA_GI': None,
                            'REFSEQ_PROTEIN_ACCESSION': 'refseq',
                            'REFSEQ_DNA_ACCESSION': None,
                            'REFSEQ_GENOMIC_DNA_ACCESSION': None,
                            'TREMBL': None,
                            'ECOGENE': 'ecogene',
                            'RGD': 'rgd',
                            'ENSEMBL': 'ensembl',
                            'ECOCYC': 'ecocyc',
                            'MIM': 'mim',
                            'HPRD': 'hprd',
                            'FLYBASE': 'flybase',
                            'CGNC': 'cgnc',
                            'GRID_LEGACY': None,
                            'VECTORBASE': 'vectorbase',
                            'SGD': 'sgd',
                            'IMGT/GENE-DB': 'imgt_gene_db',
                            'ZFIN': 'zfin',
                            'RATMAP': 'ratmap',
                            'WORMBASE': 'wormbase',
                            'SWISSPROT': 'uniprot',
                            'MAIZEGDB': 'maizegdb',
                            'MGI': 'mgi',
                            'TAIR': 'tair',
                            'ANIMALQTLDB': None,
                            'HGNC': 'hgnc',
                            'SYNONYM': None,
                            'OFFICIAL_SYMBOL': None,
                            'SYSTEMATIC_NAME': None,
                            'BIOGRID': None}
        self.source2relation = {'REFSEQ_PROTEIN_ACCESSION':'is_a',
                                'SWISSPROT':'is_a',
                                'ENSEMBL':'is_a',
                                'MIM':'is_a',
                                'HPRD':'is_a'}
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        if entry[0] == 'BIOGRID_ID':
            self.ready = True
            return
        if not self.ready:
            return
        [biogrid_id, identifier, source, organism] = entry
        aid = self.database.__getAbstractIdByIdDict__({'biogrid':biogrid_id})
        if aid is None:
            aid = self.database.addSpecies({'biogrid':biogrid_id}, [], 'biogrid')
        mysource = self.source2name[source]
        if mysource is None or mysource == '':
            return
        if source not in self.source2relation:
            self.database.addExternalItem(aid, {mysource:identifier})
        else:
            relation = self.source2relation[source]
            raid = self.database.__getAbstractIdByIdDict__({mysource:identifier})
            if raid is None:
                raid = self.database.addSpecies({mysource:identifier}, [], 'biogrid')
            self.database.__addRelationshipByAbstractId__('biogrid', relation, raid, aid)


class DruggableRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        self.druggable_aid = database.addSpecies({"druggable":"1"},["druggable domain"],"druggable")
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        [interpro_domain_id] = entry
        interpro_aid = self.database.addSpecies({"interpro":interpro_domain_id},[],"druggable")
        self.database.__addRelationshipByAbstractId__("interpro", "is_a", interpro_aid, self.druggable_aid)
        

class MatadorRead(TabLineRead):
    def __init__(self,database,file_content, size_known=False):
        self.druggable_aid = database.addSpecies({"druggable":"2"},["druggable"],"matador")
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        [pubchem_compound_id, chemical_name, atc, protein_id, protein_name, mesh_id, uniprot_id, protein_score, protein_annotation, mesh_score, mesh_annotation, matador_score, matador_annotation] = entry
        if pubchem_compound_id == 'chemical_id':
            return
        # building the resource_id_list
        resource_id_list = []
        if protein_id.startswith('9606.ENS'):
            resource_id_list.append( ["ensembl", protein_id[5:]] )
        elif protein_id.startswith('D') and len(protein_id)==7:
            # MeSH ID terms
            pass
        for uid in uniprot_id.split(): # add uniprot ids without the "_HUMAN" tag behind the ID
            resource_id_list.append( ["uniprot", uid] )
            if "_" in uid:
                resource_id_list.append( ["uniprot", uid.split("_")[0]] )
        # adding the actual species
        found_species = False
        for [resource, rid] in resource_id_list:
            entry_aid = self.database.__getAbstractIdByIdDict__({resource:rid})
            if entry_aid != None:
                self.database.__addRelationshipByAbstractId__("matador", "is_a", entry_aid, self.druggable_aid)
                found_species = True
                break
        # if nothing has been found just add the first reference
        if not found_species and resource_id_list:
            [resource, rid] = resource_id_list[0]
            entry_aid = self.database.addSpecies({resource:rid},[protein_name],"matador")
            self.database.__addRelationshipByAbstractId__("matador", "is_a", entry_aid, self.druggable_aid)
        

class DrugcardRead(Read):
    def __init__(self,database,file_content, size_known=False):
        Read.__init__(self,database,file_content, size_known=size_known)
    def read(self,file_content):
        entry = {}
        current_header = None
        for line in file_content.split("\n"):
            if line.startswith("#BEGIN_DRUGCARD"):
                entry = dict(id=[line.split()[1]])
            elif line.startswith("#END_DRUGCARD"):
                if entry != {}:
                    self.safeProcessEntry(entry)
                    self.reportReadEntry()
            elif line=="":
                continue
            elif line.startswith("#"):
                current_header = line[2:-1]
                entry[current_header] = []
            else:
                entry[current_header].append(line)
        self.safeProcessEntry(entry)

class DrugbankRead(DrugcardRead):
    def __init__(self,database,file_content, size_known=False):
        self.druggable_aid = database.addSpecies({"druggable":"3"},["druggable"],"drugbank")
        self.created = 0
        self.db_translate = {'PubChem_Compound':'pubchem_compound',
                             'SwissProt':'uniprot',
                             'PubChem_Substance':'pubchem_substance',
                             'KEGG_Drug':'kegg_drug',
                             'KEGG_Compound':'kegg_compound',
                             'CAS_Registry_Number':'cas'}
        DrugcardRead.__init__(self,database,file_content, size_known=size_known)
        print("Drugbank created %s entries"%self.created)
    def processEntry(self,entry):
        entry = dict([(x, entry[x]) for x in entry if entry[x][0] != 'Not Available'])
        id = entry['id'][0]
        links = dict([(x[:-3],entry[x][0]) for x in list(entry.keys()) if x.endswith('_ID') and not x.startswith('Drug_Target')])
        new_links = dict([(self.db_translate[x], links[x]) for x in links if x in self.db_translate])
        names = []
        [[names.append(y) for y in entry[x]] for x in 'Brand_Names Generic_Name Chemical_IUPAC_Name'.split() if x in entry]
        # add drug
        drug_aid = self.database.__getAbstractIdByIdDict__(new_links)
        if drug_aid is None:
            self.created = self.created + 1
            drug_aid = self.database.addSpecies(new_links,names,"drugbank")
        # add targets
        drug_target_indices = list(set([x.split('_')[2] for x in entry if x.startswith('Drug_Target_')]))
        for dtindex in drug_target_indices:
            starter = "Drug_Target_%s_"%dtindex
            this_entry = dict([(x[len(starter):], entry[x]) for x in entry if x.startswith(starter)])
            links = dict([(x[:-3],this_entry[x][0]) for x in list(this_entry.keys()) if x.endswith('_ID')])
            if "SwissProt" not in links:
                continue
            links = {'uniprot':links['SwissProt']}
            entry_aid = self.database.__getAbstractIdByIdDict__(links)
            if entry_aid is None:
                self.created = self.created + 1
                entry_aid = self.database.addSpecies(links,this_entry['Name'],"drugbank")
            self.database.__addRelationshipByAbstractId__("drugbank", "is_a", entry_aid, self.druggable_aid)
            self.database.__addRelationshipByAbstractId__("drugbank", "is_druggable_by", entry_aid, drug_aid)


class StitchActionsRead(TabLineRead):
    def __init__(self,database, file_content, size_known=False):
        global tempvar
        tempvar = []
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        global tempvar
        if entry[2] == 'inhibition':
            if entry[4] == '0':
                # protein, chemical
                protein = entry[0]
                chemical = entry[1]
            else:
                chemical = entry[0]
                protein = entry[1]
            #protein = '.'.join(protein.split('.')[1:])
            tempvar.append([protein, chemical])
class StringProteinAliasesRead(TabLineRead):
    def __init__(self,database, file_content, size_known=False):
        self.ensembl_re = re.compile('^ENS[A-Z]*[FPTG]\\d{11}')
        self.proteins = {}
        for x in tempvar:
            self.proteins[x[0]] = True
        self.druggable_aid = database.addSpecies({"druggable":"4"},["druggable"],"stitch")
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
    def processEntry(self,entry):
        if entry[3].endswith('AC'):# we ignore the IDs at this point or entry[3].endswith('ID'):
            protein_id = entry[0]+'.'+entry[1]
            if protein_id in self.proteins:
                uniprot_id = entry[2]
                if self.ensembl_re.match(entry[1]):
                    source = 'ensembl'
                else:
                    source = 'locus'
                entry_aid = self.database.__getAbstractIdByIdDict__({"stitch":protein_id})
                if entry_aid is None:
                    entry_aid = self.database.addSpecies({"stitch":protein_id, source:entry[1]},[protein_id],"stitch")
                    self.database.__addRelationshipByAbstractId__("stitch", "is_a", entry_aid, self.druggable_aid)
                # add uniprot id and link
                uniprot_aid = self.database.__getAbstractIdByIdDict__({"uniprot":uniprot_id})
                if uniprot_aid is None:
                    uniprot_aid = self.database.addSpecies({"uniprot":uniprot_id},[],"stitch")
                self.database.__addRelationshipByAbstractId__("stitch", "is_encoded_by", uniprot_aid, entry_aid)
class StitchChemicalSourcesRead(TabLineRead):
    def __init__(self,database, file_content, size_known=False):
        self.source_conversion = {"PC"   :"pubchem_compound",
                                  "KEGG" :"kegg_compound",
                                  "ChEBI":"chebi"}
        self.chemicals = {}
        for x in tempvar:
            self.chemicals[x[1]] = True
        TabLineRead.__init__(self,database,file_content, size_known=size_known)
        self.add_inhibitions()
    def processEntry(self,entry):
        if entry[0] in self.chemicals:
            [cid, source, sid] = entry
            if source not in self.source_conversion:
                return
            caid = self.database.addSpecies({"stitch":cid},[cid],"stitch")
            said = self.database.addSpecies({self.source_conversion[source]:sid},[],"stitch")
            # add relation between stitch ID and the other id
            if said not in self.database.__getRelatedAbstractIdsByAbstractId__(caid):
                self.database.__addRelationshipByAbstractId__("stitch", "is_version_of", said, caid)
    def add_inhibitions(self):
        for [protein, chemical] in tempvar:
            entry_aid = self.database.__getAbstractIdByIdDict__({"stitch":protein})
            drug_aid  = self.database.__getAbstractIdByIdDict__({"stitch":chemical})
            if entry_aid is None or drug_aid is None:
                continue
            self.database.__addRelationshipByAbstractId__("stitch", "is_druggable_by", entry_aid, drug_aid)

class BindingDBSDFRead(LargeSdfRead):
    def __init__(self, database, file_content, size_known=False):
        self.int_at_end = re.compile("^.* \\d+$")
        self.druggable_aid = database.addSpecies({"druggable":"5"},["druggable"],"bindingdb")
        LargeSdfRead.__init__(self, database, file_content, size_known=size_known)
    def processEntry(self, entry):
        if 'BindingDB monomerid' not in entry:
            return
        bindingdb_id = entry['BindingDB monomerid'][0]
        baid = self.database.__getAbstractIdByIdDict__({"bindingdb":bindingdb_id})
        if baid is None:
            baid = self.database.addSpecies({"bindingdb":bindingdb_id}, [], "bindingdb")
        if 'BindingDB Display Name' in entry:
            bindingdb_name = entry['BindingDB Display Name'][0]
            self.database.addName(baid, [bindingdb_name])
        target_numbers = set([x.split()[-1] for x in entry if ' ' in x and self.int_at_end.match(x)])
        for target_number in target_numbers:
            target_name = entry['TARGET Biomolecule '+target_number][0]
            uniprot_id = entry['UniProtKB Accession Number '+target_number][0]
            target_aid = self.database.__getAbstractIdByIdDict__({"uniprot":uniprot_id})
            if target_aid is None:
                target_aid = self.database.addSpecies({"uniprot":uniprot_id}, [target_name], "bindingdb")
            self.database.__addRelationshipByAbstractId__("bindingdb", "is_a", target_aid, self.druggable_aid)
            self.database.__addRelationshipByAbstractId__("bindingdb", "is_druggable_by", target_aid, baid)

class BindingDBTXTRead(TabLineRead):
    def __init__(self, database, file_content, size_known=False):
        TabLineRead.__init__(self, database, file_content, size_known=size_known)
    def processEntry(self, entry):
        [bindingdb_id, chebi_id] = entry
        chebi_aid = self.database.__getAbstractIdByIdDict__({"chebi":chebi_id})
        if chebi_aid is None:
            baid = self.database.addSpecies({"chebi":chebi_id, "bindingdb":bindingdb_id}, [], "bindingdb")
        else:
            #self.database.addExternalItem(chebi_aid, {"bindingdb":bindingdb_id})
            baid = self.database.addSpecies({"chebi":chebi_id, "bindingdb":bindingdb_id}, [], "bindingdb")




def readLocal(filename,database,ReadClass,size_known=False):
    """A fast test function to integrate databases from local flatfiles."""
    f = open(filename, "r")
    content = "".join(f.readlines())
    f.close()
    ReadClass(database, content, size_known)

class LargeFile:
    def __init__(self, filename):
        self.filename = filename
        self.file = open(self.filename, "r")
    def close(self):
        self.file.close()
    def newline(self):
        read = ""
        done = False
        while True:
            # last line
            if done:
                yield read
                break
            while "\n" not in read and not done:
                this_read = self.file.read(1024)
                if len(this_read) <3:
                    # last line found
                    done = True
                read += this_read
            i = read.find("\n")
            line = read[:i]
            read = read[i+1:]
            yield line
            
        

def readLargeLocal(filename,database,ReadClass,size_known=False):
    """A test function to integrate databases from local flatfiles which are too large for the normal reading."""
    
    f = LargeFile(filename)
    ReadClass(database, f, size_known)
    f.close()




if __name__=="__main__":
    from . import database
    d = database.Database()
    import datetime, os

    from optparse import OptionParser
    parser = OptionParser("%prog [options]\n\n libSBAnnotation option parser\n no -p option assumes that you have the raw database files in this diretory")
    parser.add_option("-t", "--test", dest="use_test_databases", help="also include hidden database files", action="store_true", default=False )
    parser.add_option("-b", "--big", dest="use_big_databases", help="use the big versions of the databases", action="store_true", default=False )
    parser.add_option("-d", "--druggability", dest="use_druggability", help="include druggability information", action="store_true", default=False )
    parser.add_option("-p", "--path", dest="path", help="path to the raw dabase files" )
    (options, args) = parser.parse_args()
    path= '.'
    if options.path: path = options.path
    start = datetime.datetime.now()
    print("Starting reading")

    if options.use_big_databases:
        ncbi_prepath = "/ncbi_taxonomy/"
        size = 511111
    else:
        ncbi_prepath = "/ncbi_taxonomy/small_"
        size = False


    #readLargeLocal(path+"/BIOGRID-IDENTIFIERS-3.1.74.tab.txt",d,BioGridRead)
    #d.dump()
    #d.load()
    #import sys
    #sys.exit(1)




    readLocal(path+ncbi_prepath+"nodes.dmp", d, NcbiTaxonomyNodesRead, size)
    readLocal(path+ncbi_prepath+"names.dmp", d, NcbiTaxonomyNamesRead, size)
    print("Finished reading NCBI Taxonomy")

    readLocal(path+"/gene_ontology.1_2.obo",d,GeneOntologyRead)
    print("Finished reading GeneOntology")

    readLocal(path+"/chebi.obo",d,ChebiRead)
    print("Finished reading ChEBI")

    readLocal(path+"/compound",d,KeggCompoundRead)
    print("Finished reading KEGG Compound")

    readLocal(path+"/drug",d,KeggDrugRead)
    print("Finished reading KEGG Drug")

    readLocal(path+"/enzyme",d,KeggEnzymeRead)
    print("Finished reading KEGG Enzyme")

    readLocal(path+"/reaction",d,KeggReactionRead)
    print("Finished reading KEGG Reaction")

    readLocal(path+"dbxref.tab",d,SGDCrossRefRead)
    print("Finished reading SGD X Ref")

    readLocal(path+"/interpro.xml",d,InterproRead)
    print("Finished reading Interpro")

    #readLocal(path+"telomere.nodes", d, AllysonNodesRead, size)
    #readLocal(path+"telomere.edges", d, AllysonEdgesRead, size)
    #print "Finished reading Allyson"

    if options.use_test_databases:
        #readLocal(path+"/genome", d, KeggGenomeRead)
        #print "Finished reading KEGG Genome"

        #readLocal(path+"/reactome",d,ReactomeRead)
        #print "Finished reading Reactome"

        readLocal(path+"/SBO_OBO.obo", d, SboRead)
        print("Finished reading Systems Biology Ontology")

        if options.use_big_databases:
            entrez_prepath = "/entrez_gene/"
            refseq_size = False
            go_size = False
            info_size = 5482124
            #entrez_prepath = "/entrez_gene/small_"
            #refseq_size = 670842
            #go_size = 871650
            #info_size = 347479
        else:
            entrez_prepath = "/entrez_gene/smaller_"
            refseq_size = 150210
            go_size = 140792
            info_size = 45411
        #readLocal(path+entrez_prepath+"gene2refseq", d, EntrezGene2RefSeqRead, refseq_size)
        #readLocal(path+entrez_prepath+"gene2go", d, EntrezGene2GoRead, go_size)
        #readLocal(path+"/entrez_gene/gene_refseq_uniprotkb_collab", d, EntrezGeneRefseqUniprotRead, 3930294)
        #readLocal(path+entrez_prepath+"gene_info", d, EntrezGeneInfoRead, info_size)
        #print "Finished reading Entrez Gene"

        uniprot_sizes = {"human":66509}
        organisms = "archaea bacteria fungi human invertebrates mammals plants rodents unclassified vertebrates viruses".split()
        if options.use_big_databases:
            uniprot_path = "uniprot/"
        else:
            uniprot_path = "uniprot_small/"
        for organism in uniprot_sizes:#organisms:
            size = False
            if organism in uniprot_sizes:
                size = uniprot_sizes[organism]
            readLocal(path+uniprot_path+"uniprot_trembl_"+organism+".dat", d, UniprotRead, size)
        readLocal(path+uniprot_path+"uniprot_sprot.dat", d, UniprotRead, False)
        print("Finished reading Uniprot")

        readLargeLocal(path+"/BIOGRID-IDENTIFIERS-3.1.74.tab.txt",d,BioGridRead)
        print("Finished reading BioGrid")

        #readLocal(path+"/CID-Synonym", d, PubchemCompoundRead, 42627806)
        #readLocal(path+"/CID-SID", d, PubchemCompound2SubstanceRead, 42627806)
        #readLocal(path+"/CID-Parent", d, PubchemCompoundParentRead, 42627806)
        #print "Finished reading Pubchem Compound"

    readLocal(path+"/enzyme.dat",d,IntEnzRead)
    print("Finished reading IntEnz")

    if options.use_druggability:
        readLocal(path+"/druggable",d,DruggableRead)
        print("Finished reading Druggable Domains")

        readLocal(path+"/matador.tsv",d,MatadorRead)
        print("Finished reading Matador")

        readLocal(path+"/drugcards.txt",d,DrugbankRead)
        print("Finished reading Drugbank")

        readLocal(path+"/stitch/small.actions.v2.0.tsv",d,StitchActionsRead)
        readLocal(path+"/stitch/small.protein.aliases.v8.3.txt",d,StringProteinAliasesRead)
        readLocal(path+"/stitch/small.chemical.sources.v2.0.tsv",d,StitchChemicalSourcesRead)
        print("Finished reading STITCH")

        readLocal(path+"/BDB_chebiId.txt",d,BindingDBTXTRead)
        readLargeLocal(path+"/BindingDB2D.sdf",d,BindingDBSDFRead)
        print("Finished reading BindingDB")
        
    print("Took",datetime.datetime.now()-start)

    start = datetime.datetime.now()
    print("Writing database to disk")
    d.dump()
    print("Took",datetime.datetime.now()-start)
    print("Testing the loading")
    d.load()
    print("Took",datetime.datetime.now()-start)
    

