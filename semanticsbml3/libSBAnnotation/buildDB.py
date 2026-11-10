"""
Readers for the small NCBI taxonomy extract in data/ (ncbi_taxonomy_small_*.dmp).

The original module had importers for about 30 sources (KEGG, ChEBI, Gene Ontology,
InterPro, SGD, UniProt, Entrez Gene, PubChem ...) in the formats of 2011. They were
removed in the Python 3 port: the database is now built by libSBAnnotation/update.py.
"""

from . import stderr

class2size = {}


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
