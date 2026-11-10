from . import config
"""
Classes implementing the tables of the semanticSBML database.
"""

__version__ = "1.8."


class ColumnNotFoundException(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        

class _LazyIndex(dict):
    """index of a table: column name -> {value: [row numbers]}, built when a column is first asked for
    (added in the Python 3 port; building all indexes on loading took most of the start-up time)"""
    def __init__(self, table):
        dict.__init__(self)
        self.table = table
    def __missing__(self, column):
        column_index = self.table.column_descriptions.index(column)
        index = {}
        for row_number, row in enumerate(self.table.rows):
            value = row[column_index]
            if value in index:
                index[value].append(row_number)
            else:
                index[value] = [row_number]
        self[column] = index
        return index


class BasicTable:
    """Class implementing the basic functionality of the six central tables of the database."""
    def __init__(self,column_descriptions,column_types):
        self.column_descriptions = column_descriptions
        self.column_types = column_types
        self.column_index_range = list(range(len(self.column_descriptions)))
        self.rows = []
        self._row_hashes = None
        self.row_number = 0
        self.global_index = _LazyIndex(self)
        self.config = config.Config()
    @property
    def row_hashes(self):
        if self._row_hashes is None:
            hashes = {}
            for row_number, row in enumerate(self.rows):
                hashes.setdefault(self.__hash_row__(row), []).append(row_number)
            self._row_hashes = hashes
        return self._row_hashes
    def dump(self):
        f=open(self.config.getpath("database","database_dir")+"internal_db.version_"+__version__+self.filename,"w",encoding="utf-8")
        try:
            f.write( '\n'.join(["\t".join([(y if typ=="str" else str(y)) for y,typ in zip(row,self.column_types)]) for row in self.rows]) )
        except TypeError:
            print('Error in', self.filename)
            print('TypeError in line')
            l = []
            for row in self.rows:
                l2 = []
                for y,typ in zip(row,self.column_types):
                    l2.append( y if typ=="str" else str(y) )
                try:
                    l.append("\t".join(l2))
                except TypeError:
                    print(l2)
                    print(list(zip(row,self.column_types)))
                    print(row)
                    raise
            '\n'.join(l)

        f.close()
    def load(self):
        def bool(s):
            return len(s) == 4
        # prepare entry to type dict
        entry2type_func = {}
        for index in range(len(self.column_types)):
            if self.column_types[index] != "str":
                entry2type_func[index] = eval(self.column_types[index])
        # end of conversion function to correct type
        fn=self.config.getpath("database","database_dir")+"internal_db.version_"+__version__+self.filename
        import os
        if not os.path.exists(fn): return
        f=open(fn,"r",encoding="utf-8",errors="replace")
        rows = self.rows
        conversions = list(entry2type_func.items())
        for line in f:
            entries = line.strip('\n').split('\t')
            for index, func in conversions:
                entries[index] = func(entries[index])
            rows.append(entries)
        f.close()
        self.row_number = len(rows)
        self._row_hashes = None
        self.global_index.clear()
    def __hash_row__(self,row):
        return str(row)
    def insertRow(self,row):
        # TESTS for the quality of the insertions
        #if len(row) != len(self.column_descriptions):
        #    raise Hell()
        #for x in row:
        #    if type(x) == type("") and "\n" in x:
        #        print row
        #        raise Hell()
        ## testing type of the row entry
        #for index in range(len(row)):
        #    entry = row[index]
        #    typed_entry = eval(self.column_types[index]+"(entry)")
        #    if type(entry) != type(typed_entry):
        #        print self.__class__
        #        print row
        #        print self.column_types
        #        raise Hell()
        ## END OF TESTS for quality of insertion

        if self._row_hashes is not None:
            self._row_hashes.setdefault(self.__hash_row__(row), []).append(len(self.rows))
        # adding the row
        self.rows.append(row)
        # adding it to the indexes that have been built
        for column, current_dict in self.global_index.items():
            value = row[self.column_descriptions.index(column)]
            if value not in current_dict:
                current_dict[value] = []
            current_dict[value].append(self.row_number)
        self.row_number = self.row_number + 1
    def deleteRow(self,rowIndex):
        raise NotImplementedError
    def updateRow(self,rowIndex,new_row):
        for column_index in self.column_index_range:
            if self.rows[rowIndex][column_index] != new_row[column_index]:
                current_dict = self.global_index[self.column_descriptions[column_index]]
                current_dict[self.rows[rowIndex][column_index]].remove(rowIndex)
                if new_row[column_index] not in current_dict:
                    current_dict[new_row[column_index]] = []
                current_dict[new_row[column_index]].append(rowIndex)
        if self._row_hashes is not None:
            old_hash = self.__hash_row__(self.rows[rowIndex])
            if rowIndex in self._row_hashes.get(old_hash, []):
                self._row_hashes[old_hash].remove(rowIndex)
                if not self._row_hashes[old_hash]:
                    del self._row_hashes[old_hash]
            self._row_hashes.setdefault(self.__hash_row__(new_row), []).append(rowIndex)
        self.rows[rowIndex] = new_row
    def __getRowSet__(self,column_description,value):
        if value not in self.global_index[column_description]:
            return set([])
        return set(self.global_index[column_description][value])
    def findAllRowIndices(self,description_dict):
        description_dict_keys = list(description_dict.keys())
        set_of_matching_rows = self.__getRowSet__(description_dict_keys[0],description_dict[description_dict_keys[0]])
        for key in description_dict_keys[1:]:
            value = description_dict[key]
            set_of_matching_rows = set_of_matching_rows.intersection(self.__getRowSet__(key,value))
        return list(set_of_matching_rows)
    def findAllRows(self,description_dict):
        return [self.rows[x] for x in self.findAllRowIndices(description_dict)]
    def entryInDatabase(self,description_dict):
        return len(self.findAllRowIndices(description_dict))>0
    def entryInDatabaseFaster(self,description_dict):
        row = [description_dict[desc] for desc in self.column_descriptions]
        row_hash = self.__hash_row__(row)
        if row_hash not in self.row_hashes:
            return False
        else:
            return True
    def __len__(self):
        return self.row_number
    def __str__(self):
        string = str(type(self))+"\n"
        string += ", ".join(self.column_descriptions)+"\n"
        string += "\n".join([", ".join([str(x) for x in row]) for row in self.rows])
        return string

class ExternalItem(BasicTable):
    def __init__(self,database_name):
        BasicTable.__init__(self,"data_id abstract_item is_reverse".split(),"str int bool".split())
        self.database_name = database_name
        self.filename = "external_item_"+self.database_name.replace(" ","_").replace("-","_").replace("/","_").lower()

class AbstractItem(BasicTable):
    def __init__(self):
        BasicTable.__init__(self,"type".split(),["str"])
        self.filename = "abstract_item"

class StringRepr(BasicTable):
    def __init__(self):
        BasicTable.__init__(self,"repr abstract_item is_favourite_description lower_case".split(),"str int bool str".split())
        self.filename = "string_repr"

class StringReprFromDatabase(BasicTable):
    def __init__(self):
        BasicTable.__init__(self,"repr abstract_item database lower_case".split(),"str int str str".split())
        self.filename = "string_repr_from_database"

class AbstractItemRelation(BasicTable):
    def __init__(self):
        BasicTable.__init__(self,"from_item to_item relation_type evidence".split(),"int int str str".split())
        self.filename = "abstract_item_relation"

class TakesPart(BasicTable):
    def __init__(self):
        BasicTable.__init__(self,"species_id reaction_id role stoichiometry".split(),"int int str int".split())
        self.filename = "takes_part"

class Localisation(BasicTable):
    def __init__(self):
        BasicTable.__init__(self,"abstract_item localisation general_item".split(),"int int int".split())
        self.filename = "localisation"
