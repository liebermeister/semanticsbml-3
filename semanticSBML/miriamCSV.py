#!/usr/bin/env python
import copy
import csv
import json
import io
import semanticSBML.annotate

class MiriamTable:
    columns = dict(name = "Compound",
                   urn = "MiriamAnnotations",
                   annotation_start = "MiriamID ")
    bioq = 'bio'
    qual = 'is'
    def __init__(self, list_of_lists):
        self.lol = list_of_lists
        self.header = self.lol[0] if len(self.lol) > 0 else []
        self.column_indices = dict()
        self.process_lol()
    def initialize_from_annotations(self, name_list=None, annotations=None):
        self.names = None
        if name_list is None:
            if annotations is not None:
                self.names = list(annotations.keys())
        else:
            self.names = name_list
        self.annotations = annotations
        # get the URNs from all annotations
        db2urn = {}
        for annotation_list in list(self.annotations.values()):
            for annotation in annotation_list:
                if annotation.qualifier == 'is':
                    if annotation.db not in db2urn:
                        db2urn[annotation.db] = ':'.join(str(annotation).split(':')[:-1])
        dbs = list(db2urn.keys())
        dbs.sort()
        self.lol = [[self.columns['name'], self.columns['urn']]+[self.columns['annotation_start']+db2urn[db] for db in dbs]]
        if self.names is not None:
            for name in self.names:
                urns = []
                id_list = ['']*len(dbs)
                if name in self.annotations:
                    for annotation in self.annotations[name]:
                        urns.append(["%s:%s"%(annotation.qualifier_type, annotation.qualifier), str(annotation)])
                        if annotation.qualifier == 'is':
                            db_index = dbs.index(annotation.db)
                            id_list[db_index] = annotation.id
                self.lol.append([name, json.dumps(urns)]+id_list)
    def _header_index(self, column_name):
        name = self.columns[column_name]
        if name in self.header:
            return self.header.index(name)
        else:
            return -1
    def _header_index2list(self, index):
        if index == -1:
            return None
        else:
            return [row[index] for row in self.lol[1:]]
    def process_lol(self):
        name_index = self._header_index('name')
        urn_index = self._header_index('urn')
        a = len(self.columns['annotation_start'])
        annotation_index_dict = dict([(self.header[i][a:], i) for i in range(len(self.header)) if self.header[i].startswith(self.columns['annotation_start'])])
        self.names = self._header_index2list(name_index)
        if urn_index == -1 and annotation_index_dict == {}:
            self.annotations = None
        else:
            self.annotations = []
            if urn_index != -1:
                json_urns = self._header_index2list(urn_index)
                for json_urn in json_urns:
                    _json = json.loads(json_urn)
                    this_annotations = []
                    for anno in _json:
                        if len(anno) == 1:
                            urn = anno[0]
                            bioq = self.bioq
                            qual = self.qual
                        else:
                            urn = anno[1]
                            [bioq, qual] = anno[0].split(':')
                        annotation_object = semanticSBML.annotate.Annotation(urn=urn, qualifier_type=bioq, qualifier=qual)
                        this_annotations.append(annotation_object)
                    self.annotations.append(this_annotations)
            else:
                for db, i in annotation_index_dict.items():
                    ids = self._header_index2list(i)
                    while len(ids) > len(self.annotations):
                        self.annotations.append([])
                    for j in range(len(ids)):
                        _id = ids[j]
                        if _id == "":
                            continue
                        annotation_object = semanticSBML.annotate.Annotation(urn="%s:%s"%(db, _id), qualifier_type=self.bioq, qualifier=self.qual)
                        if str(annotation_object) not in [str(x) for x in self.annotations[j]]:
                            self.annotations[j].append(annotation_object)
        
        
class MiriamCSV(object):
    delimiter = ','
    quotechar = '\"'
    def __init__(self, filename=None):
        self.lol = []
        if filename is not None:
            reader = csv.reader(open(filename,'r'), delimiter=self.delimiter, quotechar=self.quotechar)
            self.lol = [[x for x in row] for row in reader]
        self.table = MiriamTable(self.lol)
    def initialize_from_annotations(self, name_list=None, annotations=None):
        self.table.initialize_from_annotations(name_list, annotations)
    def __str__(self):
        s = io.StringIO()
        csvw = csv.writer(s, delimiter=self.delimiter, quotechar=self.quotechar)
        for row in self.table.lol:
            csvw.writerow(row)
        r = s.getvalue()
        s.close()
        return r
    def to_file(self, filename):
        f = open(filename, 'w')
        f.write(str(self))
        f.close()

class MiriamTSV(MiriamCSV):
    delimiter = '\t'

def main():
    print(MiriamCSV('miriam_test1.csv'))
    print(MiriamTSV('miriam_test2.tsv'))
    m = MiriamCSV()
    m.initialize_from_annotations(['glucose'], {'glucose':[semanticSBML.annotate.Annotation(db='KEGG Compound', id='C00080', qualifier_type='bio', qualifier='is')]})
    print(m)

if __name__ == "__main__":
    main()
