#!/usr/bin/env python

import csv, os.path, copy

class Qualifiers:
    def __init__(self):
        self.qualifiers = {}
    def load_from_file(self, filename=None):
        if filename == None:
            import libSBAnnotation.config
            filename = libSBAnnotation.config.Config().getpath('resources','qualifiers')
        self.qualifiers = {}
        csv_file = csv.reader(open(filename), delimiter=',')
        rows = [x for x in csv_file]
        for row,row_index in zip(rows,list(range(len(rows)))):
            if row_index == 0:
                continue
            columns = [x for x in row]
            row_dict = {}
            for column, column_index in zip(columns, list(range(len(columns)))):
                if column_index == 0:
                    continue
                row_dict[rows[0][column_index]] = float(column)
            self.qualifiers[rows[0][row_index]] = row_dict
        return copy.deepcopy(self.qualifiers)

relation_type_score = {'is_a'                     : .5,
                       'part_of'                  : .1,
                       'has_part'                 : .1,
                       'regulates'                : .01,
                       'positively_regulates'     : .01,
                       'negatively_regulates'     : .01,
                       'is_tautomer_of'           : .9,
                       'is_enantiomer_of'         : .01,
                       'is_conjugate_acid_of'     : .9,
                       'is_conjugate_base_of'     : .9,
                       'has_role'                 : .75,
                       'has_functional_parent'    : 0.,
                       'is_substituent_group_from': .01,
                       'has_parent_hydride'       : .9,
                       'encodes'                  : .1, # new from here
                       'is_encoded_by'            : .1,
                       'is_version_of'            : .3,
                       'hasFunction'              : .75,
                       'hasProcess'               : .25,
                       'inOrganism'               : .1,
                       'inProcess'                : .25,
                       'isLocated'                : .1,
                       'isPartOf'                 : .25,
                       'is_druggable_by'          : 0.}


if __name__ == '__main__':
    q = Qualifiers()
    print(q.load_from_file())
