#!/usr/bin/env python
import semanticSBML.annotate
import libSBAnnotation.annotate
import libSBAnnotation.database_new
database_new = libSBAnnotation.database_new.Meta_DB()

class Cofactors(object):
    _borg = {}

    def __init__(self):
        self.__dict__ = self._borg
        if 'cofactor_aids' not in self.__dict__:
            self.cofactor_aids = self.compute_cofactor_aids()

    def compute_cofactor_aids(self):
        aids = []
        db2id = {'chebi':['CHEBI:%s'%x for x in '29108 33699 15996 24875 17552 15422 16761 17489 36080 18367 16908 16027 15846 18009 24636 18361'.split()],
                 'kegg_compound':['C'+x.rjust(5,'0') for x in '1 2 3 4 5 6 7 8 9 11 17 20 23 35 44 46 76 80 144 562'.split()]}
        for db in db2id:
            for _id in db2id[db]:
                aid = database_new.database_instance.__getAbstractIdByIdDict__({db:_id})
                if aid is not None and aid not in aids:
                    aids.append(aid)
        return aids

    def get_cofactor_aids(self):
        return self.cofactor_aids

    def is_cofactor_aid(self, aid):
        return aid in self.cofactor_aids

    def is_cofactor_annotation(self, annotation):
        return self.is_cofactor_aid(database_new.database_instance.__getAbstractIdByIdDict__({annotation.db:annotation.id}))

    def is_cofactor_element_annotation(self, element_annotation):
        for annotation in element_annotation.getAnnotations():
            if self.is_cofactor_annotation(annotation):
                return True
        return False

    def is_cofactor_libsbml_element(self, libsbml_element):
        return self.is_cofactor_element_annotation(semanticSBML.annotate.ElementAnnotations(libsbml_element, True))
