#!/usr/bin/env python

import libSBAnnotation.database_new, copy
from . import stderr
d = libSBAnnotation.database_new.Meta_DB()

class DruggableAID:
    def __init__(self, aid, trace=[]):
        self.aid = aid
        self.trace = trace
    def get_related_aids(self, relation_type="all", direction="both"):
        return d.database_instance.__getRelationshipsByAbstractId__(self.aid, relation_type=relation_type, direction=direction)
    def get_related_druggable_aids(self, aid_blacklist=[], relation_type="all", direction="both"):
        relation_object = self.get_related_aids(relation_type=relation_type, direction=direction)
        daids = []
        for [raid, db2id, rtype, reverse] in relation_object.relations:
            if raid not in aid_blacklist:
                daids.append( DruggableAID(raid, self.trace+[[self.aid, rtype]]) )
        return daids

class DruggablePredictor:
    def __init__(self):
        self.druggable_domain_aid = d.database_instance.__getAbstractIdByIdDict__({"druggable":self.druggable_id})
        self.aid_black_list = [x for x in [d.database_instance.__getAbstractIdByIdDict__({"druggable":str(i)}) for i in range(10)] if x != None and x != self.druggable_domain_aid]
        self.druggable_domain_daid = DruggableAID(self.druggable_domain_aid)
        self.daids = [self.druggable_domain_daid]
        self.aids= [self.druggable_domain_aid]
        new_daids = [self.druggable_domain_daid]

        if self.relation_types == None:
            self.relation_types = [["all","both"]] * self.range

        for [rt,di] in self.relation_types:
            new_daids = self.add_radial_relationships_for(new_daids, relation_type=rt, direction=di)
    def isDruggableBy(self, aid):
        # implement function for a druggability predictor which has a non direct way of linking drugs to druggable ids
        raise NotImplementedError
    def add_radial_relationships_for(self, daids, relation_type="all", direction="both", shut_up=False):
        size = len(daids)
        current_item = 0
        if not shut_up:
            stderr.downloadProgress(size, 0, True, False, url="Adding radials")

        new_daids = []
        for daid in daids:
            current_item = current_item + 1
            if not shut_up:
                stderr.downloadProgress(size, current_item, url="Adding radials")
            these_new_daids = daid.get_related_druggable_aids(aid_blacklist=self.aids+self.aid_black_list,
                                                              relation_type=relation_type,
                                                              direction=direction)
            # add to daids, aids, new_daids
            for this_new_daid in these_new_daids:
                # resolve to database2id
                if self.requires_resource is not None:
                    d2id = d.database_instance.getAnnotationObjectByAbstractId(this_new_daid.aid).toDatabase2id()
                    if self.requires_resource not in d2id:
                        continue
                self.daids.append(this_new_daid)
                self.aids.append(this_new_daid.aid)
                new_daids.append(this_new_daid)

        if not shut_up:
            stderr.downloadProgress(False, current_item, False, True)

        return new_daids

class DruggableDomainPredictor(DruggablePredictor):
    def __init__(self):
        self.range = 2
        self.druggable_id = "1"
        self.relation_types = [["is_a", "both"]] + [["has_part", "both"]]*6
        self.requires_resource = None
        DruggablePredictor.__init__(self)
        # without extension I get 46 out of 1642 Matador Ids
        # with 1 extension 161 and with two extensions 1411
        return
        if True:
            self.add_radial_relationships_for(copy.deepcopy(self.daids))
        else:
            for daid in self.daids:
                if daid.aid not in [self.druggable_domain_aid]:
                    try:
                        if 'uniprot' not in d.database_instance.getAnnotationObjectByAbstractId(daid.aid).toDatabase2id():
                            self.add_radial_relationships_for([daid], shut_up=True)
                    except AttributeError:
                        # something with the aid
                        pass

class MatadorPredictor(DruggablePredictor):
    def __init__(self):
        self.range = 1
        self.druggable_id = "2"
        self.relation_types = None
        self.requires_resource = 'uniprot'
        DruggablePredictor.__init__(self)

class DrugbankPredictor(DruggablePredictor):
    def __init__(self):
        self.range = 1
        self.druggable_id = "3"
        self.relation_types = None
        self.requires_resource = 'uniprot'
        DruggablePredictor.__init__(self)

class StitchPredictor(DruggablePredictor):
    def __init__(self):
        self.range = 1
        self.druggable_id = "4"
        self.relation_types = [["is_a","to"],["is_encoded_by","to"]]
        self.requires_resource = None
        DruggablePredictor.__init__(self)
    def isDruggableBy(self, aid):
        # determine type of aid
        ao = d.database_instance.getAnnotationObjectByAbstractId(aid).toDatabase2id()
        if 'uniprot' in ao:
            locus_aids = d.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, relation_type="is_encoded_by")
        elif 'ensembl' in ao or 'locus' in ao:
            # move to encoded proteins
            locus_aids = [aid]
        else:
            return []
        # move to druggable stitch compounds
        sc_aids = []
        for locus_aid in locus_aids:
            [sc_aids.append(x) for x in d.database_instance.__getRelatedAbstractIdsByAbstractId__(locus_aid, relation_type="is_druggable_by")]
        # move to real compounds
        c_aids = []
        for sc_aid in sc_aids:
            [c_aids.append(x) for x in d.database_instance.__getRelatedAbstractIdsByAbstractId__(sc_aid, relation_type="is_version_of", direction="to")]
        # remove more specific compounds and keep their parents
        children = []
        for c_aid in c_aids:
            [children.append(x) for x in d.database_instance.__getRelatedAbstractIdsByAbstractId__(c_aid, relation_type="is_a", direction="to")]
        childrens_children = []
        for c_aid in children:
            [childrens_children.append(x) for x in d.database_instance.__getRelatedAbstractIdsByAbstractId__(c_aid, relation_type="is_a", direction="to")]
        c_aids = [x for x in c_aids if x not in children and x not in childrens_children]
        return c_aids

class BindingDBPredictor(DruggablePredictor):
    def __init__(self):
        self.range = 1
        self.druggable_id = "5"
        self.relation_types = None
        self.requires_resource = 'uniprot'
        DruggablePredictor.__init__(self)


def compareMatadorDD():
    ddp = DruggableDomainPredictor()
    mp = MatadorPredictor()
    found = 0
    found3 = 0

    size = len(mp.aids)
    current_item = 0
    stderr.downloadProgress(size, 0, True, False, url="Comparing Matador2DD")

    for aid in mp.aids:
        current_item = current_item + 1
        stderr.downloadProgress(size, current_item, url="Comparing Matador2DD")
        if aid in ddp.aids:
            found = found + 1
            #daid = [x for x in ddp.daids if x.aid == aid][0]
            #if len(daid.trace)>=2:
            #    print "-"*30
            #    for [raid, relation] in daid.trace:
            #        print d.database_instance.__getNameByAbstractId__(raid), d.database_instance.getAnnotationObjectByAbstractId(raid),relation
            #    print d.database_instance.__getNameByAbstractId__(aid), d.database_instance.getAnnotationObjectByAbstractId(aid)
            #    raw_input("?")
        else:
            f = False
            for aid2 in d.database_instance.__getRelationshipsByAbstractId__(aid).getRelatedAbstractIds():
                if aid2 != mp.druggable_domain_aid:
                    if aid2 in ddp.aids:
                        found3 = found3 + 1
                        #f = True
                        #daid = [x for x in ddp.daids if x.aid == aid2][0]
                        #print "-"*30
                        #for [raid, relation] in daid.trace:
                        #    print d.database_instance.__getNameByAbstractId__(raid), d.database_instance.getAnnotationObjectByAbstractId(raid),relation
                        #print d.database_instance.__getNameByAbstractId__(aid2), d.database_instance.getAnnotationObjectByAbstractId(aid2)
                        #print d.database_instance.__getNameByAbstractId__(aid), d.database_instance.getAnnotationObjectByAbstractId(aid)
                        #raw_input("?")
                        break
            #if not f:
            #    print aid, "has not been found by the DD approach"
            #    print d.database_instance.__getNameByAbstractId__(aid), d.database_instance.getAnnotationObjectByAbstractId(aid)
            #    raw_input("?")
    stderr.downloadProgress(False, current_item, False, True)
    print("Found",found,"of",len(mp.aids),"Matador Prediction in",len(ddp.aids),"AIDs")
    print(found3, "more if extended")

class EcPredictorFake:
    name = 'Functional Relation'
    def __init__(self, aids, protein_aids):
        self.aids = aids
        self.protein_aids = protein_aids
    def isDruggableBy(self, aid):
        raise NotImplementedError

class AllDruggabilityPredictors:
    def __init__(self):
        mod = __import__(self.__module__, {}, {}, [''])
        self.predictor_names = [x for x in dir(mod) if x.endswith('Predictor') and not x.startswith('Druggable')]
        self.predictors = [eval("mod.%s()"%x) for x in self.predictor_names]
    def extend_to_ec(self):
        all_aids = []
        for predictor in self.predictors:
            for aid in predictor.aids:
                if aid not in all_aids:
                    all_aids.append(aid)
        ec_aids = []
        protein_aids = []
        # extend to ec numbers
        for aid in all_aids:
            ro = d.database_instance.__getRelationshipsByAbstractId__(aid,"hasProcess")
            for [raid, d2id, foo, bar] in ro.relations:
                if 'ec' in d2id:
                    ec_aids.append(raid)
                    protein_aids.append(aid)
        self.predictors.append(EcPredictorFake(ec_aids, protein_aids))
        self.predictor_names.append(EcPredictorFake.name)
    def isDruggableAID(self, aid):
        for predictor in self.predictors:
            if aid in predictor.aids:
                return True
        return False
    def isDruggable(self, id_dict):
        aid = d.database_instance.__getAbstractIdByIdDict__(id_dict)
        if aid is None:
            return False
        return self.isDruggableAID(aid)
    def isDruggableAnnotation(self, annotation):
        return self.isDruggable({annotation.db:annotation.id})
    def isPredictedDruggableAIDBy(self, aid):
        list_of_predictors = []
        for predictor_name, predictor in zip(self.predictor_names, self.predictors):
            if aid in predictor.aids:
                list_of_predictors.append(predictor_name[:-9])
        return list_of_predictors
    def isPredictedDruggableAnnotationBy(self, annotation):
        return self.isPredictedDruggableBy({annotation.db:annotation.id})
    def isPredictedDruggableBy(self, id_dict):
        aid = d.database_instance.__getAbstractIdByIdDict__(id_dict)
        if aid is None:
            return []
        return self.isPredictedDruggableAIDBy(aid)
    def isDruggableBy(self, aid=None, db2id=None, annotation=None):
        if aid is None:
            if db2id is None:
                if annotation is None:
                    raise Exception('No arguments given')
                db2id = {annotation.db:annotation.id}
            aid = d.database_instance.__getAbstractIdByIdDict__(db2id)
            if aid is None:
                return []
        raids = d.database_instance.__getRelationshipsByAbstractId__(aid, relation_type='is_druggable_by').getRelatedAbstractIds()
        # extend by individual functions
        for predictor in self.predictors:
            try:
                [raids.append(x) for x in predictor.isDruggableBy(aid)]
            except NotImplementedError:
                # not necessarily implemented by every druggability predictor
                pass
        return raids

class DruggableAnnotation(object):
    def __new__(cls, all_druggability_predictors, annotation):
        if not all_druggability_predictors.isDruggableAnnotation(annotation):
            return None
        return super(DruggableAnnotation, cls).__new__(cls)
    def __init__(self, all_druggability_predictors, annotation):
        self.annotation = annotation
        self.aid = d.database_instance.__getAbstractIdByIdDict__({annotation.db: annotation.id})
        self.predicted_druggable_by = all_druggability_predictors.isPredictedDruggableAIDBy(self.aid)
        self.druggable_by = all_druggability_predictors.isDruggableBy(aid=self.aid)
        self.druggable_by_aos = [d.database_instance.getAnnotationObjectByAbstractId(aid) for aid in self.druggable_by]
    def __str__(self):
        return "%s is druggable with %s as stated by %s"%(str(self.annotation), ", ".join([str(a) for a in self.druggable_by_aos]), ", ".join(self.predicted_druggable_by))
        
class DruggableElement(object):
    def __new__(cls, all_druggability_predictors, element_annotations):
        druggable_annotations = [DruggableAnnotation(all_druggability_predictors, annotation) for annotation in element_annotations.getAnnotations()]
        if len([x for x in druggable_annotations if x is not None]) == 0:
            return None
        return super(DruggableElement, cls).__new__(cls)
    def __init__(self, all_druggability_predictors, element_annotations):
        self.element_annotations = element_annotations
        self.druggable_annotations = [DruggableAnnotation(all_druggability_predictors, annotation) for annotation in element_annotations.getAnnotations()]
    def __str__(self):
        lines = [self.element_annotations.id]
        for druggable_annotation in self.druggable_annotations:
            if druggable_annotation is not None:
                lines.append("\t"+str(druggable_annotation))
        return "\n".join(lines)

class DruggableModel(object):
    def __init__(self, all_druggability_predictors, model_elements_annotations, extend_reactions=False):
        self.mea = model_elements_annotations
        if extend_reactions:
            self.extend_reactions()
        self.druggable_elements = [DruggableElement(all_druggability_predictors, element) for element in reversed(self.mea.getElementAnnotations())]
    def __str__(self):
        lines = []
        for druggable_element in self.druggable_elements:
            if druggable_element is not None:
                lines.append(str(druggable_element))
        return "\n".join(lines)
    def extend_reactions(self):
        import semanticSBML.annotate
        id2ea = {}
        for ea in self.mea.getElementAnnotations():
            id2ea[ea.id] = ea
        for ea in self.mea.getElementAnnotations():
            if ea.type == 'reaction':
                for modifier_id in [x.getSpecies() for x in ea.libsbml_element.getListOfModifiers()]:
                    for a in id2ea[modifier_id].getAnnotations():
                        annotation = semanticSBML.annotate.Annotation(db=a.db, id=a.id, qualifier_type=a.qualifier_type, qualifier=a.qualifier)
                        ea.addAnnotation(annotation)

def main():
    adp = AllDruggabilityPredictors()
    adp.extend_to_ec()
    #aid = 181744
    #print adp.isDruggableAID(aid)
    #print adp.isPredictedDruggableAIDBy(aid)
    #print d.getAnnotationObjectByAbstractId(aid)
    #print adp.isDruggable({'uniprot':'P35354'})
    #print adp.isPredictedDruggableBy({'uniprot':'P35354'})
    import semanticSBML.annotate
    import libsbml
    import os
    path = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/'
    count = 0
    for file_name in [path+x for x in os.listdir(path) if x.endswith('.xml')]:
        doc = libsbml.readSBML(file_name)
        mea = semanticSBML.annotate.ModelElementsAnnotations(doc.getModel(), suppress_errors=True)
        dm = DruggableModel(adp, mea, True)
        count += len([x for x in dm.druggable_elements if x is not None])
    print(count, "druggable elements found")
    
    

if __name__ == "__main__":
    #compareMatadorDD()
    sp = StitchPredictor()
    aid = sp.aids[-1]
    sp.isDruggableBy(aid)
