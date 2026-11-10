#!/usr/bin/env python

"""Class implementing the wishes of Ron to query the libSBAnnotation in such a way, that it will return related resource entries
   having a certain similarity to the initially provided one. This will make use of the similarity measure implemented in the
   clustering."""

import libSBAnnotation.database_new, libSBAnnotation.database, libSBAnnotation.annotate, libModelSimilarity.cluster
import math

class ExtendedRelationObject:
    """Fake relation object for current programming purposes.

    @type  abstract_id: int
    @param abstract_id: abstract ID of the ontology entry to start from
    @type  database2id: dict
    @param database2id: dictionary of database name string to id string for the entry to start from
    @type  relations: list
    @param relations: For each related entry it includes a list of abstract ID, database2id dictionary, and the similarity float
    """
    def __init__(self, abstract_id, database2id=None):
        self.abstract_id = abstract_id
        self.database2id = database2id
        self.relations = []
    def addRelation(self,related_abstract_id,related_database2id,
                    relation_score):
        """Adds a related entry to the object."""
        # check whether such a row is already in there
        if not related_abstract_id in self.getRelatedAbstractIds():
            # related abstract id is unknown
            self.relations.append( [related_abstract_id,related_database2id,
                                    relation_score] )
        else:
            index = self.getRelatedAbstractIds().index(related_abstract_id)
            if self.relations[index][2] < relation_score:
                self.relations[index][2] = relation_score
    def getRelatedAbstractIds(self):
        """Returns a list of abstract IDs related to the query entry"""
        return [x[0] for x in self.relations]
    def __repr__(self):
        string  = "RelationObject:\n"
        string += "\n".join( [",".join([str(x) for x in row])
                              for row in self.relations] )
        return string
    def toURNAndScore(self):
        list_of_urns_and_scores = []
        for [aid, db2id, similarity] in self.relations:
            urns = []
            if db2id == None:
                continue
            for db,id in db2id.items():
                if db in libSBAnnotation.database.internal_to_miriam_names:
                    db = libSBAnnotation.database.internal_to_miriam_names[db]
                try:
                    urn = str(libSBAnnotation.annotate.Annotation(db, id))
                    urns.append(urn)
                except libSBAnnotation.annotate.AnnotationCreationError:
                    pass
            list_of_urns_and_scores.append([urns, similarity])
        return list_of_urns_and_scores

class ExtendedMeta_DB(object):
        '''
        this class is empty because it is a singleton class
        everything important is done in the L{Local} class
        http://en.wikipedia.org/wiki/Singleton_pattern#Python
        '''

        def __init__(self):
            self.meta_db = libSBAnnotation.database_new.Meta_DB()
            self.database_instance = self.meta_db.database_instance
            for method in dir(self.meta_db):
                if method.startswith("get") or method.startswith("search"):
                    self.__dict__[method] = eval("self.meta_db."+method)
        def __determine_height(self, aid, max_depth=20):
            parents = list(self.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, "is_a", 1, "to"))
            if len(parents) == 0 or max_depth == 0:
                return 0
            else:
                max_dist = -1
                for parent_aid in parents:
                    this_dist = self.__determine_height(parent_aid, max_depth-1)
                    if this_dist > max_dist:
                        max_dist = this_dist
                return max_dist+1
        def getRelationWithThreshold(self,db='',id='',threshold=.1, use_depth_for_similarity=False):
            """
            Returns an ExtendedRelationObject containing all weighted, related ontology entries which have a similarity larger than threshold
            to the entry described by db:id. The parameter use_depth_for_similarity determines whether the depth in the ontology is used
            to reweight the similarity. If this parameter is False, only the distance in the ontology will be used.

            @type  db: string
            @param db: MIRIAM compliant name for the web resource
            @type  id: string
            @param id: ID of an entry in the MIRIAM compliant web resource
            @type  threshold: float
            @param threshold: minimal similarity of a result in the ExtendedRelationObject
            @type  use_depth_for_similarity: bool
            @param use_depth_for_similarity: if this is True than the depth of an entry in the ontology is influencing the similarity
            @rtype : ExtendedRelationObject
            @return: Similar entries with the similarity
            """
            aid = self.database_instance.__getAbstractIdByIdDict__(libSBAnnotation.database.AnnotationObject({db:id}).toInternalCompliantDatabase2id())
            if aid == None:
                return None
            if use_depth_for_similarity:
                height = self.__determine_height(aid)
            else:
                height = 1.
            threshold = math.pow(threshold, height)
            aid2similarity = {aid:1.}
            finished_aids = []
            aids_queue = [aid]
            extended_relation_object = ExtendedRelationObject(aid,{db:id})
            while len(aids_queue)>0:
                # get new aid
                current_aid = aids_queue.pop(0)
                current_sim = aid2similarity[current_aid]
                relation_object = self.database_instance.__getRelationshipsByAbstractIdWithoutCaching__(current_aid)
                for [new_aid, new_db2id, rtype, reverse] in relation_object.relations:
                    # no conversion to MIRIAM compliant names here
                    # new_db2id = dict([[libSBAnnotation.database.internal_to_miriam_names[key], val] if key in libSBAnnotation.database.internal_to_miriam_names else [key, val] for key, val in new_db2id.iteritems()])
                    relation_score = libModelSimilarity.cluster.relation_type_score[rtype]
                    new_sim = relation_score*current_sim
                    if new_sim >= threshold:
                        if new_aid not in finished_aids and new_aid not in aids_queue:
                            aids_queue.append(new_aid)
                        if new_aid not in aid2similarity or aid2similarity[new_aid] < new_sim:
                            aid2similarity[new_aid] = new_sim
                        extended_relation_object.addRelation(new_aid, new_db2id, math.pow(new_sim,(1./height)))
                finished_aids.append(current_aid)
            return extended_relation_object

if __name__=='__main__':
    d = ExtendedMeta_DB()
    print(dir(d))
    for tax_id in '1 2'.split():
        print(d.getName('taxonomy',tax_id))
        print(d.getRelationWithThreshold('taxonomy',tax_id))
        #print d.getRelationWithThreshold('Taxonomy',tax_id).toURNAndScore()
