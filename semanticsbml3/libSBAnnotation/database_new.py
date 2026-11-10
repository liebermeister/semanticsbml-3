#!/usr/bin/env python

# import external modules
#import libsbml
import urllib.request, urllib.error, urllib.parse, sys, os, os.path
from . import database

# import internal modules
from . import config

##################################################
##             UNDER CONSTRUCTION               ##
##################################################

db_dir = None

class Meta_DB(object):
    '''
    this class is empty because it is a singleton class
    everything important is done in the L{Local} class
    http://en.wikipedia.org/wiki/Singleton_pattern#Python
    '''

    __instance = None

    def __new__(cls, *args, **kargs): 
        if cls.__instance is None:
            cls.__instance = Meta_DB.Local(*args, **kargs)
        return cls.__instance

    class Local:
        '''
        build the db
        downloading the databases from the internet.
        '''
        def __init__(self, suppress_build=False):
            '''
            download SBO flatfile and build database
            @param url     : url of a database flatfile
            @type url      : string
            @param filename: name of the internal database file
            @type filename : string
            '''
            #---------------------------------------------------
            #set ui functions to their defaults

            from . import stderr
            self.status_fkt = stderr.info
            self.yes_or_no_fkt = stderr.yes_or_no
            self.download_status_fkt = stderr.downloadProgress
            #---------------------------------------------------
            self._cfg = config.Config()
            self.base_dir = self._cfg.getpath("database","database_dir")
            #---------------------------------------------------
            #---------------------------------------------------
            # building the database
            preliminary_database_instance = database.Database()
            self.database_instance = None
            if not suppress_build:
                self.database_instance = preliminary_database_instance.load()
            if self.database_instance == None:
                # unable to load the database from disk.
                # The original tool now started to download and build it from
                # the internet; most of those sources no longer exist, so we stop.
                self.database_instance = preliminary_database_instance
                if not suppress_build:
                    raise DatabaseError('No annotation database found in %s. Build it with the command\n'
                                        '    semanticsbml build-database\n'
                                        '(downloads about 90 MB from EBI and geneontology.org, takes a few minutes).' % self.base_dir)
            self.__updateRelationFunctions__()
            #---------------------------------------------------
            # caches
            self.db_id_to_name = {}
            self.name_to_dblinks = {}

        def forceReload(self):
            self.database_instance = self.database_instance.load()
            self.__updateRelationFunctions__()

        def buildDB(self, *args, **kwargs):
            # (port) the original downloaded and imported about 30 sources here (KEGG, InterPro,
            # SGD ...), from addresses that no longer exist; replaced by libSBAnnotation/update.py
            raise DatabaseError('The annotation database is built with the command  semanticsbml build-database')

        def getListOfImportedDatabases(self):
            """Returns a list of imported databases. The names are identical with the names in the download sources."""
            return self.database_instance.getListOfImportedDatabases()

        def getAnnotationObject(self,urn=''):
            """
            @type  urn: string
            @param urn: MIRIAM compliant URN for an entry in a web resource
            """
            from . import annotate
            urn_converter = annotate.Annotation(urn=urn)
            return self.database_instance.getAnnotationObjectByIdDict(database.AnnotationObject({str(urn_converter.db):urn_converter.id}).toInternalCompliantDatabase2id())
        def getReactionAISetByListOfIdDicts(self,*args):
            return self.database_instance.getReactionAISetByListOfIdDicts(*args)
        def getAnnotationObjectByAbstractId(self,abstract_id):
            return self.database_instance.getAnnotationObjectByAbstractId(abstract_id)
        def getReactionDescription(self, db, id):
            aid = self.database_instance.__getAbstractIdByIdDict__({db:id})
            return self.database_instance.__getNameByAbstractReactionId__(aid)
        def searchSuperNames(self, query='', fuzziness=.5):
            """
            Convenience function doing the same as searchEntries except for the fact that exact substring search is used instead of difflib similarity.

            @param fuzziness
            
            @return: list of matching entry groups [{db1:idx,db2:idy},{db1:idz,db3:ida}]
            @rtype: [{db:id}]
            """
            annotation_object_list = []
            for super_string in self.database_instance.getSuperstringNames(query, fuzziness):
                entry_list = self.searchEntries(super_string, 1., 0)
                # just take every first annotation object
                if entry_list:
                    annotation_object_list.append(entry_list[0])
            return annotation_object_list
        def searchEntries(self, query='', fuzziness=1., depth=1):
            '''
            Seach a database entry by its name
            For fuzziness = 1. the search is exact
            For fuzziness != 1. and depth > 0 the search is extended by is_a relations
            on the basis of abstract ids
            @return: list of matching entry groups [{db1:idx,db2:idy},{db1:idz,db3:ida}]
            @rtype: [{db:id}]
            '''
            if fuzziness == 1.:
                # (the original answered from a cache of names shown earlier, which made the
                #  result of a search depend on what had been looked up before; removed)
                annotation_object_list = self.database_instance.getAnnotationObjectListByName(query)
                if annotation_object_list == None:
                    return []
                else:
                    return [annotation_object.toMiriamCompliantDatabase2id() for annotation_object in annotation_object_list if annotation_object != None]
            else:
                close_names = self.database_instance.getCloseNames(query,fuzziness)
                # make names unique
                unique_close_names = []
                [unique_close_names.append(x) for x in close_names if not x in unique_close_names]
                abstract_ids = [self.database_instance.__getAbstractIdByNameList__([name]) for name in unique_close_names]
                # if searching in depth if activated
                related_abstract_ids = []
                [[related_abstract_ids.append(x) \
                for x in list(self.database_instance.__getRelatedAbstractIdsByAbstractId__(aid, relation_type="is_a", depth=depth))] for aid in abstract_ids]
                [abstract_ids.append(x) for x in related_abstract_ids]
                # make them unique on the basis of abstract ids
                unique_abstract_ids = []
                [unique_abstract_ids.append(aid) for aid in abstract_ids if aid not in unique_abstract_ids]
                # test if we have annotation objects
                unique_abstract_ids = [x for x in unique_abstract_ids if self.database_instance.getAnnotationObjectByAbstractId(x) != None]
                return [self.database_instance.getAnnotationObjectByAbstractId(x).toMiriamCompliantDatabase2id() for x in unique_abstract_ids]
        def getRelationWithType(self,urn='',source=None,relation_type='all',direction="both"):
            """
            @type  urn: string
            @param urn: MIRIAM compliant URN for an entry in a web resource
            @type  source: string
            @param source: internal name of the source from which relations should be returned, None if all sources should be queried
            @type  relation_type: string
            @param relation_type: type of the relations which should be included into the results
            @type direction: string
            @param direction: can be both, from, or to and adds a restriction to the direction of the relations returned
            @rtype : list of lists of URNs
            @return: each list represents an ontology entry for which a list of URNs is returned
            """
            from . import annotate
            urn_converter = annotate.Annotation(urn=urn)
            aid = self.database_instance.__getAbstractIdByIdDict__(database.AnnotationObject({str(urn_converter.db):urn_converter.id}).toInternalCompliantDatabase2id())
            if aid == None:
                return []
            relation_object = self.database_instance.__getRelationshipsByAbstractIdWithoutCaching__(aid, relation_type, direction, source)
            return [[self.database_instance.getAnnotationObjectByAbstractId(aid).toMiriamCompliantDatabase2id(), rtype, rdir] for [aid, d2id, rtype, rdir] in relation_object.relations]
        def getRelationFromSource(self,urn='',source=None,relation_type='all',direction="both"):
            """
            @type  urn: string
            @param urn: MIRIAM compliant URN for an entry in a web resource
            @type  source: string
            @param source: internal name of the source from which relations should be returned, None if all sources should be queried
            @type  relation_type: string
            @param relation_type: type of the relations which should be included into the results
            @type direction: string
            @param direction: can be both, from, or to and adds a restriction to the direction of the relations returned
            @rtype : list of lists of URNs
            @return: each list represents an ontology entry for which a list of URNs is returned
            """
            from . import annotate
            urn_converter = annotate.Annotation(urn=urn)
            aid = self.database_instance.__getAbstractIdByIdDict__(database.AnnotationObject({str(urn_converter.db):urn_converter.id}).toInternalCompliantDatabase2id())
            if aid == None:
                return []
            relation_object = self.database_instance.__getRelationshipsByAbstractIdWithoutCaching__(aid, relation_type, direction, source)
            list_of_related_id_dicts = [self.database_instance.getAnnotationObjectByAbstractId(x).toDatabase2id() for x in relation_object.getRelatedAbstractIds()]
            # make the returned web resource names MIRIAM compliant
            list_of_related_id_dicts = [dict([[database.internal_to_miriam_names[key], val] if key in database.internal_to_miriam_names else [key, val] for key, val in iddict.items()]) for iddict in list_of_related_id_dicts]
            # get the URNs
            list_of_list_of_urns = []
            for iddict in list_of_related_id_dicts:
                list_of_urns = []
                for key, val in iddict.items():
                    try:
                        urn = str(annotate.Annotation(key, val))
                        list_of_urns.append(urn)
                    except annotate.AnnotationCreationError:
                        pass
                list_of_list_of_urns.append(list_of_urns)
            return list_of_list_of_urns
        def getRelation(self,db='',id='',relation_type='all',depth=1,direction="both"):
            """
            @param direction: can be both, from, or to
            @type direction: string"""
            list_of_related_id_dicts = self.database_instance.getListOfRelatedIdDictsByIdDict({db:id},relation_type,depth,direction)
            return [dict([[database.internal_to_miriam_names[key], val] if key in database.internal_to_miriam_names else [key, val] for key, val in iddict.items()]) for iddict in list_of_related_id_dicts]
        def getRelationIS(self,db='',id=''):
            '''
            
            @return: dictionary of db id pairs that have the is relationship with the inserted pair
            @rtype: {db:id}
            '''
            return self.database_instance.getAnnotationObjectByIdDict({db:id}).toMiriamCompliantDatabase2id()
        def getRelationCHILD(self,db='',id='',depth=1):
            '''
            return all childred of the inserted pair up to a specified recursion depth
            @param depth: recursion depth (we need some convention here like -1 is infinite or 0)
            @rtype: [{db:id}]
            '''
            return self.getRelation(db,id,"is_a",depth,direction="from")
        def getRelationPARENT(self,db='',id='',depth=1):
            '''
            see L{getRelationCHILD}
            '''
            return self.getRelation(db,id,"is_a",depth,direction="to")
        def getName(self,db='',id='', try_direct_name=False):
            '''
            @return: return a human readable string representation of the db id pair
            @rtype: str
            '''
            key = db+":"+id
            if key not in self.db_id_to_name:
                self.db_id_to_name[key] = self.database_instance.__getNameByIdDict__({db:id}, try_direct_name)
            self.name_to_dblinks[self.db_id_to_name[key]] = {db:id}
            return self.db_id_to_name[key]
        def getListOfAllRelationTypes(self):
            """Returns a list of strings denoting all the relations in the database."""
            return ["all"]+self.database_instance.getListOfAllRelationTypes()
        def getListOfLinkedDatabases(self):
            """Returns a list of strings denoting all the relations in the database."""
            instances = ["_".join(x.split("_")[2:]) for x in dir(self.database_instance) if x.startswith("external_item") and eval("len(self.database_instance."+x+".rows)>0")]
            instances = [database.internal_to_miriam_names[x] if x in database.internal_to_miriam_names else x for x in instances]
            return instances
        def __updateRelationFunctions__(self):
            """Dynamically creates new getRelationX functions for all relation_types in the database."""
            relation_types = self.database_instance.getListOfAllRelationTypes()
            for relation_type in relation_types:
                function_name = "getRelation"+relation_type.upper()
                if function_name not in dir(self):
                    def wrap(self,db,id,depth=1,direction="both",relation_type=relation_type):
                        return eval("self.getRelation(db,id,relation_type,depth,direction)")
                    wrap.__name__ = function_name
                    setattr(self.__class__,wrap.__name__,wrap)
    
class DatabaseError(Exception):

    def __init__(self,message):
        self.errorMessage = message

    def __repr__(self):
        return self.errorMessage

    def __str__(self):
        return self.errorMessage


