#!/usr/bin/env python

# import external modules
#import libsbml
import urllib.request, urllib.error, urllib.parse, sys, os, os.path
from . import database, buildDB

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
            #TODO add at least one backup resource for each of the files!
            self.downloadable_databases = {"Taxonomy":{"download_sources":["Small Taxonomy Nodes","Small Taxonomy Names"], "required_ram":18, "download_size":1},
                               "ChEBI":{"download_sources":["ChEBI"], "required_ram":270, "download_size":13},
                               "GeneOntology":{"download_sources":["GeneOntology"], "required_ram":210, "download_size":16},
                               "KEGG Compound":{"download_sources":["KEGG Compound"], "required_ram":200, "download_size":30},
                               "KEGG Drug":{"download_sources":["KEGG Drug"], "required_ram":90, "download_size":18},
                               "KEGG Reaction":{"download_sources":["KEGG Reaction"], "required_ram":70, "download_size":4},
                               "KEGG Enzyme":{"download_sources":["KEGG Enzyme"], "required_ram":130, "download_size":26},
                               "KEGG Genome":{"download_sources":["KEGG Genome"], "required_ram":12, "download_size":1},
                               "Reactome":{"download_sources":["Reactome"], "required_ram":33, "download_size":2},
                               "Systems Biology Ontology":{"download_sources":["Systems Biology Ontology"], "required_ram":8, "download_size":1},
                               "Entrez Gene":{"download_sources":["Entrez Gene Refseq","Entrez Gene Go","Entrez Gene Uniprot","Entrez Gene Info"], "required_ram":10000, "download_size":260},
                               "Uniprot":{"download_sources":["Uniprot Human"], "required_ram":0, "download_size":50},
                               "SGD":{"download_sources":["SGD"], "required_ram":0, "download_size":3},
                               "Interpro":{"download_sources":["Interpro"], "required_ram":0, "download_size":75}
                               }
            self.first_build_databases = ["libsbml_resources", "listofresources", "element_score_matrix",
                              "MIRIAM_resources", "SBO_resources", "qualifiers"]
            self.default_databases = ["Taxonomy", "ChEBI", "GeneOntology", "KEGG Compound", "KEGG Reaction", "Systems Biology Ontology", "Interpro"]
            self.order_of_download_sources = ["libsbml_resources", "listofresources", "element_score_matrix",
                              "MIRIAM_resources", "SBO_resources", "qualifiers",
                              "Small Taxonomy Nodes", "Small Taxonomy Names",
                              "ChEBI", "GeneOntology", "KEGG Compound",
                              "KEGG Enzyme", "KEGG Drug", "KEGG Reaction", "KEGG Genome",
                              "Reactome", "Systems Biology Ontology",
                              "Entrez Gene Refseq", "Entrez Gene Go", "Entrez Gene Uniprot", "Entrez Gene Info",
                              "Uniprot Human", "SGD", "Interpro"
                              ]
            self.download_sources ={
            'Uniprot Human':{
                'web':[{'url':'ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_human.dat.gz'}],
                'integration_class':'UniprotRead',
                'belongs_to_database':'Uniprot'
                },
            'Entrez Gene Refseq':{
                'web':[{'url':'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz'}],
                'integration_class':'EntrezGene2RefSeqRead',
                'belongs_to_database':'Entrez Gene'
                },
            'Entrez Gene Go':{
                'web':[{'url':'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'}],
                'integration_class':'EntrezGene2GoRead',
                'belongs_to_database':'Entrez Gene'
                },
            'Entrez Gene Uniprot':{
                'web':[{'url':'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz'}],
                'integration_class':'EntrezGeneRefseqUniprotRead',
                'belongs_to_database':'Entrez Gene'
                },
            'Entrez Gene Info':{
                'web':[{'url':'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'}],
                'integration_class':'EntrezGeneInfoRead',
                'belongs_to_database':'Entrez Gene'
                },
            'Systems Biology Ontology':{
                'web':[{'url':'http://www.ebi.ac.uk/sbo/exports/Main/SBO_OBO.obo'}],
                'integration_class':'SboRead',
                'belongs_to_database':'Systems Biology Ontology'
                },
            'Reactome':{
                'web':[{'url':'http://www.reactome.org:8000/cgi-bin/remoteattsearch?DB=gk_current&SUBMIT=Search&QUERY_CLASS=PhysicalEntity&ATTRIBUTE%230=stableIdentifier&OPERATOR%230=IS+NOT+NULL&OUTPUTINSTRUCTIONS=stableIdentifier.identifier,_displayName,referenceEntity.identifier,referenceEntity.crossReference.identifier'}],
                'integration_class':'ReactomeRead',
                'belongs_to_database':'Reactome'
                },
            'Small Taxonomy Nodes':{
                'web':[{'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/ncbi_taxonomy_small_nodes.dmp'}],
                'integration_class':'NcbiTaxonomyNodesRead',
                'belongs_to_database':'Taxonomy'
                },
            'Small Taxonomy Names':{
                'web':[{'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/ncbi_taxonomy_small_names.dmp'}],
                'integration_class':'NcbiTaxonomyNamesRead',
                'belongs_to_database':'Taxonomy'
                },
            'ChEBI':{
                'web':[{'url':'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'}],
                'integration_class':'ChebiRead',
                'belongs_to_database':'ChEBI'
                },
            'GeneOntology':{
                'web':[{'url':'http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo'}],
                'integration_class':'GeneOntologyRead',
                'belongs_to_database':'Gene Ontology'
                },
            'KEGG Compound':{
                'web':[{'url':'ftp://ftp.genome.jp/pub/kegg/ligand/compound/compound'}],
                'integration_class':'KeggCompoundRead',
                'belongs_to_database':'KEGG Compound'
                },
            'KEGG Drug':{
                'web':[{'url':'ftp://ftp.genome.jp/pub/kegg/medicus/drug/drug'}],
                'integration_class':'KeggDrugRead',
                'belongs_to_database':'KEGG Drug'
                },
            'KEGG Reaction':{
                'web':[{'url':'ftp://ftp.genome.jp/pub/kegg/ligand/reaction/reaction'}],
                'integration_class':'KeggReactionRead',
                'belongs_to_database':'KEGG Reaction'
                },
            'KEGG Enzyme':{
                'web':[{'url':'ftp://ftp.genome.jp/pub/kegg/ligand/enzyme/enzyme'}],
                'integration_class':'KeggEnzymeRead',
                'belongs_to_database':'KEGG Enzyme'
                },
            'KEGG Genome':{
                'web':[{'url':'ftp://ftp.genome.jp/pub/kegg/genes/genome'}],
                'integration_class':'KeggGenomeRead',
                'belongs_to_database':'KEGG Genome'
                },
            'SGD':{
                'web':[{'url':'http://downloads.yeastgenome.org/chromosomal_feature/dbxref.tab'}],
                'integration_class':'SGDCrossRefRead',
                'belongs_to_database':'SGD'
                },
            'Interpro':{
                'web':[{'url':'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz'}],
                'integration_class':'InterproRead',
                'belongs_to_database':'Interpro'
                },
            'libsbml_resources':{
                'location':self._cfg.getpath('resources','libsbml_resources'),
                'web':[{'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/libsbmlresources.xml'}]
                },
            'listofresources':{
                'location':self._cfg.getpath('resources','listofresources'),
                'web':[{'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/listofresources.xml'}]
                },
            'element_score_matrix':{
                'location':self._cfg.getpath('resources','element_score_matrix'),
                'web':[{'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/elementscore.txt'}]
                },
            'qualifiers':{
                'location':self._cfg.getpath('resources','qualifiers'),
                'web':[{'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/qualifiers.csv'}]
                },
            'MIRIAM_resources':{
                'location':self._cfg.getpath('resources','MIRIAM_resources'),
                'web':[
                    {'url':'http://www.ebi.ac.uk/miriam/main/XMLExport', 'post_data':'getNameParam=MIRIAM_resources&submit=Export!'},
                    {'url':'http://semanticsbml.svn.sourceforge.net/viewvc/semanticsbml/trunk/resources/MIRIAM_resources.xml'}]
                },
            'SBO_resources':{
                'location':self._cfg.getpath('resources','SBO_resources'),
                'web':[{'url':'http://www.ebi.ac.uk/sbo/exports/Main/SBO_OBO.obo'}]                             }
            }

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
                    raise DatabaseError('Annotation database not found in %s' % self.base_dir)
            self.__updateRelationFunctions__()
            #---------------------------------------------------
            # caches
            self.db_id_to_name = {}
            self.name_to_dblinks = {}

        def forceReload(self):
            self.database_instance = self.database_instance.load()
            self.__updateRelationFunctions__()

        def __resourceNameToDownloadableFiles__(self, name):
            return_names = []
            for resource_name in self.download_sources:
                if 'belongs_to_database' not in self.download_sources[resource_name]:
                    continue
                if self.download_sources[resource_name]['belongs_to_database'] == name:
                    return_names.append( resource_name )
            return return_names

        def buildDB(self,resource_name=None,force_download=False):
            '''
            @param resource_name: specify a certain resource that should be build
            @type resource_name: str
            '''


            # make sure that the folder for the database has been created
            if not os.path.isdir(self.base_dir):
                for dir_name in [os.sep.join(self.base_dir.split(os.sep)[:-2]),
                         self.base_dir]:
                    if not os.path.isdir(dir_name):
                        os.mkdir(dir_name)
            # determine if we need to write the database to the hard_disk
            database_has_changed = False
            #---------------------------------------------------
            resource_order = self.order_of_download_sources
            resources = self.download_sources
            #---------------------------------------------------
            #if not self.yes_or_no_fkt("The internal database of semanticSBML has to be (re)built. This only works if you are connected to the internet.\n\nAre you connected to the Internet?"): raise DatabaseError('No Internet Connection')
            if resource_name:
                if type(resource_name) == type([]):
                    self.status_fkt("Building Database Resources: %s \n"%", ".join(resource_name))
                    resources = {}
                    for rn in resource_name:
                        try:
                            resources[rn] = self.download_sources[rn]
                        except KeyError:
                            # rn not in download_sources
                            if rn in [self.download_sources[x]['belongs_to_database'] for x in self.download_sources if 'belongs_to_database' in self.download_sources[x]]:
                                # translate into filenames
                                for file_name in self.__resourceNameToDownloadableFiles__(rn):
                                    resources[file_name] = self.download_sources[file_name]
                            else:
                                raise DatabaseError('Database Resource %s Unknown!' % rn)
                else:
                    self.status_fkt("Building Single Database Resource %s \n"%resource_name)
                    resources = {}
                    try:
                        resources = {resource_name:self.download_sources[resource_name]}
                    except KeyError:
                        # resource_name not in download_sources
                        if resource_name in [self.download_sources[x]['belongs_to_database'] for x in self.download_sources if 'belongs_to_database' in self.download_sources[x]]:
                            # translate into filenames
                            for file_name in self.__resourceNameToDownloadableFiles__(resource_name):
                                resources[file_name] = self.download_sources[file_name]
                        else:
                            raise DatabaseError('Database Resource %s Unknown!' % resource_name)
            else:
                # minimal selection of resources
                resources = [resource_name for resource_name in resources if resource_name in self.first_build_databases]
                if self.yes_or_no_fkt("No resources selected for building process\nDo you want to use the default selection?\nOtherwise you have to select resources manually?"):
                    for resource_name in self.order_of_download_sources:
                        if 'belongs_to_database' not in self.download_sources[resource_name]:
                            # default downloadable file, does not need to be appended in this step
                            continue
                        bigger_database_name = self.download_sources[resource_name]['belongs_to_database']
                        if bigger_database_name in self.default_databases:
                            resources.append(resource_name)
                else:
                    # add more resources to the selection
                    for resource_name in self.download_sources:
                        if resource_name in resources:
                            continue
                        if self.yes_or_no_fkt("Do you want to build resource %s"%resource_name):
                            resources.append(resource_name)
            #---------------------------------------------------
            # remove those resources already present on the disk
            present_databases = self.getListOfImportedDatabases()
            if not force_download:
                # removing resources from the big database
                resources = [resource_name for resource_name in resources
                         if not ( 'belongs_to_database' in self.download_sources[resource_name] and
                              self.download_sources[resource_name]["belongs_to_database"] in present_databases)
                         ]
                # removing things with a location already present on the disk
                resources = [resource_name for resource_name in resources
                         if not ( 'location' in self.download_sources[resource_name] and
                              os.path.exists(self.download_sources[resource_name]["location"]))
                         ]
            resources = dict([[name,self.download_sources[name]] for name in resources])
            
            download_files = [self.download_sources[resource_name]["web"][-1]["url"] for resource_name in resources]
            # skip the function if there are no files to download
            if len(download_files) == 0:
                resources = {}
            if resources != {}:
                if not self.yes_or_no_fkt("""
CREATING LOCAL DATABASE

The tools of the semanticSBML suite need a local database of names
with links to external databases to operate. In order to build
this local database you need to be authorized to access some
external databases. At the moment, KEGG and Gene Ontology are
mandatory, so please check their license agreements. 

It is possible that this tool may not be able to download the
necessary files automatically, e.g. if you need to use a proxy server
to connect to the internet. In this case please download the files
yourself and start this tool again from the directory where they are
stored. You would need to retrieve the following files:

  %s

Can I continue?

The database files will be created in "%s".

This process will take some time, so please be patient.""" % ('\n  '.join(download_files), self.base_dir)):
                    return False

            
            # --------------------------------------------------
            errorstr=''
            for name in resource_order:
                if name not in resources:
                    continue
                resource = resources[name]
                i=0
                while i<len(resource['web']):
                    post_data=None
                    if 'post_data' in resource['web'][i]: post_data = resource['web'][i]['post_data']
                    try:
                        internet_contents = self._download(resource['web'][i]['url'],post_data)
                    except urllib.error.URLError as e:
                        err="Error Getting %s Because:\n\t%s\nIf Available: Trying Next Resource...\n"%(resource['web'][i]['url'],e)
                        errorstr+=err
                        self.status_fkt(err)
                        i+=1
                        continue
                    except ContentTooShortError as e:
                        err="Error Getting %s Because:\n\t%s\nIf Available: Trying Next Resource...\n"%(resource['web'][i]['url'],e)
                        errorstr+=err
                        self.status_fkt(err)
                        i+=1
                        continue
                    except IOError as e:
                        err="Error Getting %s Because:\n\t%s\nIf Available: Trying Next Resource...\n"%(resource['web'][i]['url'],e)
                        errorstr+=err
                        self.status_fkt(err)
                        i+=1
                        continue
                    else:
                        if internet_contents:
                            breaking = False
                            if "location" in resource:
                                open(resource['location'],'w').write(internet_contents)          #'SBO_OBO.obo','w')
                                breaking = True
                            if "integration_class" in resource:
                                # process the contents
                                self.status_fkt("Processing data for database "+name+".",False,True)
                                IntegrationClass = buildDB.__dict__[resource['integration_class']]
                                IntegrationClass(self.database_instance,internet_contents)
                                database_has_changed = True
                                self.database_instance.addToListOfImportedDatabases(name)
                                breaking = True
                            if breaking:
                                break
                        else:
                            err='Error: No content found'
                            errorstr+=err
                            self.status_fkt(err)
                            i+=1
                            continue
                else:
                    self.status_fkt( "Could Not Get %s\n%s\n"%(name,errorstr), False, True)
            # if database_has_changed: # we want the database to be written anyways, maybe some of the structure of the database has been updated
            # write the new database contents to disk
            self.database_instance.dump()
            self.__updateRelationFunctions__()
            self.status_fkt("Finished Building The Database\n",True)
            self._cfg.write()
            return True

        def getListOfImportedDatabases(self):
            """Returns a list of imported databases. The names are identical with the names in the download sources."""
            return self.database_instance.getListOfImportedDatabases()

        def getListOfDownloadableDatabases(self):
            return list(self.downloadable_databases.keys())
        def checkConnection(self):
            '''
            TODO this did not work for me!
            Function checkConnection checks whether the computer is connected to the internet
            '''
            if os.system("ping www.google.com -c 2"): return True
            else: return False
        def _download(self,url,post_data=None):
            '''
            Function download is downloading and returning the database flatfiles.
            @param url     : url of a database flatfile
            @type url      : string            
            @param post_data : data for POST requests e.g. {'submit':'Download'}
            @type post_data : dict
            @return: downloaded page
            @rtype: str
            '''
            #-----------------------------------------------------
            #create a realistic user agent
            headers = {
            'User-Agent' : 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)',
            'Accept' :
            'text/xml,application/xml,application/xhtml+xml,text/html;q=0.9,text/plain;q=0.8,image/png,*/*;q=0.5',
            'Accept-Language' : 'fr-fr,en-us;q=0.7,en;q=0.3',
            'Accept-Charset' : 'ISO-8859-1,utf-8;q=0.7,*;q=0.7'
            }
            #-----------------------------------------------------
            #check if we got an internet connetction
            #self.status_fkt('Checking connection ...')
            #if not self.checkConnection(): raise DatabaseError('No Internet Connection')
            #-----------------------------------------------------
            self.status_fkt("Loading %s\n"%url)
            #-----------------------------------------------------
            reqObj = urllib.request.Request(url, post_data, headers)
            fp = urllib.request.urlopen(reqObj)
            headers = fp.info()
            ##    This function returns a file-like object with two additional methods:
            ##    * geturl() -- return the URL of the resource retrieved
            ##    * info() -- return the meta-information of the page, as a dictionary-like object
            ##Raises URLError on errors.
            ##Note that None may be returned if no handler handles the request (though the default installed global OpenerDirector uses UnknownHandler to ensure this never happens). 
            #-----------------------------------------------------
            #read & write fileObj to filename
            filename = 'outfile'
            outcontents=''
            result = filename, headers
            bs = 1024*8
            size = -1
            read = 0
            blocknum = 0
            #-----------------------------------------------------
            size_known = "content-length" in headers
            if size_known: 
                size_known = int(headers["Content-Length"])
                size = int(headers["Content-Length"])
            self.download_status_fkt(size_known,0,True,url=url)
            #-----------------------------------------------------
            while 1:
                block = fp.read(bs)
                if block == "":
                    break
                read += len(block)
                outcontents+=block
                blocknum += 1
                #-------------------------------------------
                self.download_status_fkt(size_known,(blocknum*bs),False,False)
            #-----------------------------------------------------
            self.download_status_fkt(size_known,size_known,False,True) 
            fp.close()
            del fp
            #-----------------------------------------------------
            # raise exception if actual size does not match content-length header
            if size >= 0 and read < size:
                raise ContentTooShortError("retrieval incomplete: got only %i out of %i bytes" % (read, size), result)

            #-----------------------------------------------------
            # taking care of gzipped content
            if url.endswith(".gz"):
                import io, gzip
                stream = io.StringIO(outcontents)
                zfile = gzip.GzipFile(fileobj=stream)
                outcontents = zfile.read()
            #-----------------------------------------------------

            return outcontents
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

class ContentTooShortError(DatabaseError):
    def __init__(self,message):
        self.errorMessage = 'Bla: %s'% message


if __name__ == "__main__":
    db = Meta_DB()
    db.buildDB()
    #print db.getRelationFromSource("Taxonomy","1",source="taxonomy")
