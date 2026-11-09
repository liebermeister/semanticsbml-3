#!/usr/bin/env python
'''
B{MIRIAM annotation manipulation} this module provides classes for the manipulation of MIRIAM annotations of an SBML model or single SBML elements
'''
import os
import sys
import xml
import libSBAnnotation.annotate,libSBAnnotation.database_new

#external libaries
import libsbml, re, urllib.request, urllib.parse, urllib.error
from . import xmlobject
#semanticSBML libraries
import libSBAnnotation.config


#////////////////////////////////////////////////////
ldb = None
nldb = None
#////////////////////////////////////////////////////
#for fkt in 'uri2name  name2uri database2elements elements2databases name2altURIs database2regex name2www get_db_from_uri'.split():
#get functions from moved annotate module
Annotation     = libSBAnnotation.annotate.Annotation

redownload         = libSBAnnotation.annotate.redownload
checkAndGetXml     = libSBAnnotation.annotate.checkAndGetXml

uri2name       = libSBAnnotation.annotate.uri2name#maps ever known (also deprecated) to the db name
name2uri       = libSBAnnotation.annotate.name2uri#this will map the db name to the correct uri
database2elements  = libSBAnnotation.annotate.database2elements
elements2databases = libSBAnnotation.annotate.elements2databases#maps libsbml element types to db
name2altURIs       = libSBAnnotation.annotate.name2altURIs#maps db name to deprecated uris
database2regex     = libSBAnnotation.annotate.database2regex
name2www       = libSBAnnotation.annotate.name2www#maps db name to weblink in fto
get_db_from_uri    = libSBAnnotation.annotate.get_db_from_uri
#////////////////////////////////////////////////////
AnnotationError     = libSBAnnotation.annotate.AnnotationError
AnnotationDatabaseError = libSBAnnotation.annotate.AnnotationDatabaseError
AnnotationCreationError = libSBAnnotation.annotate.AnnotationCreationError
#////////////////////////////////////////////////////
#build lookup dictionaries
#------------------------------------------------------------------------------------
bio_qualifier_type2name   = {}
model_qualifier_type2name = {}
bio_qualifier_name2type   = {}
modelqual_name2type       = {}
modelqual_name2desc       = {}
bioqual_name2desc         = {}
bio_qualifiers        = []
model_qualifiers      = []

#////////////////////////////////////////////////////
x = checkAndGetXml('resources','libsbml_resources')
#////////////////////////////////////////////////////
try:
    for bqr in x.root.BioQualifier[0].qualifier:
        try:
            bioqual_name2desc[bqr.name]=bqr.description
            bio_qualifiers.append(bqr.name)
            bio_qualifier_name2type[bqr.name] = getattr(libsbml,bqr.var)
            bio_qualifier_type2name[getattr(libsbml,bqr.var)] = bqr.name
        except AttributeError:
            pass #new qualifier but old libsbml
    for mqr in x.root.ModelQualifier[0].qualifier:
        try:
            modelqual_name2desc[mqr.name] = mqr.description
            model_qualifiers.append(mqr.name)
            model_qualifier_type2name[getattr(libsbml,mqr.var)] = mqr.name
            modelqual_name2type[mqr.name] = getattr(libsbml,mqr.var)
        except AttributeError:
            pass #new qualifier but old libsbml
except KeyError as e:
    msg = '%e<br>Old %s needs update!<br> Trying to get it ...'%(e,libSBAnnotation.config.Config().getpath('resources','libsbml_resources'))
    sys.stderr.write(msg)
    redownload(msg,'libsbml_resources')
    raise AnnotationDatabaseError("Ignore the Line Above<br>Everything Went Fine ... <br>Please Restart the Process<br>Yes this look a bit messy, we are working on it!")
#------------------------------------------------------------------------------------

class Annotator:
    '''convenience class for L{ModelElementsAnnotations}, L{ElementAnnotations} and L{Annotation}''' 

    def __init__(self,libsbml_model=None):
        #deprecated_module_remove_code_this_is_supposed_to_crash
        self.libsbml_model=libsbml_model

    def getElementAnnotations(self,libsbml_model=None):
        return ModelElementsAnnotations(self.libsbml_model).getElementAnnotations()

    def getNumNotAnnotatedElements(self,libsbml_model):
        return ModelElementsAnnotations(libsbml_model).getNumNotAnnotatedElements()

    def addAnnotationLink(self,element,db,id,qualifier=None):
        ElementAnnotations(element).addAnnotation(Annotation(db,id,qualifier))

    def addAnnotationLinkAutomatic(self,element):
        ElementAnnotations(element).addAnnotationAutomatic()

    def remAnnotationLink(self,element,db,id):
        ElementAnnotations(element).remAnnotation(db,id)

    def remAnnotation(self,element):
        pass
    def getAnnotationSuggestions(self,element,query=''):
        if query:
            return ElementAnnotations(element).getSuggestions()
        else:
            ae=ElementAnnotations(element)
            suggestions=[]
            for q in ae.getQuerys():
                suggestions+=ae.getSuggestions(q)
            return suggestions
    def getAnnotationQuerys(self,element):
        return ElementAnnotations(element).getQuerys() 

class ModelElementsAnnotations:
    '''
    Represents all L{ElementAnnotations}s (each contains a list of L{Annotation}s)
    functions that operate on multiple L{ElementAnnotations}s
    '''
    def __init__(self,libsbml_model,suppress_errors=None, libsbml_document = None):
        '''
        Read the SBML file elements that could be annoated.
        Elements are stored as L{ElementAnnotations}s
        @param suppress_errors: 
        @type suppress_errors: bool 
        '''
        self.libsbml_document = libsbml_document
        self.libsbml_model = libsbml_model if not libsbml_document else libsbml_document.getModel()
        self.suppress_errors = suppress_errors
        self._annotationElements =  []
        self._readElementAnnotations(self.libsbml_model)

    def getNumNotAnnotatedElements(self):
        '''
        return the number of elements that have no MIRIAM annotations
        @return: number of not MIRIAM annotated elemtns
        @rtype: int
        '''
        return len([x for x in self._annotationElements if not x.isAnnotated()])    

    def getElementAnnotations(self):
        '''
        return a list of L{ElementAnnotations} that can contain MIRIAM annotations
        @return: return all annotatable elements of the model
        @rtype: [L{ElementAnnotations}]
        '''
        return self._annotationElements

    def remElementAnnotations(self,element):#TODO not finished yet must delte from sbml too
        '''
        remove an element
        '''
        self._annotationElements.pop(element)
    
    def _readElementAnnotations(self,libsbml_model):
        '''
        go through all MIRIAM annotatabel elements in a libsbml model
        @type libsbml_model: L{libsbml.Model}
        @param libsbml_model: libsbml Model instance
        '''
        self._annotationElements.append(ElementAnnotations(libsbml_model,self.suppress_errors))
        error_str = ''
        for type in 'CompartmentTypes Compartments SpeciesTypes Species Parameters Reactions Events Rules'.split():
            for e in getattr(libsbml_model,'getListOf'+type)():
                try:
                    self._annotationElements.append(ElementAnnotations(e,self.suppress_errors))
                except AnnotationError as e:
                    error_str += str(e)+"\n"
        if error_str and not self.suppress_errors:
            raise AnnotationError(error_str)
        self._annotationElements.reverse()

class ElementAnnotations:
    '''
    this class represents the MIRIAM annotations of a libsbml element
    '''

    def __init__(self,libsbml_element=None,suppress_errors=None):
        '''
        read MIRIAM annotations, type, id, metaid and name from libsbml element
        not all libsbml elements are allowed as input

        @param libsbml_element: libsbml element from eg listOfSpecies()
        @type libsbml_element: libsbml.species, libsmbl.compartment OR ...
        @raise AnnotationError: error if wrong type is inserted
        @param suppress_errors: if is set all read errors will be suppressed 
        @type suppress_errors: bool
        '''
        #-----------------------------------------------------------------
        global nldb
        if not nldb: nldb = libSBAnnotation.database_new.Meta_DB()
        #-----------------------------------------------------------------
        for t in [libsbml.SpeciesType,
                  libsbml.Species,
                  libsbml.CompartmentType,
                  libsbml.Compartment,
                  libsbml.Reaction,
                  libsbml.Compartment,
                  libsbml.AssignmentRule,libsbml.RateRule,libsbml.AlgebraicRule,
                  libsbml.Event,
                  libsbml.Model,
                  libsbml.Parameter,
                  libsbml.UnitDefinition]:
            if isinstance(libsbml_element, t):
                break
        else:
            raise AnnotationError('ElementAnnotations initialized with wrong type, you inserted %s' % libsbml_element)
        if not libsbml_element.isSetMetaId():#FIXME may not be uique
            if type(libsbml_element) in [libsbml.AssignmentRule,libsbml.RateRule,libsbml.AlgebraicRule,libsbml.Event]:
                suffix = '_noid_meta'
            else:
                suffix = '_meta'
            libsbml_element.setMetaId(libsbml_element.getId()+suffix)
        
        self.libsbml_element = libsbml_element
        #-----------------------------------------------------------------

        #sbml things we need to know (since rules get identical ids as their species, this leads to functional
        #problems. we therefore extend their id to "ID+ [~ Rule]")
        #-----------------------------------------------------------------
        self.type   = libsbml_element.getElementName()
        if self.type == 'assignmentRule':
            self.id = libsbml_element.getId() +' [Assignment Rule]'
        elif self.type == 'rateRule':
            self.id = libsbml_element.getId() + ' [Rate Rule]'  
        elif self.type == 'algebraicRule':
            self.id = libsbml_element.getId() + ' [Algebraic Rule]'
        else:
            self.id = libsbml_element.getId()
        self.metaid     = libsbml_element.getMetaId()
        self.name   = libsbml_element.getName()
        self.repr   = self.name or self.id
        #-----------------------------------------------------------------
        self.suppress_errors=suppress_errors
        self.error_log  =''
        #-----------------------------------------------------------------
        self._readAnnotations()


    def checkMeOut(self):
        '''TODO
        this function is not the way to go
        what we need to do is to create a baseclass 
        then inherit from this so we can do a isinstance, but this is also not a good solution, using isinstance i usually discouraged
        '''
        return 'MIRIAM'


    def search_id_by_name_exact(self,name):
        return self.getSuggestions(name,exact_search=True)


    def unsetAnnotations(self):
        '''
        delete all MIRIAM annotations of this element
        '''
        self._removeResources(lambda resource: True)
        self._readAnnotations()

    def getAnnotations(self,filter=None):
        '''
        @return: List of L{Annotation} class instances
        @rtype:[Annotation]
        '''
        if filter:
            ret=[]
            for db in filter:
                db=get_db_from_uri(db)
                for a in self.annotations:
                    if a.db==db:
                        ret.append(a)
            return ret
        #self._readAnnotations()#!!! this was very bad, never do this here since it seg faults libsbml!!!
        return self.annotations 

    def addAnnotationAutomatic(self):
        '''check if the elements is already annotated with L{isAnnotated} if not check if a annotation link can be found using L{getQuerys} and L{getSuggestions} and if only one link is found add the annotation with L{addAnnotation}'''

        if self.isAnnotated():
            raise AnnotationError('Element %s already annotated' % self.name)

        else:
            suggestions=[]
            anFound=False#this is now a flag it was supposed to filter too many result by name differences but this does not make sense any more
            query = self.name or self.id
            #print 'query: ',query#DEBUG

            if not query:#lame this will never happen, id is required :(
                raise AnnotationError('no search query found')

            for an in self.getSuggestions(query,True):
                anFound=True
                #elif tmp_name!=name:
                #   raise AnnotationError('Too many suggestions found for %s' % self.name)
                suggestions.append(an)

            if not anFound:
                raise AnnotationError('No suggestions found for %s ' % query)
            else:
                for an in suggestions:
                    self.addAnnotation(an)

    def addAnnotation(self,annotation, suppress_errors=False):
        '''
        add a MIRIAM annotation with qualifier to the sbml element
        check if a CVTerm with the right qualifier already exists and add annotation there or 
        create a new CVTerm
        @param annotation: instance of an L{Annotation} class
        @type annotation: L{Annotation}

        '''
        new_cvt = libsbml.CVTerm()
        qtypes  = {'model':0,'bio':1,'unknown':2}                             #TODO move to xml
        new_cvt.setQualifierType(qtypes[annotation.qualifier_type])

        if annotation.qualifier_type=='bio':
            new_cvt.setBiologicalQualifierType(bio_qualifier_name2type[annotation.qualifier])
        elif annotation.qualifier_type=='model': 
            new_cvt.setModelQualifierType(modelqual_name2type[annotation.qualifier])
        elif annotation.qualifier_type=='unknown':
            pass
        else: 
            raise AnnotationCreationError('Unknown qualifier type in annotation %s'%annotation)

        new_cvt.addResource(str(annotation))
        self.libsbml_element.addCVTerm(new_cvt)
        self._readAnnotations()
        #self.annotations.append(annotation)



    def modAnnotation(self,annotation):
        '''
        modify the qualifier of an annotation
        @param annotation: L{Annotation} class instance with new qualifier but old link
        @type annotation: L{Annotation}
        '''
        self.remAnnotation(annotation)
        #print 'modifying2 to ',annotation.db, annotation.id, annotation.qualifier_type,annotation.qualifier#DEBUG
        self.addAnnotation(annotation)

    def remAnnotation(self,annotation):
        '''
        delete a MIRIAM annotation from libsbml element and resynchronize the annotations with libsbml element
        @param annotation: L{Annotation} class instance 
        @type annotation: L{Annotation}
        '''
        spellings = set([annotation.safe_str()] + annotation.getResourceStrings())
        try:
            for alturi in name2altURIs[annotation.db]:
                spellings.add(str('%s#%s'%(alturi,annotation.id)))
        except KeyError:
            pass
        self._removeResources(lambda resource: resource in spellings)
        self._readAnnotations()

    def _removeResources(self,remove_this):
        '''
        remove all resources for which the function remove_this(resource string) is true.
        A CVTerm without resources must not stay in the element (current libSBML
        crashes when writing it), therefore the list of CVTerms is rebuilt.
        '''
        kept = []
        changed = False
        for n in range(self.libsbml_element.getNumCVTerms()):
            cvt = self.libsbml_element.getCVTerm(n).clone()
            res = cvt.getResources()
            for resource in [res.getValue(m) for m in range(res.getLength())]:
                if remove_this(resource):
                    cvt.removeResource(resource)
                    changed = True
            if cvt.getResources().getLength():
                kept.append(cvt)
        if changed:
            self.libsbml_element.unsetCVTerms()
            for cvt in kept:
                self.libsbml_element.addCVTerm(cvt)

    def _readAnnotations(self):
        '''
        read all MIRIAM annotations (from CVTerms) in the libsbml element
        '''
        self.annotations = []
        error_str        = ''
        for n in range(self.libsbml_element.getNumCVTerms()):
            cvt = self.libsbml_element.getCVTerm(n)
            res = cvt.getResources()
            for m in range(res.getLength()):
                resource = res.getValue(m)
                db = ''
                id = ''
                try:
                    db,id = libSBAnnotation.annotate.split_resource(resource)
                except ValueError as e:
                    error_str+='Annotation resource \'%s\' does not follow a known scheme (identifiers.org URL, urn:miriam:DB:ID or URI#ID). ' % resource
                    continue
                try:
                    if self.libsbml_element.getCVTerm(n).getQualifierType()==1:
                        self.annotations.append(Annotation(db,id,self.libsbml_element.getCVTerm(n).getQualifierType(),bio_qualifier_type2name[self.libsbml_element.getCVTerm(n).getBiologicalQualifierType()],True))
                    elif self.libsbml_element.getCVTerm(n).getQualifierType()==0:
                        self.annotations.append(Annotation(db,id,self.libsbml_element.getCVTerm(n).getQualifierType(),model_qualifier_type2name[self.libsbml_element.getCVTerm(n).getModelQualifierType()],True))
                    elif self.libsbml_element.getCVTerm(n).getQualifierType()==2:
                        self.annotations.append(Annotation(db,id,self.libsbml_element.getCVTerm(n).getQualifierType(),None,True))
                except AnnotationCreationError as e:
                    error_str+=str(e)+' '
                except KeyError as e:
                    msg = 'Unknown Qualifier Number: %s,<br> Trying to get new %s ...'%(e,libSBAnnotation.config.Config().getpath('resources','libsbml_resources'))
                    redownload(msg,'libsbml_resources')
                    raise AnnotationDatabaseError("%s<br><b>Please Restart semanticSBML</b><br>"%msg)
        self.annotations.reverse()
        if error_str:
            if self.suppress_errors!= None:
                self.error_log=error_str
            else:
                raise AnnotationCreationError(error_str[:-1])

    def getReadErrors(self):
        return self.error_log

    def getQuerys(self):
        '''
        concatenate the name, id and metaid values of the libsbml element and return them
        @return: list of query strings
        @rtype: [str]
        '''
        return list(set([self.id,self.metaid,self.name]))

    def isAnnotated(self):
        '''
        return Trueif the elements contains MIRIAM annotations
        @return: is MIRIAM annotated
        @rtype: C{bool}
        '''
        for n in range(self.libsbml_element.getNumCVTerms()):
            if self.libsbml_element.getCVTerm(n).getResources().getLength():
                return True
        else:
            return False

    def getSuggestions(self,query=None,exact_search=False):#TODO incomplete due to database
        '''
        get a list of suggested L{Annotation}s
        if no query is specified the function L{getQuerys} is used
        @param query: search words
        @type query: str
        @param exact_search: search only for exact matches
        @type exact_search: C{bool}
        '''
        querys = None
        if query: querys=[query]
        else: querys=self.getQuerys()
        for q in querys:
            db_links=None
            if not query or exact_search:#do an exact search
                db_links = nldb.searchEntries(q.lower())
                if len(db_links)>1:#only return one group hits
                    db_links=[]
            else:
                db_links = nldb.searchEntries(q.lower(),libSBAnnotation.config.Config().getfloat('annotate','search_precision'))
            for link_group in db_links:
                for db,id in link_group.items():
                    try:
                        type = self.type
                        if 'rule' in type.lower(): type = 'rule' 
                        elif 'species' in type.lower(): type = 'species' 
                        if db in libSBAnnotation.config.Config().get('annotate','annotate_%s_with'%type).split(','):
                            yield Annotation(db,id,qualifier="is") # unknown as a qualifier made the model crash during the toSBML method, setting to "is" instead
                    except AnnotationCreationError as e:#our database is messed up, filter it here whatever, show must go on
                        pass
                        #sys.stderr.write( 'xWarning: '+str(e))




if __name__ == "__main__":

    document= libsbml.readSBMLFromString(open('/home/tlubitz/Desktop/sboclass/semanticsbml/trunk/models/test01.xml','r').read())
    model   = document.getModel() 
    ae=None
    '''
    try:
        ae  = ModelElementsAnnotations(model)   
    except AnnotationError,e:
        print 'model intitalization failed: ',e
    else:
        for elem in ae.getElementAnnotations():
            print 'ELEMENT: ',elem.isAnnotated(),' ',elem.name

            for an in elem.getAnnotations():
                print '\t\tAnnotation: ',an.db,an.id,an.qualifier,an.qualifier_type,str(an)
            try:
                elem.addAnnotationAutomatic()
            except AnnotationError,e:
                print '\t\t\tAutomatic Annotation Failed: ',e
        print ae.getNumNotAnnotatedElements()
    anno    = Annotation('Gene Ontology','GO:1234567','is')
    '''
    prick = ModelElementsAnnotations(model)
    print(prick.getElementAnnotations())


