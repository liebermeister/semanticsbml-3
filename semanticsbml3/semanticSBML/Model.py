"""
Interface to all functions in the semanticSBML core for one document actions

Example:
    m = semanticSBML.Model("bla.xml")
    m.set_callback(output)
    m.annotate()
    m_merged = m.merge(other_model)
    m_merged.get_libsbml_model()
"""

import sys
import types
import libsbml
import os
from . import annotate
import libSBAnnotation.database_new

#!/usr/bin/env python

class LibSbmlError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class NeedFilenameError(Exception):
    pass


class Model(object):

    def __new__(cls, document, callback = None, *args, **kwargs):
        '''
        initialize the model from either
            - string of the filename
            - a other model object
            - a libsbml SBMLDocument
            - None, create a new model
        '''
        if type(document) == cls:
            model = document
            if callback:
                model.set_callback(callback)
            return model

        model = object.__new__(cls)
        model.set_callback(callback)
        model.filename = None

        #if None was passed create a new empty model
        if document == None:
            model.document = libsbml.SBMLDocument()
            model.document.createModel()
        #TODO if a filehandle??? not sure here
        elif hasattr(document,"read"):
            try:
                model.filename = document.filename
            except AttributeError:
                pass
        #TODO shouldn't this be in an else case? if the above failes this will be done anyway????
            model.document = libsbml.readSBMLFromString(document.read())
            model._check_document(document.name)
        #if a filename was passed 
        elif type(document) == type(str()):#is document a string
            model.filename = document
            model.document = libsbml.readSBML(document)
            model._check_document(document)
        #if all fails assume a libsbml.SBMLDocument was passed
        else:
            model.document = document


        #######print model.document.toSBML()

        #check if the above really returned a SBMLDocument
        if not isinstance(model.document, libsbml.SBMLDocument):
            raise TypeError

        try:
            # libsbml3.0 compliance
            model.document.setLevelAndVersion(2, 3)
        except AttributeError:
            model.document.setLevel(2)
            model.document.setVersion(3)

        return model

    def __init__(self,*args):
        '''initialize the L{Annotator} and the L{Meta_DB} and the L{SboDatabase}'''
        self._annotator     = annotate.Annotator()
        self._ldb       = libSBAnnotation.database_new.Meta_DB()
        self.miriams = None
        self.sbos    = None

    def save_as(self, filename):
        '''
        save model as sbml file on the hd using the save function
        '''
        self.filename = filename
        (x,extension)=os.path.splitext(self.filename)
        #if extension.lower()!='.xml':
        #   self.filename += '.xml'
        self.save()

    def save(self):
        '''
        save the model on the hd
        '''
        if self.filename is None:
            raise NeedFilenameError
        import codecs
        codecs.open(self.filename,'w', encoding='utf-8').write(str('<?xml version="1.0" encoding="UTF-8"?>\n%s'%self.document.toSBML()))

    def _check_document(self, name):
        if self.document.getNumErrors() > 0:
            emessages=''
            for n in range(self.document.getNumErrors()):
                emessages+=self.document.getError(n).getMessage()
            raise LibSbmlError('Opening %s returned the following error: %s' % (name,emessages))


    def check_consistency(self):#TODO ever used !!!???
        if sys.platform == 'win32':
            return None
        self.document.checkConsistency()
        return self._get_errors()

    def convert(self):
        print(model.document.toSBML())

    def _get_errors(self):
        '''check for libsbml errors'''
        ret = dict()
        for t in ('Warning', 'Error', 'Fatal'):
            n = getattr(self.document, 'getNum' + t + 's')()
            if n:
                getter = getattr(self.document, 'get' + t)
                ret[t.lower() + 's'] = [getter(i).getMessage() for i in range(n)]
        return ret

    #
    # accessor methods
    # 
    def set_callback(self, callback):
        '''
        deprecated function to set a callback object
        '''
        self.callback = callback

    def get_callback(self):
        '''
        @return: callback object
        '''
        return self.callback

    def get_libsbml_document(self):
        '''
        @return: libsbml document instance
        '''
        return self.document
    def get_libsbml_model(self):
        '''
        @return: libsbml model instance
        '''
        return self.document.getModel()

    def get_model_as_svg(self):
        '''
        return a graphical representation of the model in the svg format
        this function needs dot to be installed (graphviz)
        @return: graph representation in svg
        @rtype: str
        '''
        import semanticSBML.sbml2dot
        #from sbml2dot import sbml2svg
        return semanticSBML.sbml2dot.sbml2svg(self.get_libsbml_document())


    def check(self):
        '''
        run a semantic check on the model
        @return: check results object
        @rtype: [L{CheckResult}]
        '''
        from .sbmlcheck import Checker
        checker = Checker(self)
        return checker.run()

    def addAnnotationLinkAutomatic(self,element,callback=None):
        '''
        add annotations automatically if not links are found throw an error
        '''
        annotate.ElementAnnotations(element).addAnnotationAutomatic()

    def addAnnotationLink(self,element,db,id,qtype=None,qual=None):
        '''
        add an annotation link (Database Identifier) to a specified element 
        if not found throw L{InvalidAnnotationError}
        '''
        an=None
        if qtype:   an = annotate.Annotation(db,id,qualifier_type=qtype,qualifier=qual) 
        else:       an = annotate.Annotation(db,id,'unknown')
        annotate.ElementAnnotations(element,True).addAnnotation( an )

    def addSboTerm(self,element,id):
        '''
        add SBO term to a specified element 
        if not found throw L{InvalidAnnotationError}
        '''
        #print 'in Model'
        from . import sbo
        sbo.SboElement(element).setSbo(id)
    def delAnnotationLink(self,element,db,id):
        '''
        delete an annotation link, see L{addAnnotationLink}
        '''
        annotate.ElementAnnotations(element).delAnnotation( annotate.Annotation(db,id) )

    def delSboTerm(self,element,id):
        '''
        delete SBO Term
        '''
        from . import sbo
        sbo.SboElement(element).unsetSboTerm()                

    def getNumNotAnnotatedElements(self):
        '''
        @return: number of elements that are not annotated
        '''
        return str(annotate.ModelElementsAnnotations(self.get_libsbml_model(),'').getNumNotAnnotatedElements())

    def getElementAnnotations(self,error_log=None):
        '''
        return the list of annotation elements of this model
        @return: list of sbml elemtns that can be annotated
        '''
        from . import sbo
        tuplelist = []
        #TODO this is not a good solution since we don't get indiviual errors every time, could also be a problem because it does not update properly
        if not self.miriams:
            self.miriams    = annotate.ModelElementsAnnotations(self.get_libsbml_model(),error_log).getElementAnnotations()
        if not self.sbos:
            self.sbos       = sbo.ModelSbo(self.get_libsbml_model()).getElementsSbo()
        tuplelist  = [self.miriams,self.sbos]
        return tuplelist


    def getAnnotationSuggestions(self,element,query=None):
        '''
        return a list of suggested annotations
        @return: suggestion for the annotation of sbml elements
        @rtype: [(str_name,(str_db,str_id)]
        '''
        return annotate.ElementAnnotations(element).getSuggestions(query)
    def getAnnotationQuerys(self,element):
        '''
        return a list of of query strings that can be used to search MIRIAM annotations in the internal database
        @return: list of search query suggestions
        @rtype: [str]
        '''
        #o=[]
        #for n in self._annotator.get_names_from_element(element):
        #   o+=n
        #return o
        return annotate.ElementAnnotations(element).getQuerys()
    def getAnnotation(self,element):
        '''
        return the MIRIAM annotation of the inserted libsbml element
        @return: list of databasenames and identifiers found in the element
        @rtype: [annotate.Annotation]
        '''
        returnlist = []
        miriams    = annotate.ElementAnnotations(element).getAnnotations()
        returnlist = miriams# + sbos
        return returnlist
    def issetAnnotation(self,element):
        '''
        this function should not be used but L{getAnnotationStatus}
        @return: does the element contain annoations
        @rtype: bool
        '''
        return annotate.ElementAnnotations(element).isAnnotated()
    def id2str(self,db,id):
        '''
        return a string representation of the database id couple
        @return: str representation of the identifier
        @rtype: str
        '''
        return self._ldb.getName(db,id)
    def id2str_compartment(self,id):
        '''
        see L{id2str}
        '''
        return self._ldb.getName("GeneOntology", id)
    def id2str_reaction(self,id):
        '''
        see L{id2str} but for reactions
        '''
        return self._ldb.getReactionDescription("kegg_reaction", id)
    def compartmentIds(self):
        '''
        @return: list of all compartment ids
        @rtype: [str]
        '''
        return [x[0] for x in self._ldb.database_instance.external_item_go.rows]
    def databases(self,typeCode=None):
        '''
        return the list of databases used in the internal database
        if no input return all databases 
        else return the databases that match the libsbml object type
        @return: list of database names
        @rtype: [str]
        '''
        return list(self._ldb.downloadable_databases.keys())

    #let's not include merge here as it needs at least 2 documents as input 
    #def merge(self, other_model):
    #    pass






    # Functions for SBO Terms, using semanticSBML/sbo.py

    def getElementSbos(self):
        '''
        return the list of sbo elements of this model
        @return: list of sbo elements
        '''
        #        return annotate.ModelElementsAnnotations(self.get_libsbml_model()).getElementAnnotations()
        from . import sbo
        return sbo.ModelSbo(self.get_libsbml_model()).getElementsSbo()

    def getSbo(self,element):
        '''
        return SBO Term of the inserted libsbml element
        @return: SBO Term
        @rtype : string
        '''
        from . import  sbo
        return sbo.SboElement(element).getSbo()

    def hasSetSbo(self,element):
        '''
        Function isSetSbo returns bool-value: has element set an SBO-Term?
        @return: has element SBO term?
        @rtype : bool
        '''
        from . import sbo
        return sbo.SboElement(element).hasSboTerm()
