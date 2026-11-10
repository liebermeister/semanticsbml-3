#!/usr/bin/env python

# import external modules
import libsbml, re, urllib.request, urllib.parse, urllib.error, difflib, sys, os

# import internal modules
from . import config, annotate
from libSBAnnotation import database_new



def checkSBO(section,option):
    '''
    check if SBO-flatfile exists
    input is the same as get for from config.Config
    @param section: the section in the config
    @type section: str
    @param option: the config entry
    @type option: str
    '''
    filepath = config.Config().getpath(section,option)
    if not os.path.exists(filepath):
        raise annotate.AnnotationDatabaseError(
            'No annotation database found (%s is missing). Build it with the command\n'
            '    semanticsbml build-database' % filepath)

checkSBO('resources','SBO_resources')

# define some lists, holding the ids, names, definitions, etc., of the Sbo Terms
idlist   = []
namelist = []
deflist  = []
isalist  = []
typelist = []
infolist = []
sbmlterm = None     # Testing purpose
db_dir   = None

class SboDatabase:

    __instance = None
    def __init__(self):
        if SboDatabase.__instance is None:
            SboDatabase.__instance = SboDatabase.__database()
        self.__dict__['_db_instance']= SboDatabase.__instance


    def __getattr__( self, attr ):
        return getattr(self.__instance, attr)

    def __setattr__( self, attr, value):
        return setattr( self.__instance, attr, value)


        
    class __database:
        '''
        Class __database is downloading the flatfile of Sbo including all Sbo Terms
        and further information. It builds 4 lists containing the ids, names,
        definitions and types of the Sbo Terms.
        '''
        
        def __init__(self):
            '''
            download SBO flatfile and build database
            '''
            global db_dir
            db_dir = config.Config().getpath('database', 'database_dir')
            
            try:
                data = open(db_dir+'SBO_OBO.obo','r')
            except IOError:
                hardfile     = open(db_dir+'SBO_OBO.obo','w')
                sboinput     = internetfile.read()
                hardfile.write(sboinput)
                hardfile.close()


            self.buildDatabase()

        def buildDatabase(self):
            '''
            Function buildDatabase is running over the downloaded Sbo flatfile and
            generates the lists idlist, namelist, deflist and typelist for further
            analyses.
            '''
            global idlist, namelist, deflist, isalist, typelist, db_dir

            x = 0
            f = open(db_dir+'SBO_OBO.obo','r')
              
            for i,line in enumerate(f):
                if line.startswith('[Term]'):
                    x = 1
                if line.startswith('id:'):
                    #idlist.append(line[4:len(line)].rstrip('\n'))
                    try:
                        idlist.append(int(line[8:len(line)].rstrip('\n')))
                    except:
                        idlist.append(line)
                if line.startswith('name:'):
                    namelist.append(line[6:len(line)].rstrip('\n'))
                if line.startswith('def:'):
                    deflist.append(line[5:len(line)].rstrip('\n'))
                if line.startswith('is_a:') and x == 1:
                    isalist.append(int(line[10:17]))
                    x = 2
                if line.startswith('is_a') and x == 2:
                    x = 0
                if line.startswith('is_o'):
                    idlist.pop()
                    namelist.pop()
                    deflist.pop()                                        
                if line.startswith('\n') and x!=0 and len(isalist)==len(idlist)-1:
                    isalist.append('')
                    x = 0
                if line.startswith('\n') and len(deflist)==len(idlist)-1:
                    deflist.append('')
                    x = 0

            #deflist.append('')
            isalist.append('')

            #print "id ", idlist[len(idlist)-5]
            #print "name ", namelist[len(namelist)-5]
            #print "def ", deflist[len(deflist)-5]
            #print "is a ", isalist[len(isalist)-5]

            for i,elem in enumerate(isalist):
                typ = self.get_type_by_isa(elem)
                if typ == '':
                    typ = self.get_type_by_isa(isalist[i])
                typelist.append(typ)

            #self.remove_obsoletes()
            self.make_fancy_deflist()
            self.make_fancy_namelist()

            self.id2name = dict( zip(idlist, namelist) )
            self.id2def  = dict( zip(idlist, deflist)  ) 
            self.id2type = dict( zip(idlist, typelist) )


        def make_fancy_deflist(self):
            '''
            Function make_fancy_deflist is removing all the ugly symbols from the definitions in the SBO-
            flatfile. Removed are '[]', '<*>', '\\' and the quotation marks at beginning and end of definition.
            Definitions that are too long to be handed over to a widget are cut down to a length of 200 and a
            remark is added to contact the SBO-Website for whole definition.
            '''
            

            global deflist,idlist

            for i,singleDef in enumerate(deflist):
                noMathMl     = deflist[i].split('<')[0]
                noBrackets   = noMathMl.split(' [')[0]
                noSlashes    = noBrackets.replace('\\','')
                noleftQuote  = noSlashes.lstrip('"')
                norightQuote = noleftQuote.rstrip('"')
                deflist[i]   = norightQuote

                if len(deflist[i]) > 100:
                        deflist[i] = deflist[i][:100] + '... (for further definition please see <a href="http://www.ebi.ac.uk/sbo/main/entry.do?sboId=SBO:%.7d">SBO-Website</a>)'%(idlist[i])


        def make_fancy_namelist(self):
            '''
            Function make_fancy_namelist is removing all backslashes that occur in the names of the SBO-flatfile.
            '''
            
            global namelist
            for i,name in enumerate(namelist):
                noSlashes = namelist[i].replace('\\','')
                namelist[i]   = noSlashes

        def remove_obsoletes(self):
            '''
            Function remove_obsoletes removes all obsolete SBO-Terms, e.g. "obsolete parameter, SBO:0000006"
            '''
            global namelist, deflist, idlist, isalist, typelist
            removelist = []

            for i,name in enumerate(namelist):
                if 'obsolete' in name:
                    removelist.append(i)

            removelist.reverse()

            print(removelist)
            del namelist[0]

            '''
            for i in removelist:
                namelist.remove[i]
                deflist.remove(i)                
                idlist.remove(i)                    
                isalist.remove(i)                    
                typelist.remove(i)
            '''

            
        def search_subtrees(self,sboterm):
            '''
            searches for every kind of inhibitors (gives back the subbranch of inhibitors)
            '''
            global isalist,idlist
            subtree_ids = [sboterm]

            for i, elem in enumerate(idlist):
                if isalist[i] == sboterm:
                        subtree_ids.append(idlist[i])
                        self.search_subtrees(idlist[i])

            return subtree_ids


        def get_name_for_id(self,element):
            '''
            Function get_name_for_id receives an id and returns the corresponding
            name.
            @param element: a Sbo identifier (e.g. 14)
            @type element : int
            @return       : Sbo Name
            @rtype        : string
            '''

            ident = element
            
            if type(element) == type('katze'):
                if element.startswith('SBO'):
                    ident = int(element[4:])
                else:
                    ident = int(element)
                    
            try:
                return self.id2name[ident]
            except:
                return 'unknown'

            return
                

        def get_definition_for_id(self,element):
            '''
            Function get_definition_for_id receives an id and returns the corresponding
            definition.
            @param element: a Sbo identifier (e.g. 14)
            @type element : int
            @return       : Sbo definition
            @rtype        : string
            '''
            global idlist, deflist

            if str(element).startswith('SBO'):
                term = int(element[4:])
            else:
                term = element

            for i,x in enumerate(idlist):
                if x == term:
                    return deflist[i]

            return 'unknown'

        def get_type_by_isa(self,element):
            '''
            Function get_type_by_isa returns the type for an entity. The parental
            types are predefined, all children are going over a recursive call of
            the function until they reach a parental type.
            @param element: an is_a relation
            @type element : string
            @return       : type
            @rtype        : string
            '''

            global isalist, idlist

            typ = ''

            if element == 235:
                typ = 'species'
            elif element == 2:
                typ = 'parameter'
            elif element == 231:
                typ = 'reaction'
            elif element == 64:
                typ = 'kinetic_law'
            else:
                for i,elem in enumerate(idlist):
                    if elem == element:
                        typ = self.get_type_by_isa(isalist[i])
            return typ

        def get_type_for_id(self,element):
            '''
            Function get_definition_for_id receives an id and returns the corresponding
            definition.
            @param element: a Sbo identifier (e.g. 14)
            @type element : int
            @return       : type of given Sbo id
            @rtype        : string
            '''
            global typelist

            ident = element

            if type(element) == type('katze'):
                if str(element).startswith('SBO'):
                    ident = int(element[4:])
                else:
                    ident = int(element)
                    
            try:
                return self.id2type[ident]
            except:
                return 'unknown'

        
        def search_id_by_name_exact(self,element):
            '''
            Function search_id_by_name_exact receives a name and returns a corresponding identifier.
            It performs an exact search for the name.
            @param element: a name (e.g. species name, reaction name, etc.)
            @type element : string
            @return       : id corresponding to the given name
            @rtype        : integer
            '''

            global idlist, namelist

            for i, x in enumerate(namelist):
                if element == str(x):
                    return str(idlist[i])

            return ''


class ModelSbo:
    '''
    Class for adding/modifying Sbo Terms of a libsbml element
    '''

    def __init__(self, libsbml_model=None):
        '''
        Reads the SBML file elements with Sbo Term
        @param libsbml_model: a libsbml-model taken from an xml-file
        @type libsbml_model : libsbml-model
        '''
        self._sboElements = []
        self._readElementsSbo(libsbml_model)
        ldb               = SboDatabase()

    def _readElementsSbo(self,libsbml_model):
        '''
        Goes through all elements with Sbo Term in libsbml model and saves
        them in list _sboElements.
        @param libsbml_model: a libsbml-model
        @type libsbml_model : libsbml-model
        '''

        global sbmlterm                            # Testing purpose
        #sbmlterm = libsbml.Reaction()              # Testing purpose

        for typ in 'Compartments Species Parameters Reactions Events Rules'.split():
            for e in getattr( libsbml_model, 'getListOf' + typ )():
                #if not (SboElement(e)).id == -1:
                self._sboElements.append(SboElement(e))
                #else:
                 #   pass
               
        #self._sboElements.append(sbmlterm)         # Testing purpose
        self._sboElements.reverse()
        self._sboElements.append(SboElement(libsbml_model))
        

    def getElementsSbo(self):
        '''
        returns a list of SBO elements in the model
        '''        
        return self._sboElements


    def remSboTerm(self,libsbml_element):
        '''
        Unsets given SboTerm by using unsetSboTerm-Function of SboElement()
        @param libsbml_model: a libsbml-model
        @type libsbml_model : libsbml-model
        '''
        self.libsbml_element = libsbml_element
        self.id              = libsbml_element.getSBOTerm()

        for i,elem in enumerate(self._sboElements):
            if elem == self.libsbml_element:
                UnsetTerm            = SboElement(elem)
                self._sboElements[i] = UnsetTerm.unsetSboTerm(elem)


    '''
    def addSboAutomatic(self):
        """ add SBO terms automatically (search xml values, id, name ...)"""
        # mal sehen ob wichtig ...
        raise Exception('Not implemented yet')
    '''



class SboElement:
    '''
    Class for adding SBO Terms to single elements
    '''

    def __init__(self, libsbml_element, sbo_term=None):
        '''
        Initializes the attributes for a given libsbml_element.
        @param libsbml_element: an element from a libsbml-model
        @type libsbml_element : libsbml-element
        @param sbo_term       : an optional sbo-term
        @type sbo_term        : string
        '''
        self._sbodb          = SboDatabase()
        self.libsbml_element = libsbml_element
        #sbml things we need to know (since rules get identical ids as their species, this leads to functional
        #problems. we therefore extend their id to "ID + [~ Rule]")     
        self.type            = self.libsbml_element.getElementName()
        if self.type == 'assignmentRule':
            self.libsbml_id = libsbml_element.getId() + ' [Assignment Rule]'
        elif self.type == 'rateRule':
            self.libsbml_id = libsbml_element.getId() + ' [Rate Rule]'  
        elif self.type == 'algebraicRule':
            self.libsbml_id = libsbml_element.getId() + ' [Algebraic Rule]'
        else:
            self.libsbml_id = libsbml_element.getId()   
        #self.libsbml_id      = self.libsbml_element.getId()
        self.getAnnotations = self.getSboTerm

        if sbo_term:
            self.setSbo(sbo_term)
        else:
            self.id          = libsbml_element.getSBOTerm()
            
        self.db              = 'SBO'
        self.metaid          = ''
        self._update()


    def _update(self):
        '''
        updates information for SBO element
        '''
        self.name        = self._sbodb.get_name_for_id(self.id)
        self.definition  = self._sbodb.get_definition_for_id(self.id)
        self.type        = self._sbodb.get_type_for_id(self.id)
        self.qualifier   = 'is'
       

    def checkMeOut(self):
        '''
        returns the information that an element is a SBO element (corresponding function in annotate.py for MIRIAM)
        '''        
        return 'SBO'

    def getQuerys(self):
        '''
        concatenate the name, id and metaid values of the libsbml element and return them
        @return: list of query strings
        @rtype: [str]
        '''
        return self.name
        #return [self.id,self.qualifier,self.name]

    def getName(self):
        '''
        Returns name of libsbml-element
        '''
        return self.name
    
    def getSboTerm(self):
        '''
        Returns Sbo term of libsbml-element
        '''
        if self.id == None or self.id == -1:
            return None
        else:
            return self.id

    def __str__(self):
        sbo = self.getSboTerm()
        if not sbo:
            return 'not set'
        else:
            return 'SBO:%.7d'%sbo


    def search_id_by_name_exact(self,element):
        '''
        as above function, only worse
        '''
        return self._sbodb.search_id_by_name_exact(element)

    def getDefinition(self):
        '''
        Returns definition of libsbml-element
        '''
        return self.definition

    def getType(self):
        '''
        Returns type of libsbml-element
        '''
        return self.type

    def getMoreGeneral(self,sboterm):
        '''
        receives SBO term and searches for a more general SBO term
        '''
        global idlist,typelist,isalist

        returnlist    = []
        ids_and_isas  = dict(zip(idlist,isalist))

        try:
            idMoreGeneral = ids_and_isas[sboterm]
        except:
            idMoreGeneral = ''
        
        if idMoreGeneral == '':
            returnlist.append('No more general SBO term found')
            returnlist.append('')
            returnlist.append('')
            returnlist.append('')
            return returnlist

        nameMoreGeneral = self._sbodb.get_name_for_id(idMoreGeneral)
        typeMoreGeneral = self._sbodb.get_type_for_id(idMoreGeneral)
        defMoreGeneral  = self._sbodb.get_definition_for_id(idMoreGeneral)

        returnlist.append(idMoreGeneral)
        returnlist.append(nameMoreGeneral)
        returnlist.append(typeMoreGeneral)
        returnlist.append(defMoreGeneral)        

        return returnlist

    def getMoreSpecific(self,sboterm):
        '''
        receives SBO term and searches for more specific SBO terms
        '''
        global idlist,namelist,typelist,deflist,isalist
        itemlist   = []
        returnlist = []

        for i,ident in enumerate(isalist):
            if ident == sboterm:
                itemlist.append(idlist[i])
                itemlist.append(namelist[i])
                itemlist.append(typelist[i])
                itemlist.append(deflist[i])
                returnlist.append(itemlist)
            itemlist = []

        if returnlist == []:
            returnlist.append('No more specific SBO terms found')

        return returnlist

    def unsetSboTerm(self):
        '''
        Unsets the attributes of a given libsbml-element
        @param libsbml_element : a libsbml-element
        @type libsbml_element  : libsbml-element
        @return libsbml_element: libsbml-element with unset attributes
        @rtype libsbml_element : libsbml-element
        '''
        self.id     = None
        self.libsbml_element.unsetSBOTerm()
        self.metaid = ''
        self._update()


    def setSbo(self,sboterm):
        '''
        Unsets the attributes of a given libsbml-element
        @param libsbml_element : a libsbml-element
        @type libsbml_element  : libsbml-element
        @param sboterm         : an Sbo term
        @type sboterm          : string
        @return libsbml_element: libsbml-element with new attributes
        @rtype libsbml_element : libsbml-element
        '''

        if type(sboterm) != type(1):
            raise Exception('Wrong type, required int') 

        #print sboterm

        #if not self.name == 'unknown':
        self.id     = sboterm
        self.libsbml_element.setSBOTerm(sboterm)
        self.metaid = 'SBO'
        self._update()

    def findByKeyword(self,name,fuzzyness):
        '''
        Function findByKeyword receives a keyword and starts a search over all names, whether the keyword appears in it.
        @param name : a keyword
        @type name  : string
        @return infolist : list of matching names (including id, type, definition)
        @rtype infolist : list of 4-tuple-lists
        '''

        global idlist, namelist, typelist, deflist
        infolist      = []
        doublechecker = []
 
        if type(name) != type('name'):
            name = str(name)
        try:
                name_parts = name.lower().rsplit(' ')
        except:
                name_parts = name

    
        for i, x in enumerate(namelist):
            singlewords = x.lower().rsplit(' ')
            for name_part in name_parts:
                    infos = difflib.get_close_matches(name_part.lower(),singlewords,2,fuzzyness)
            if infos != []:
                infolist.append([idlist[i], namelist[i], typelist[i], deflist[i]])
                doublechecker.append(idlist[i])
                
                
        for i, x in enumerate(deflist):
            singlewords = x.lower().rsplit(' ')
            infos = difflib.get_close_matches(name.lower(),singlewords,2,fuzzyness)
            if infos != []:
                if not idlist[i] in doublechecker:
                    infolist.append([idlist[i], namelist[i], typelist[i], deflist[i]])

        if infolist == []:
            infolist.append('No suggestions found in SBO database')

        return infolist

    def isAnnotated(self):
        '''
        Returns 'True' if the element has SBO Term, else otherwise
        @return: has SBO Term
        @rtype : bool
        '''
        return self.libsbml_element.isSetSBOTerm()


if __name__ == "__main__":

     # TESTING Sbo
     global f
     global ldb
     ldb = SboDatabase()                          # generate instance of class SboDatabase
     #print ldb.get_name_for_id(1)                # insert id, get name
     #print ldb.getSuggestions('events')          # insert exact name, get id+name+definition
     #print ldb.getSuggestions('obslete')         # insert approximate name, get id+name+definition (fuzzy)
     #print ldb.get_type_by_isa('SBO:0000196')


     # TESTING ModelSbo
     document  = libsbml.readSBMLFromString(open('/home/tlubitz/Desktop/sboclass/semanticsbml/trunk/models/sbo.xml','r').read())
     sbmlmodel = document.getModel()
     model     = ModelSbo(sbmlmodel)         # generate model
     #model.remSboTerm(sbmlterm)             # removes Sbo term of given sbmlterm
                                             # sbmlterm is no part of actual file: libsbml.Reaction() --> generated above
     #print model.getElementsSbo()

     # TESTING SboElement
     #Element = SboElement(sbmlterm)          # generate element
     #print Element.getSbo()                  # set a given sbmlterm to sboterm 14 (SBO:0000014) plus name, def, etc.

# TODO: If we change the getSBOTerm-Output (e.g. 14) AT ONCE to the rightly defined 'SBO:0000014', then we can omit the triple or quadruple new definition in the different functions
