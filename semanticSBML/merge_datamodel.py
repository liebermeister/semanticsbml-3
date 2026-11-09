'''
Datamodel used for the merging
'''
#//////////////////////////////////////////////////////////////////////////////////
#external libaries
import libsbml
import sys 
#semanticSBML libraries
import libSBAnnotation.config
from . import annotate
from . import merge_logic
from . import sbo
from . import xmlobject

#//////////////////////////////////////////////////////////////////////////////////
class MergeError(Exception):

    def __init__(self,message):
        self.errorMessage = message

    def __repr__(self):
        return self.errorMessage

    def __str__(self):
            return self.errorMessage

class MergeBioQuantityLocationError(MergeError):
    def __init__(self,message):
        self.errorMessage = 'Error in Merge_BioQuantity: %s'% message
#//////////////////////////////////////////////////////////////////////////////////
class MergeBase:
    '''
    The MergeBase provides basic values that are needed to handle libsbml elements
        - the libsbml element itself
        - the model representation (L{OriginModel}) belonging to the element
        - the element type e.g. reaction
        - a string representation of the element (element name or id)
    '''
    def __init__(self,libsbml_element,merge_model):
        self.libsbml_element    = libsbml_element
        if not merge_model:
            self.libsbml_model = None
        else:
            self.libsbml_model = merge_model.libsbml_model
        self.type        = libsbml_element.getElementName()
        self.repr        = libsbml_element.getName() or libsbml_element.getId()
        self.merge_model = merge_model
#//////////////////////////////////////////////////////////////////////////////////
class SubBioentity():
    '''
    The external identity of a L{OriginElement} is defined by this class.
    It stores MIRIAM annotations that links the element to external databases.
    Functions for getting and removing L{Annotation}s are provided. 
    '''
    def __init__(self,libsbml_element): 
        #instance variable
        self.miriam_annotation_element=annotate.ElementAnnotations(libsbml_element,True)
        #class methods
        self.remAnnotation=self.miriam_annotation_element.remAnnotation#remove a single L{Annotation} from the element
        self.getMiriamAnnotations=self.miriam_annotation_element.getAnnotations
        self.unsetAnnotations=self.miriam_annotation_element.unsetAnnotations
        self.addAnnotation=self.miriam_annotation_element.addAnnotation
        
        self.sbo_annotation_element = sbo.SboElement(libsbml_element)
        self.getSboTerm = self.sbo_annotation_element.getSboTerm
        self.setSboAnnotation = self.sbo_annotation_element.setSbo

#//////////////////////////////////////////////////////////////////////////////////
class SubBioquantity(MergeBase):
    def __init__(self,libsbml_element,merge_model):
        MergeBase.__init__(self,libsbml_element,merge_model)
        #class method
        self.searchUnit= lambda x: self.merge_model.global_unit_defs[x]
    def getBioEntities(self):
        raise NotImplementedError
    def setBioEntities(self,entities):
        pass
    def getQuantityType(self):
        return ''
    def setQuantityType(self,type):
        pass
    def getUnit(self,refresh=False):
        raise NotImplementedError
    def setUnit(self,unitdef):
        if self.type != 'reaction':#kinda destroys the sense of this interface but its more visual here
            self.unit=unitdef
    def getLocation(self):
        raise NotImplementedError
    def setLocation(self,location):
        self.location=location
#//////////////////////////////////////////////////////////////////////////////////
class SubBioquantitySpecies(SubBioquantity):
    ####################################
    def getQuantityType(self):
        #Sarah Keating said only this is needed
        if self.libsbml_element.getHasOnlySubstanceUnits():
            return 'Amount'
        else:#it is a conentration  
            return 'Concentration'
    ####################################
    def setQuantityType(self,type):
        if type=='Concentration':
            self.libsbml_element.setHasOnlySubstanceUnits(False)
        if type =='Amount':
            self.libsbml_element.setHasOnlySubstanceUnits(True)
    ####################################
    def getUnit(self):
        if self.getQuantityType == 'Concentration': 
            substance=self.merge_model.global_unit_defs['substance']
            su=self.libsbml_element.getSubstanceUnits()
            if su:
                substance   = self.searchUnit(su)
            volume = self.location.getUnit(True)#TODO TODO TODO this is not a solution, why can we not deepcopy the unit??
            vo=self.libsbml_element.getSpatialSizeUnits()
            if vo:
                volume      = self.searchUnit(vo)
            
            return UnitDef(unit_defs=(substance,volume))
        else:
            units=self.libsbml_element.getUnits()
            if units:
                return self.searchUnit(units) 
            else:
                return self.merge_model.global_unit_defs['substance']
    ####################################
    def getLocation(self):
        return self.merge_model.id2merge_element[self.libsbml_element.getCompartment()]
            
#-----------------------------------------------
class SubBioquantityReaction(SubBioquantity):
    def __init__(self,libsbml_element,merge_model):
        SubBioquantity.__init__(self,libsbml_element,merge_model)
        self.entity_types=['Reactant','Modifier','Product']

    ####################################
    def getBioEntities(self):#good
        '''
        @rtype entities: {'Reactant':[],'Modifier':[],'Product':[]}
        '''
        ret={'Reactant':[],'Modifier':[],'Product':[]}
        #try:
        for part in self.entity_types:
            for i in range(getattr(self.libsbml_element,'getNum%ss'%part)()):
                ret[part].append(self.merge_model.id2merge_element[getattr(self.libsbml_element,'get%s'%part)(i).getSpecies()])
        #except KeyError,e:
        #   tmp=self.merge_model.id2merge_element.keys()
        #   tmp.sort()
        #   print 'grsslsl ',tmp,e
        #   sdlfk
        return ret

    ####################################
    def setBioEntities(self,entities):#TODO check
        '''
        @type entities: {'Reactant':[],'Modifier':[],'Product':[]}
        '''
        for type in self.entity_types:
            list_of = getattr(self.libsbml_element,'getListOf'+type+'s')()
            list_of.clear()
            for me in entities[type]:
                species_reference=getattr(self.libsbml_element,'create'+type)()
                if me.getTuple():
                    species_reference.setSpecies(me.getTuple().getMergedElement().getId())
                else:
                    species_reference.setSpecies(me.getId())

    ####################################
    def getUnit(self):
        time=self.merge_model.global_unit_defs['time']
        return UnitDef(unit_defs=(self.merge_model.global_unit_defs['substance'],time))

    ####################################
    def getLocation(self):  
        res=[]
        #tmp=self.merge_model.id2merge_element.keys().sort()
        #tmp.sort()
        #print 'id2merge_element ',tmp
        for type in self.entity_types:
            tmp_res=set()
            if getattr(self.libsbml_element,'getNum'+type+'s')()>0:
                for n in range(getattr(self.libsbml_element,'getNum'+type+'s')()):
                    tmp_res.add(self.merge_model.id2merge_element[getattr(self.libsbml_element,'get'+type)(n).getSpecies()].getLocation())
            res.append(list(tmp_res))
        return res 

#------------------------------------------------               
class SubBioquantityCompartment(SubBioquantity):
    #####
    def __init__(self,libsbml_element,merge_model):
        SubBioquantity.__init__(self,libsbml_element,merge_model)
        self.qtypes=['length', 'area', 'volume']
    
    ####################################
    def getQuantityType(self):
        dim2type=dict(zip(list(range(1,4)),self.qtypes))
        try:
            return dim2type[self.libsbml_element.getSpatialDimensions()]
        except KeyError:
            return ''
    ####################################
    def setQuantityType(self,type):
        type2dim=dict(zip(self.qtypes,list(range(1,4))))
        try:
            self.libsbml_element.setSpatialDimensions(type2dim[type])
        except KeyError:
            raise MergeError('Unknown Quantity Type')
    ####################################
    def getUnit(self):
        units=self.libsbml_element.getUnits() 
        if not units: 
            dim2type=dict(zip(list(range(1,4)),self.qtypes))
            spacialDimension=self.libsbml_element.getSpatialDimensions()
            return self.merge_model.global_unit_defs[dim2type[spacialDimension]]
        else:
            return self.searchUnit(units) 
    ####################################
    def getLocation(self):
        if self.libsbml_element.isSetOutside():
            if self.libsbml_element.getOutside():
                return self.merge_model.id2merge_element[self.libsbml_element.getOutside()]
            else:
                raise MergeBioQuantityLocationError('requeue location since it is linking to a yet unknown location')
#------------------------------------------------               
class SubBioquantityParameter(SubBioquantity):
    ####################################
    def getUnit(self):
        unit=self.libsbml_element.getUnits()
        if not unit:
            #sys.stderr.write( 'Warning: parameter %s unit not set!'%self.libsbml_element.getId())
            return None 
        else:
            return self.searchUnit(unit)
#//////////////////////////////////////////////////////////////////////////////////
class SubModelstatement(MergeBase):
    def __init__(self,libsbml_element,merge_model):
        MergeBase.__init__(self,libsbml_element,merge_model)
        self.statement=None
    def getFast(self):
        raise NotImplementedError
    def setFast(self,val):
        raise NotImplementedError
    def getConstant(self):
        if self.type in ['species','parameter','compartment']:
            return self.libsbml_element.getConstant()
        else:
            raise NotImplementedError
    def setConstant(self,val):
        if self.type in ['species','parameter','compartment']:
            if isinstance(val,str):
                if val=='True': val=True 
                elif val=='False': val=False 
                else: raise MergeError(' setConstant failed because % was inserted'%val)
            if not isinstance(val, bool):val=bool(val)
            self.libsbml_element.setConstant(val)
            self.constant=val
        else:
            raise NotImplementedError

    def getCharge(self):
        raise NotImplementedError
    def setCharge(self,val):
        raise NotImplementedError
    def getReversible(self):
        raise NotImplementedError
    def setReversible(self,val):
        raise NotImplementedError
    def getSize(self):
        raise NotImplementedError
    def setSize(self,val):
        raise NotImplementedError
    def getQuantity(self):
        raise NotImplementedError
    def setQuantity(self,val):
        raise NotImplementedError
    def getStatement(self):
        '''
        overwritten by reaction subclass
        '''
        if self.statement:
            return self.statement#must return this since only later on statements will be collected
        math_obj=None
        if self.libsbml_element.getId() in self.merge_model.id2rule:
            #if self.statement:
            #   sys.stderr.wirte('ooops big error, I wasnt sure this will happen, more than one rule for this element',self.libsbml_element.getId())
            math_obj = self.merge_model.id2rule[self.libsbml_element.getId()]
        
        #search for attached function definitions
        if math_obj:
            self._attatchFunctionDefinition(math_obj)
            return math_obj

    def setStatement(self,math_obj):
        self.statement = math_obj
        
    def getInitialAssignment(self):
        if self.libsbml_element.getId() in self.merge_model.id2assignment:
            ia=self.merge_model.id2assignment[self.libsbml_element.getId()]
            self._attatchFunctionDefinition(ia)
            return ia

    def setInitialAssignment(self,math_obj):
        self.initialAssignment=math_obj

    def _attatchFunctionDefinition(self,math_obj):
        for en in merge_logic.BioRelations().getMathElementNames(math_obj.getMath()):
            if en in self.merge_model.id2function_definition:
                try:
                    math_obj.functionDefinitions.append(self.merge_model.id2function_definition[en].clone())
                except AttributeError as e:
                    math_obj.functionDefinitions=[self.merge_model.id2function_definition[en].clone()]
#//////////////////////////////////////////////////////////////////////////////////
class SubModelstatementSpecies(SubModelstatement):
    ####################################
    def getCharge(self):
        return self.libsbml_element.getCharge()
    ####################################
    def setCharge(self,val):
        if not isinstance(val,float):val=float(val)
        return self.libsbml_element.setCharge(val)
        self.charge=val
    ####################################
    def getQuantity(self):
        if self.libsbml_element.isSetInitialConcentration():
            return ('Concentration',self.libsbml_element.getInitialConcentration())
        if self.libsbml_element.isSetInitialAmount():
            return ('Amount',self.libsbml_element.getInitialAmount())
    ####################################
    def setQuantity(self,val):
        if not isinstance(val,tuple):
            raise OriginModelStatementError('ModelStatement quantity for a species must be a tuple (type,value)')
        type=val[0]
        value=val[1]
        if not isinstance(value, float):value=float(value)
        if type=='Concentration':
            self.libsbml_element.setInitialConcentration(value)
        elif type=='Amount':
            self.libsbml_element.setInitialAmount(value)
        else:
            nononothisisnotright
        self.quantity=val

#------------------------------------------------               
class SubModelstatementReaction(SubModelstatement):
    ####################################
    def getFast(self):
        return self.libsbml_element.getFast()
    ####################################
    def setFast(self,val):
        if not isinstance(val, bool):val=bool(val)
        self.libsbml_element.setFast(val)
        self.fast=val

    ####################################
    def getReversible(self):
        return self.libsbml_element.getReversible()
    ####################################
    def setReversible(self,val):
        if not isinstance(val, bool):val=bool(val)
        self.libsbml_element.setReversible(val)
        self.reversible=val

    ####################################
    def getStatement(self):
        math_obj = self.libsbml_element.getKineticLaw()
        #search for attached function definitions
        if math_obj:
            self._attatchFunctionDefinition(math_obj)
            return math_obj
    ####################################
    def setStatement(self,math_obj):
        #get the correct rule or kinetic law from other model
        self.statement = math_obj
        self.libsbml_element.setKineticLaw(math_obj)

#------------------------------------------------               
class SubModelstatementCompartment(SubModelstatement):
    ####################################
    def getSize(self):
        if self.libsbml_element.isSetSize():
            return self.libsbml_element.getSize()
    ####################################
    def setSize(self,val):
        if not isinstance(val, float):val=float(val)
        self.libsbml_element.setSize(val)

            
#------------------------------------------------               
class SubModelstatementParameter(SubModelstatement):
    ####################################
    def getQuantity(self):
        if self.libsbml_element.isSetValue():
            return ('Value',self.libsbml_element.getValue())
    ####################################
    def setQuantity(self,val):
        if not isinstance(val, float):val=float(val)
        return self.libsbml_element.setSize(val[1])

#//////////////////////////////////////////////////////////////////////////////////
class OriginElement(MergeBase):
    '''
    abstraction of a biological entity for merging

    wrap class that consists of three sublcass of the base type
        - SubBioentity
        - SubBioquantity
        - SubModelstatement
    '''
    def __init__(self,libsbml_element,merge_model):

        #.....................................

        MergeBase.__init__(self,libsbml_element,merge_model)
        merge_model.libsbml_model.getSBMLDocument().setLevelAndVersion(2,4,False)
        #.....................................
        #libsbml metadata in every element TODO kill stuff not needed
        self.name   = libsbml_element.getName()
        self.id     = libsbml_element.getId()
        self.newid  = self.id
        #global unique idenfier
        self.mid=self.id
        #self.mid='%s-%s'%(self.id,str(self.merge_model.mid))
        #.....................................
        #link to OriginElements that show a similarity
        self._merge_tuple_candidates=set()
        self.tuple=None
        #.....................................
        #signal for messgages to outside funcitons  
        self.sig_tuple_changed=None
        #.....................................
        #defined by 
        #use command dispatcher pattern to get the right subclass
        self._bioentity         = SubBioentity(self.libsbml_element)
        self._bioquantity   = globals()['SubBioquantity%s'%self.type.title()](self.libsbml_element,merge_model)
        self._modelstatement    = globals()['SubModelstatement%s'%self.type.title()](self.libsbml_element,merge_model)
        #create a dictionary of all attributes so they can be directly accessed
        self._allattributes = dict([(x,self._bioentity) for x in dir(self._bioentity) if not '_' == x[0]]+\
                 [(x,self._bioquantity) for x in dir(self._bioquantity) if not '_' in x]+\
                 [(x,self._modelstatement) for x in dir(self._modelstatement) if not '_' in x])
        #.....................................


    def __getattr__(self, n):   
        '''
        realize wrap pattern!
        make methods of the sub objects
            - SubBioquantity
            - SubBioentity
            - SubModelstatement
        directly accessible
        '''
        if n in self._allattributes:
            return getattr(self._allattributes[n], n)
        else:
            raise AttributeError(n)

    def getDepends_onElements(self):
        return self.merge_model.getElementsDependsOn(self)

    def getDependingElements(self):
        return self.merge_model.getDependingElements(self)

    def connect_view_functions(self,tuple_changed_fkt):
        '''
        set an external function that should be executed when the tuple this element belongs to is changed
        @type tuple_changed_fkt: function pointer
        '''
        #print 'e: connect view fkt'#DEBUG
        self.sig_tuple_changed=tuple_changed_fkt

    def addTupleCandidate(self,element,force=False):
        '''
        add and element to the list of similar elements
        @type element: OriginElement
        '''
        neighbous=[]
        if self.tuple and not force:neighbous=self.tuple.getElements()
        if self.merge_model!=element.merge_model and not (element in neighbous):
            self._merge_tuple_candidates.add(element)

    def addTupleCandidates(self,elements):
        '''
        update the list of similar elements
        @type elements: [OriginElement]
        '''
        for e in elements:
            self.addTupleCandidate(e)

    def getTupleCanidates(self):
        '''
        returns a list of similar elements, elements that have matching MIRIAM annotations
        @rtype: [OriginElement]
        '''
        return list(self._merge_tuple_candidates)
    def getTuple(self):
        '''
        returns the L{MergeTuple} the element instance belongs to, if it does not belong to an element None is returnend
        @rtype: MergeTuple or None
        '''
        return self.tuple
    def setTuple(self,tuple):
        '''
        set the L{MergeTuple} the element belongs to
        if an external function of a view is set, call the function so the view can be updated
        @type tuple: MergeTuple
        '''
        #global vars
        #print 'tuple set',tuple,self.mid#DEBUG
        old_tuple=self.tuple
        self.tuple=tuple
        if self.sig_tuple_changed and self.tuple!=old_tuple:
            #TODO bad hack we need a global switch here to check wich merge phase we are!!!!!!!!!!!!!!!!!
            #add elements to update queue
            merge_logic.BioRelations().elements_update_queue=merge_logic.BioRelations().elements_update_queue.union(self.getDependingElements())
            #update the external view function
            #print 'sending signal tuple changed'#DEBUG
            self.sig_tuple_changed(self)
        #print self.mid,'tuple set:',self.getTuple()
    def setId(self,id):
        '''
        set the "new" id, this identifier is globally unique in the merged model
        @type id: str
        '''
        self.newid=id
    def getId(self):
        return self.newid
#//////////////////////////////////////////////////////////////////////////////////
#class MergedElement(OriginElement):
#   pass
#//////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////
#class MergeTuple:
#   pass
#//////////////////////////////////////////////////////////////////////////////////
class OriginModel:
    '''
    abstraction of a model for merging
    a model contains biological entities that can be merged
    '''
    def __init__(self,libsbml_model,mid):
        '''
        '''
        #--------------------------------------
        self.libsbml_model = libsbml_model
        self.mid           = mid
        self.repr          = self.libsbml_model.getName() or self.libsbml_model.getId()
        #--------------------------------------
        self.global_unit_defs={}
        self._getUnits()
        #--------------------------------------
        #get assignment/rate rules
        self.id2rule = self._getRules()
        #--------------------------------------
        #get initital assignments
        self.id2assignment={}
        for asm in list(self.libsbml_model.getListOfInitialAssignments()):
            self.id2assignment[asm.getId()]=asm
        #--------------------------------------
        self.id2function_definition={}
        for fd in list(self.libsbml_model.getListOfFunctionDefinitions()):
            self.id2function_definition[fd.getId()]=fd
        self._buildElements()

    def _buildElements(self):
        #////////////////////////////////////////////////////////////////////////////
        #////////////////////////////////////////////////////////////////////////////
        #--------------------------------------
        self._merge_elements=[]
        #process compartments first since they are neede for units etc
        compartment_queue=list(self.libsbml_model.getListOfCompartments())
        while compartment_queue:
            element = compartment_queue.pop(0)
            try:
                self._merge_elements.append(OriginElement(element,self))
            except MergeBioQuantityLocationError:
                compartment_queue.append(element)
        #get all merge elements
        for listname in ("Species", "Parameters", "Reactions"):#TODO must extend the list (), outsource to XML file
            for element in getattr(self.libsbml_model, 'getListOf' + listname)():
                self._merge_elements.append(OriginElement(element,self))
        #--------------------------------------
        self.id2merge_element=dict([(x.id, x) for x in self._merge_elements])
        #print 'id2merge_element: ',self.id2merge_element.keys()
        #--------------------------------------
        #build mapping of libsmbl-element-type to merge-elements
        self._type2merge_elements={'species':[],'parameter':[],'reaction':[],'compartment':[]}
        for me in self._merge_elements:
            self._type2merge_elements[me.type].append(me)
        #--------------------------------------
        #build element dependency dicts
        self.element_depends_on = {}
        self.elements_depending_on = {}

        for e in self.getElements():
            self.element_depends_on[e] = self._dependingElements(e)
            for de in self.element_depends_on[e]:
                if de in self.elements_depending_on: self.elements_depending_on[de].add(e)
                else: self.elements_depending_on[de]=set([e])
        for e in [x for x in self.getElements() if not x in list(self.elements_depending_on.keys())]:
            self.elements_depending_on[e]=set([])

    def getElements(self):
        return self._merge_elements

    def getElementsByType(self,type):
        type2compare_types = {#TODO write this into an xml file
                'species':['species','parameter'],
                'parameter':['parameter','species','compartment'],
                'reaction':['reaction'],
                'compartment':['compartment','parameter'],
                }
        res=[]
        for t in type2compare_types[type]:
            res+=self._type2merge_elements[t]
        return res
    
    def _getRules(self):
        out={}
        for rule in list(self.libsbml_model.getListOfRules()):
            if isinstance(rule, libsbml.AssignmentRule) or isinstance(rule, libsbml.RateRule):
                out[rule.getVariable()]=rule
            elif isinstance(rule, libsbml.AlgebraicRule):#TODO wirte this in merger as check at the beginning
                raise MergeError('AlgebraicRules are currently not supported.')
        return out
    
    def _getUnits(self):
        '''
        build a list of units defined in the model
        complete list with standard units
        this funtion takes no input and creates no output but fills self.global_unit_defs
        '''
        #----------------------------------------
        #look up units in the model
        for ud in list(self.libsbml_model.getListOfUnitDefinitions()):
            self.global_unit_defs[ud.getId()]=UnitDef(ud)
        #----------------------------------------
        #get standard units, and fill if not already set, check for errors
        x = xmlobject.XMLFile(path=libSBAnnotation.config.Config().getpath('resources','libsbml_resources'))
        for unit in x.root.DefaultUnitDefinitions[0].unit:
            if unit.id not in self.global_unit_defs:
                self.global_unit_defs[unit.id]=UnitDef(default_unit_def_name=unit.id)
        for unit in x.root.knownUnits[0].kinds.split():
            if unit not in self.global_unit_defs:
                self.global_unit_defs[unit]=UnitDef(default_unit_def_name=unit)

    def _dependingElements(self,ielement):
        '''
        get a list of elements this element depends on
        '''
        #print 'get dependin on',self.type,self.mid#DEBUG

        out=[]
        #--------------------------------
        #get location elements
        if ielement.type == 'reaction':
            for els in ielement.getLocation()+list(ielement.getBioEntities().values()):#get reaction participants and their locations
                for e in els: 
                    if e: out.append(e)
        else:
            try:
                if ielement.getLocation(): out.append(ielement.getLocation())
            except NotImplementedError:
                pass
        #--------------------------------
        #get elements in mathematical statements
        if ielement.getStatement():
            out += [self.id2merge_element[n] for n in merge_logic.BioRelations().getMathElementNames(ielement.getStatement().getMath()) if n in self.id2merge_element]

        return set(out)
        #print 'a',[e.mid for e in self._depends_on]#DEBUG

    def getDependingElements(self,element):
        '''
        get a list of elements depending on the inserted element
        @rtype: set(L{OriginElement}) 
        '''
        return [self.id2merge_element[ret] for ret in self.elements_depending_on[element] if isinstance(ret,str)]+[ret for ret in self.elements_depending_on[element] if not isinstance(ret,str)]

    def getElementsDependsOn(self,element):
        '''
        get e list of elements the inserted element depends on
        @rtype: set(L{OriginElement})
        '''
        return [self.id2merge_element[ret] for ret in self.element_depends_on[element] if isinstance(ret,str)]+[ret for ret in self.element_depends_on[element] if not isinstance(ret,str)]
#//////////////////////////////////////////////////////////////////////////////////

class MergedModel(OriginModel):
    def __init__(self,libsbml_model,origin_models):
        OriginModel.__init__(self,libsbml_model,'M')
        self.id2merge_element={}
        self.global_unit_defs={}
        self.id2rule={}
        for m in origin_models: #did,doc in enumerate(origin_models):
            #model = doc.getModel()
            #m     = OriginModel(model,did)
            self.id2merge_element.update(m.id2merge_element)
            self.global_unit_defs.update(m.global_unit_defs)
            self.id2rule.update(m.id2rule)

    def refresh(self):
        OriginModel.__init__(self,self.libsbml_model,'M')
    def _buildElements(self):
        pass
    def addElement(self,element):
        self.id2merge_element[element.getId()]=element
    def remElement(self,element):
        self.id2merge_element.pop(element.getId())
#//////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////



#//////////////////////////////////////////////////////////////////////////////////
class UnitDef:
    '''
    unit definition
    can contain multiple L{Unit}s
    '''
    def __init__(self,libsbml_unit_definition=None,default_unit_def_name='',unit_defs=[]):
        '''
        either a libsbml.UnitDefinition is inserted or the name of a default unit, 
        the default units are looked up in the libsbmlresource

        @type libsbml_unit_definition: L{libsbml.UnitDefinition}
        @type default_unit_def_name: str
        @type unit_items= [Unit]
        '''
        self.id=''
        self.unit_items=[]
        self.unit_defs=()
        self.type=''
        self.libsbml_ud=None

        if libsbml_unit_definition:#construct unit by libsbml.UnitDefinition
            self.type='libsbml'
            self.libsbml_ud=libsbml_unit_definition
            self.id=libsbml_unit_definition.getId()
            for unit_item in list(libsbml_unit_definition.getListOfUnits()):
                self.unit_items.append(Unit(libsbml.UnitKind_toString(unit_item.getKind()),unit_item.getScale() or 0.0 ,unit_item.getMultiplier() or 1.0,unit_item.getExponent()))
        elif default_unit_def_name:#construct unit by default unit def name
            self.type='default'
            x = xmlobject.XMLFile(path=libSBAnnotation.config.Config().getpath('resources','libsbml_resources'))
            for xml_unit in x.root.DefaultUnitDefinitions[0].unit:#predefined units
                if xml_unit.id==default_unit_def_name:
                    try:
                        exponent=xml_unit.exponent
                    except:
                        exponent=1.0
                    self.unit_items.append(Unit(kind=xml_unit.kind,exponent=exponent))
                    break
            else:#predefined kinds
                self.unit_items.append(Unit(default_unit_def_name))
        elif unit_defs:#contruct unit by
            self.type='composed'
            self.unit_defs=unit_defs
        else:
            raise MergeError('UnitDef was not properly initialized')
        
    def toHTML(self):
        '''
        return a html string representation of the unit definition
        '''
        seperator='*'
        res=''

        if self.unit_defs and self.unit_items:
            return self.unit_defs[0].toHTML() + '*' + self.unit_defs[1].unit_items[0].toHTML(-1.0)
        else:
            for unit in self.unit_items:
                res+=unit.toHTML()+seperator
            return res[:-1]

    def toStr(self):
        return str(self)
    #@toStr
    def __str__(self):
        '''
        return a string representation of the unit
        '''
        res=''
        seperator=' * '
        if self.unit_defs:
            return self.unit_defs[0].toStr()+' * '+self.unit_defs[1].unit_items[0].toStr(-1.0)
        else:
            for unit in self.unit_items:
                res +='%s%s'%(unit.toStr(),seperator)
            return res[:-3]

#//////////////////////////////////////////////////////////////////////////////////
class Unit:
    def __init__(self,kind,scale=0.0,multiplier=1.0,exponent=1.0):
        known_units=xmlobject.XMLFile(path=libSBAnnotation.config.Config().getpath('resources','libsbml_resources')).root.knownUnits[0].kinds.split()
        if kind in known_units:
            self.kind=kind
        else:
            raise MergeError('Unknown unit kind \'%s\' '%kind)
        self.scale=float(scale)
        self.multiplier=float(multiplier)
        self.exponent=exponent
    def toStr(self,exponent_mulitiplier=1.0):
        return 'Kind: %s | Scale: %s | Mulitiplier: %s | Exponent %s'%(self.kind,self.scale,self.multiplier,self.exponent*exponent_mulitiplier)
    def toHTML(self,exponent_mulitiplier=1.0):
        res=''
        if self.multiplier!=1.0:
            res+='%s*'%self.multiplier
        res+=self.kind
        if self.scale!=0.0:
            res+='*10<sup>%s</sup>'%self.scale
        if self.exponent!=1.0 or exponent_mulitiplier!=1:
            res='(%s)<sup>%s</sup>'%(res,self.exponent*exponent_mulitiplier)
        return res 
