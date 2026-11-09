#

'''
B{Merging of Multiple SBML Models}
this module enables the merging of sbml models
to merge two models an internal representation of the model is created
the elements of the internal representatin are compared
a list of matching elements is created
interfaces for modyfieing the match list are provided
matching elements can be in conflict
interfaces for conflict resolution are provided
if there are no conflicts or all conflicts are resolved the internal representation of the model is translated back into a new sbml model

'''


#python libraries
import copy, re, os, sys
from . import xmlobject
#external libaries
import libsbml 
#semanticSBML libraries
import libSBAnnotation.database
import libSBAnnotation.config
import semanticSBML.annotate
from . import Model
from . import sbmlmisc
from .merge_datamodel import OriginElement, OriginModel, MergedModel, MergeError, MergeBioQuantityLocationError
from .merge_logic import BioRelations 


#/////////////////////////////////////////////
#GLOBAL VARIABLES
setStatusMessage=None#module function to set status messages of an external view
#/////////////////////////////////////////////

g_known_units = None
g_bio_qualifier_list=[]
def _getLibsbmlResouces():
    global g_known_units
    try:
        g_known_units = xmlobject.XMLFile(path=libSBAnnotation.config.Config().getpath('resources','libsbml_resources')).root.knownUnits[0].kinds.split()
    except:
        sys.stderr.write('error opening file (missing?) ',libSBAnnotation.config.Config().getpath('resources','libsbml_resources'))
        import sys
        sys.exit(1)
    else:
        bq = xmlobject.XMLFile(path=libSBAnnotation.config.Config().getpath('resources','libsbml_resources')).root.BioQualifier
        for q in bq[0].qualifier:
            g_bio_qualifier_list.append(q.name) 

class MergeTupleError(MergeError):
        def __init__(self,message):
                self.errorMessage = 'Error in element tuple: %s'% message

class MergedElement(OriginElement):
    '''
    the merged elements represent the combined information of a tuple-element
    it is derived from a L{OriginElement} and thus contains the same information
    it also has functions that create dictionaries that contain conflicting values of the tuple-elements
    all values of this element will have updated ids if the prepare function was called, that is ids
    that only point to merged elements and elements that will not be merged
    '''
    def __init__(self,tuple,merged_model,models=None):
        '''
        @param tuple: tuple form wich the element is constructed
        @param models: use a specific model to construct this element, the models in this have a priority for their choice
        @type models: [L{OriginModel}]
        '''
        self._granny        = BioRelations()
        self.tuple          = tuple
        self.merge_elements = tuple.getElements()
        #--------------------------------------
        #copy an element of tuple as basis
        merge_element2clone=None
        types=[e.type for e in self.merge_elements]
        mixed_tuple= ('parameter' in types) and (('species' in types) or ('compartment' in types))#if parameter and a species/compartment are merged this is a mixed tuple
        if models:#a order of preferred models was specified
            if mixed_tuple:#if mixed tuple drop all simple elements
                models = [m for m in models if tuple.getElement(m).type != 'parameter']
            for m in models:
                if tuple.getElement(m):
                    merge_element2clone=tuple.getElement(m)#if a model is specified clone the element of the specified model
                    break
            else:
                #print 'models set but no el found with one of the models, switching to random model',models
                merge_element2clone=self.merge_elements[0]

        else:
            if mixed_tuple:
                for e in self.merge_elements:#prefere the more complex element
                    if e.type=='species' or e.type=='compartment':
                        merge_element2clone = e
                        break
            else:
                merge_element2clone=self.merge_elements[0]
        #--------------------------------------
        #initialize the 
        cloned_libsbml_element=merge_element2clone.libsbml_element.clone()
        doc = cloned_libsbml_element.getSBMLDocument()
        #self._granny.setNewId(cloned_libsbml_element,merged_model) 
        OriginElement.__init__(self,cloned_libsbml_element,merged_model)

        #set unique id TODO need to register that globally
        self.mid = self.id+'_merge'

        self.annoation_conflict_suggestions = None
        for type in 'bioEntity bioQuantity bioQLocation modelStatement'.split():
            for prop in 'inConflict conflictResolved'.split():
                setattr(self,'%s_%s'%(type,prop),False)

    def dependingElements(self):
        '''
        overwrite depending Elements since this only concerncs L{OriginElement}s 
        '''
        #print 'md e'
        pass

    def setId(self,id):
        '''set the local id and the id of the libsbml_element'''
        self.libsbml_element.setId(id)
        self.id=id

    def getId(self):
        '''return the local id'''
        return self.id
    def resolve(self):
        '''
        resolve inter class problems
        check for general problems and resolve them
        '''
        error_str=''
        #check for amount concentration problem
        if self.getQuantityType()=='Amount' and self.getQuantity()[0]=='Concentration':
            error_str='The element type is set to amount but an initial concentration value is set, this is not possible.'

        if error_str: raise MergeError(error_str)


    def resolve_bioQuantity(self,quantity_type,unit,location,participants):
        '''
        try to resove all bioQunatity problems  
        '''
        #resolve quantity type
        if quantity_type: self.setQuantityType(quantity_type)
        #resolve unit
        if unit: self.setUnit(unit)
        #resolve location
        if location: self.setLocation(location)
        #resolve bioEntities (participants of reactions conflict)
        if participants: 
            self.setBioEntities(participants)

        self.bioQLocation_conflictResolved=True
        self.bioQuantity_conflictResolved=True
    def resolve_modelStatement(self,simple_attributes,statement,initialAssignment):
        '''
        resolve conflicts in all Merge_ModelStatement values by setting it to the values given as input
        '''
        if simple_attributes:
            for atr,val in simple_attributes.items():
                getattr(self,'set'+atr.capitalize())(val)
        if statement: 
            self.setStatement(statement)
        if initialAssignment: self.setInitialAssignment(initialAssignment)
        self.modelStatement_conflictResolved=True

    def resolve_bioEntity(self,sbo_term=None):
        dbs_is=set()
        for an in self.getMiriamAnnotations():
            if an.qualifier=='is':
                if set([an.db]).issubset(dbs_is):#replace with database comparison
                    self.bioEntity_conflictResolved=False
                    raise MergeError('You are trying to merge two MIRIAM annotations with qualifier \'is\' pointing to the database %s.'%an.db)
                else:
                    dbs_is.add(an.db)

        if sbo_term:
            self.setSboAnnotation(sbo_term.getSboTerm())
        self.bioEntity_conflictResolved=True

    def prepare(self,find_conflics=True):
        #prepare self
        #----------------------------------------------------
        #modify math and other pointers
        if self.getStatement(): 
            #print 'translating math'
            self.setStatement(self._granny.updateMathElementIDs(self.getStatement(),self.merge_model))

        #TODO why was this block needed, why does it make problems now, DELETE?
        if self.type == 'reaction' and False:
            new={'Reactant':[],'Modifier':[],'Product':[]}
            for t in list(new.keys()): 
                for e  in self.getBioEntities()[t]:
                    new[t].append(e)
            #print 'x1 ',new
            #print 'x2 ',self.getBioEntities()
            self.setBioEntities(new)
        #find conflicts
        self.attribute_value_choice={}
        self.statement_choice=[]
        self.initialAssignment_choice=[]
        self.location_choice=[]
        self.bioEntities_choice={'Reactant':[],'Modifier':[],'Product':[]}
        self.bioQuantityType_choice=[]
        self.unit_choice=[]
        self.annoation_choice={}
        self.sbo_choice=[]
        #----------------------------------------------------
        if find_conflics:
            #print 'preparing merge'#DEBUG
            self.prepareModelStatement()
            self.prepareBioQuantity()#order is important quantiy then location

            self.prepareLocation()
            #reset unwanted things
            self.unsetAnnotations()
            self.prepareAnnotations()

    def prepareModelStatement(self):
        simple_attributes=['quantity','charge','reversible','fast','constant','size']
        initial_values={}
        choice_values={}
        statements=[]
        initialAssignments=[]
        for e in self.merge_elements:
            for atr in simple_attributes:
                try:
                    if not atr in initial_values:
                        initial_values[atr]=getattr(e,'get'+atr.title())()
                        choice_values[atr]=[getattr(e,'get'+atr.title())()]
                    elif initial_values[atr]!=getattr(e,'get'+atr.title())():
                        choice_values[atr].append(getattr(e,'get'+atr.title())())
                        self.modelStatement_inConflict=True
                        #if atr=='quantity':print 'conflict for atr %s in elem %s '%(atr,self.repr)#DEBUG
                except NotImplementedError:
                    pass
            if e.getStatement():

                math=self._granny.updateMathElementIDs(e.getStatement(),e.merge_model)
                if e.type!='reaction':#let rule point to merged element
                    math.setId(self.getId())
                    #TODO whats wrong now?
                    pass
                if not statements:
                    statements.append(math)
                elif math.getFormula() not in [x.getFormula() for x in statements]:#TODO comparision with math comparison code from copasi
                    statements.append(math)
                    self.modelStatement_inConflict=True
            if e.getInitialAssignment():
                math=self._granny.updateMathElementIDs(e.getInitialAssignment(),e.merge_model)
                if e.type!='reaction':#let initialAssigment point to merged element
                    math.setId(self.getId())
                if not initialAssignments:
                    initialAssignments.append(math)
                elif math.getFormula() not in [x.getFormula() for x in initialAssignments]:
                    initialAssignments.append(math)
                    self.modelStatement_inConflict=True

        for atr in simple_attributes:
            if atr in choice_values and len(choice_values[atr])<2:
                choice_values.pop(atr)
        self.attribute_value_choice=choice_values
        self.statement_choice=statements
        self.initialAssignment_choice=initialAssignments

    def prepareLocation(self):
        location_choice=None
        if self.type=='reaction':#reactands and products may be in differnent comapartments
            for n in range(1):
                locations=[]#TODO TODO TODO not good since multiple locations vs one location in e.g. Product will not be recognized
                new=True
                for e in self.merge_elements:
                    for x in e.getLocation()[n]:
                        if x.getTuple():x=x.getTuple().getMergedElement()
                        if new: locations.append(x)
                        if x.mid not in [x.mid for x in locations]:
                            self.bioQLocation_inConflict=True
                            break
                    new=False
        elif self.type in ['species','compartment']:
            location_choice=[]
            for e in self.merge_elements:
                if e.getLocation():
                    #print 'merging location'#DEBUG
                    if e.getLocation().getTuple():
                        if not location_choice: location_choice.append(e.getLocation().getTuple().getMergedElement())
                        elif e.getLocation().getTuple().getMergedElement() not in location_choice:
                            location_choice.append(e.getLocation().getTuple().getMergedElement())
                            self.bioQLocation_inConflict=True
                            #print 'location conflict'DEBUG
                    elif e.getLocation() not in location_choice: 
                        if not location_choice: location_choice.append(e.getLocation())
                        else:location_choice.append(e.getLocation())
                        self.bioQLocation_inConflict=True
                        #print 'location conflict1'#DEBUG
            location_choice=[x for x in location_choice if x]
        #print 'wer ',self.id,self.bioQLocation_inConflict#DEBUG
        self.location_choice=location_choice

    def prepareBioQuantity(self):
        unit_conflict=False
        participant_conflict=False
        type_conflict=False

        bqtype=self.getQuantityType()
        chooseTypes=[]
        unit_choice=[]
        bioEntities_choice={'Reactant':[],'Modifier':[],'Product':[]}
        unit=''
        #print 'merging bq',tuple.merged_element.repr
        for e in self.merge_elements:
            #Unit
            uid=''
            curunit=e.getUnit() or []
            uid=str(curunit)
            if not unit:
                unit=uid
                unit_choice.append(curunit)
            elif unit!=uid:
                unit_choice.append(curunit)
                unit_conflict = True
            #bioQuantityType e.g. amount,concentration or volume,area,point
            if not e.getQuantityType() in chooseTypes:
                chooseTypes.append(e.getQuantityType())
                type_conflict = True
            #if e.getQuantityType()!=bqtype:
            #   print 'bq in conflict type',self.bioQuantityType,bqtype
            #   chooseTypes=[e.getQuantityType(),bqtype]#TODO the max number of bqtypes is 2 that is why we can do that but it may have to change
            #   self.bioQuantity_inConflict=True
            #bioEntities
            if self.type == 'reaction':
                for type,belist in e.getBioEntities().items():
                    elist=[]
                    for et in belist: 
                        if et.getTuple():
                            elist.append(et.getTuple().getMergedElement())
                        else:
                            elist.append(et)
                    if not bioEntities_choice[type]: bioEntities_choice[type]=[elist]
                    else:
                        for prev_elist in bioEntities_choice[type]:
                            if set(prev_elist)!=set(elist):
                                bioEntities_choice[type].append(elist)
                                participant_conflict = True
                                #print 'participants in conflict',self.id#DEBUG
        self.unit_choice=unit_choice
        if len(chooseTypes)<2: 
            chooseTypes=[]
            type_conflict = False
        self.bioQuantity_inConflict=(unit_conflict or type_conflict or participant_conflict)
        self.bioQuantityType_choice=chooseTypes
        for type in ['Reactant','Modifier','Product']:
            if len(bioEntities_choice[type])<2: 
                bioEntities_choice[type]=None
        self.bioEntities_choice=bioEntities_choice

    def prepareAnnotations(self):
        '''
        copy all annotations of all merge tuple merge elements into this elements
        search for conflicts
            - same id different qualifier
            - same db differnet id, qualifier is in both cases
        @return: return a list of annations that have conflicting qualifiers
        @rtype: [Annotation]
        '''
        allAnnotations={}
        retAn={}
        db2id={}
        for e in self.merge_elements:
            for an in e.getMiriamAnnotations():
                #check for double is annotation to same database
                if an.qualifier=='is':#TODO not very suffisticated yet since qualifier type is ignored
                    if an.db in db2id:
                        if db2id[an.db]!=an.id:
                            #print 'found same db but different id with qualifier is in',self.repr#DEBUG 
                            self.bioEntity_inConflict=True
                    else:
                        db2id[an.db]=an.id
                try:
                    (qt,q)=allAnnotations[str(an)][0]#TODO poor solution, fials if there are more than 2 different qualifiers
                    #print 'testing:',allAnnotations[str(an)],str(an),an.qualifier_type,an.qualifier
                    if qt!=an.qualifier_type or q!=an.qualifier:
                        allAnnotations[str(an)].append((an.qualifier_type,an.qualifier))
                        self.bioEntity_inConflict=True
                        retAn[str(an)]=allAnnotations[str(an)]
                except KeyError as e :
                    allAnnotations[str(an)]=[(an.qualifier_type,an.qualifier)]
                    self.addAnnotation(an)
        self.annoation_choice = retAn
        sbo_terms = list(set([e.sbo_annotation_element for e in self.merge_elements if not e.getSboTerm() == None]))
        #print 'sbo_terms ',sbo_terms#DEBUG
        if sbo_terms:
            self.setSboAnnotation(sbo_terms[0].getSboTerm())
        if len(sbo_terms) > 1:
            self.sbo_choice = sbo_terms
            self.bioEntity_inConflict = True
        #print 'sbo_choice ',[e.getSboTerm() for e in self.sbo_choice]

    def checkAnnotations(self):
        '''
        function is needed when annotations are deleted to see if there are still conflicts, 
        this is not the same as conflict resolution
        maybe merge this function with L{prepareAnnotations} since it is basically a stripped copy of it
        this function just changes the internal state so no input and output is needed
        '''
        allAnnotations={}
        retAn={}
        db2id={}
        for e in self.merge_elements:
            for an in e.getMiriamAnnotations():
                #check for double is annotation to same database
                if an.qualifier=='is':#TODO not very suffisticated yet since qualifier type is ignored
                    if an.db in db2id:
                        if db2id[an.db]!=an.id:
                            #print 'found same db but different id with qualifier is in',self.repr#DEBUG
                            self.bioEntity_inConflict=True
                    else:
                        db2id[an.db]=an.id
                try:
                    (qt,q)=allAnnotations[str(an)][0]#TODO poor solution, fials if there are more than 2 different qualifiers
                    if qt!=an.qualifier_type or q!=an.qualifier:
                        allAnnotations[str(an)].append((an.qualifier_type,an.qualifier))
                        self.bioEntity_inConflict=True
                except KeyError:
                    allAnnotations[str(an)]=[(an.qualifier_type,an.qualifier)]
class MergeTuple:
        '''
        collection of matching biological entites (L{OriginElement}s)
        link to a merged entity (L{MergedElement})
        functions for adding and removing biological entites
        '''
        def __init__(self,el1,el2):
                '''
                The tuple is initialized with two elements
                '''
                self._granny=BioRelations()
                self.elements={}
                self.pair_score={}
                self.pair_location_score={}
                self.type=el1.type
                self._merged_element=None
                self.sig_MergedElement_created=None
                self.addElement(el1,el2)


        def connect_view_functions(self,external_fkt):
                '''
                connect external view function so it can update if the Merged Element was created
                '''
                #print 't:connected view fkt'#DEBUG
                self.sig_MergedElement_created=external_fkt

        def setMergedElement(self,merged_model,prepare=True,models=None):
                '''
                create L{MergedElement}, this function can be used for the initial createn and the updating
                it will send a signal that the element was created/refreshed
                @param prepare: prepare the merged element, fill choice objects
                @type prepare: bool
                @rtype: L{MergedElement}
                @type model: [L{OriginModel}]
                '''
                #create a new merged element if it did not exitst of if the tuple content has changed TODO isnt this always the case?
                if not self._merged_element or prepare: self._merged_element=MergedElement(self,merged_model,models)
                merged_model.addElement(self._merged_element)
                #tuple content has changed prepare merged element (e.g. update ids in math) and detect current conflicts
                find_conflics=True
                if models: find_conflics=False #do not find conflicts if prefrerred order of models was specified
                if prepare: self._merged_element.prepare(find_conflics)
                #send a signal to external views
                if self.sig_MergedElement_created: 
                        #print 'tuple: sending signal merged e crea'#DEBUG
                        self.sig_MergedElement_created()
                #get a unique id for this element
                oldId=self._merged_element.getId()
                self._granny.setNewId(self._merged_element,merged_model)
                newId=self._merged_element.getId()

                #return the merged element
                return (oldId,newId)

        def getElements(self):
                '''
                return all L{MergedElement}s of this tuple
                @rtype: [L{OriginElement}]
                '''
                return list(self.elements.values())

        def getElement(self,model):
                '''
                get an element that belongs to a specific model
                @type model:L{OriginModel}
                @rtype: l{OriginElement} or None
                '''
                try:
                        return self.elements[model.mid]
                except KeyError:
                        return None

        def getMergedElement(self):
                '''
                get the merged element of this tuple, this can fail if L{setMergedElement} was not called before
                @rtype: L{MergedElement}
                '''
                return self._merged_element

        def addElement(self,old_element,new_element,force=False):
                '''
                @param force: adding of element is requred (eg user interaction)
                @type force: bool
                @param old_element: supposedly old element (this name is not needed anymore)
                @type old_element:L{OriginElement}
                @type new_element:L{OriginElement}
                '''
                def addNewElement(old_e,new_e):

                        self.elements[new_e.merge_model.mid]=new_e
                        self.elements[old_e.merge_model.mid]=old_e

                        #store scores
                        for element in self.getElements():
                                #if not (element.mid and new_e.mid):sdffj
                                if frozenset([element,new_e]) not in self.pair_score and not element==new_e:
                                        self.pair_score[frozenset([element.mid,new_e.mid])]=self._granny.compareElements(element,new_e)
                                        self.pair_location_score[frozenset([element.mid,new_e.mid])]=self._granny.compareLocation(element,new_e)
                        #update similarity list and scores of each element
                        allcandidates=[]
                        for element in self.getElements():
                                allcandidates+=element.getTupleCanidates()
                        for element in self.getElements():
                                element.addTupleCandidates(allcandidates)

                        #set the new tuple
                        new_e.setTuple(self)
                        old_e.setTuple(self)

                if new_element.merge_model.mid not in self.elements:
                        #print 'first elem of model in this tuple, will add'
                        if self.elements:
                                #print 'adding ',self.elements[old_element.merge_model.mid].mid,old_element.mid,' xxx ',new_element.mid
                                if self.elements[old_element.merge_model.mid].mid!=old_element.mid:
                                        print('yeiks',self.elements[old_element.merge_model.mid].mid,old_element.mid)
                        addNewElement(old_element,new_element)
                elif force:
                        #print 'add forced, will add'
                        if new_element.merge_model.mid in self.elements:
                                try:
                                        self.remElement(self.elements[new_element.merge_model.mid],True)
                                except MergeTupleError:
                                        pass#we do not want to stop here
                        addNewElement(old_element,new_element)
                else:
                        location_score=self._granny.compareLocation(new_element,old_element)
                        score=self._granny.compareElements(new_element,old_element)
                        if self.pair_location_score[frozenset([self.elements[new_element.merge_model.mid].mid,old_element.mid])]==location_score and self.pair_score[frozenset([self.elements[new_element.merge_model.mid].mid,old_element.mid])]<score:#if same location score and #if higher pair score
                                #print 'same location score, but higher similarity found'
                                #remove old element from same model
                                self.remElement(self.elements[new_element.merge_model.mid])
                                #add new element
                                addNewElement(old_element,new_element)

                        elif self.pair_location_score[frozenset([self.elements[new_element.merge_model.mid].mid,old_element.mid])]<location_score:
                                #print 'higher location score found'
                                #remove old element from same model
                                self.remElement(self.elements[new_element.merge_model.mid])
                                #add new element
                                addNewElement(old_element,new_element)
                        else:
                                #print 'not added, but updating simlist'#DEBUG
                                allcandidates=[new_element]+new_element.getTupleCanidates()
                                for element in self.getElements():
                                        allcandidates+=[element]+element.getTupleCanidates()
                                for element in self.getElements():
                                        #print allcandidates
                                        element.addTupleCandidates(allcandidates)
                                new_element.addTupleCandidates(allcandidates)

                                raise MergeTupleError('Element %s Candidate for multiple Tuples, rejected for This Tuple'%new_element.repr)

        def remElement(self,element,force=False):
                '''
                remove element form tuple  and add 
                remaining elements as candidates
                @param force: force removal but do not destroy tuple
                @param direct: user requested removal, update candidate list
                '''
                remove_e=self.elements.pop(element.merge_model.mid)

                #update candiates list
                for remainig_element in self.getElements():
                        remainig_element.addTupleCandidate(element)
                        element.addTupleCandidate(remainig_element,True)
                        self.pair_score.pop(frozenset([element.mid,remainig_element.mid]))
                        self.pair_location_score.pop(frozenset([element.mid,remainig_element.mid]))

                element.setTuple(None)  

                if len(self.elements)==1 and not force: #if only two elements left and the removing was not forced, forced means that a new elemnt will be added
                        #print 'els left: ',len(self.elements),'force: ',force

                        for laste in list(self.elements.values()):
                                self.remElement(laste)

                        #TODO destroy this instance - but how?
                        #print 'destroy this tuple'
                        raise MergeTupleError('Less than 2 elements in Tuple')


class Merger:
    '''
    class for merging multiple libSBML models
    '''
    def __init__(self,documents,callback_statusbar_fkt,popup_msg_fkt=lambda x: 0):
        '''
        @type models: [L{Model}]
        @type callback_statusbar_fkt: function pointer
        @param callback_statusbar_fkt: function of a view to display status messages send during merging
        '''
        import codecs,sys
        _reader = libsbml.SBMLReader()
        global setStatusMessage     #get gobal var to manipulate it, calling does not require this
        setStatusMessage = callback_statusbar_fkt
        #init BioRelations since it is a single class and if the init with the L{Merger} is not done here L{Merger} will miss later on
        self._granny = BioRelations(self)
        self.libsbml_documents=[]
        warnings=''
        for doc in documents:
            doc.setLevelAndVersion(2,4,False)
            try:
                codecs.open("doc_file.xml",'w', encoding='utf-8').write(str('<?xml version="1.0" encoding="UTF-8"?>\n%s'%doc.toSBML()))
            except:
                raise MergeError('Could not write temporary SBML file. %s %s'%(sys.exc_info()[0],sys.exc_info()[1]))
            ld = _reader.readSBML("doc_file.xml")
            os.remove("doc_file.xml")
            model = ld.getModel()
            lm = ld.getModel()
            #////////////////////////////////////////////////////////////
            #check for all elements we cannot handle yet
            #AlgebraicRule, Event, Constraint, SpeciesType, CompartmentType
            repr       = lm.getName() or lm.getId()
            unsupported = ['CompartmentTypes','SpeciesTypes','Events','Constraints']
            for uns in unsupported:
                if getattr(lm,'getNum'+uns)():
                    warnings += 'Model %s contained %s, this is currently not supported and will be ignored'%(repr,uns)
            for r in list(lm.getListOfRules()):
                if isinstance(r,libsbml.AlgebraicRule):
                    warnings += 'Model %s contained AlgebraicRule(s), this is currently not supported and will be ignored'%repr
            #////////////////////////////////////////////////////////////
            self.libsbml_documents.append(ld)
        if len(self.libsbml_documents)<2:
            raise MergeError('Less Than 2 Models to Merge')
        if warnings:
            popup_msg_fkt('Warnings encounterd on loading',warnings)
        #--------------------------------------------------------------- 
        #---------------------------------------------------------------
        #set unique ids and meta ids in every model
        elements_w_id=['Species', 'Compartments', 'Reactions', 'Parameters', 'FunctionDefinitions', 'InitialAssignments', 'Rules']
        oldId_newId={}
        self.merge_models    = []
        for mid,document in enumerate(self.libsbml_documents):
            model = document.getModel()
            oldId_newId[mid]=[]
            for etype in elements_w_id:
                l = getattr(model,'getListOf'+etype)()
                for e in list(getattr(model,'getListOf'+etype)()):
                    #======================
                    if not etype in ["FunctionDefinitions", 'InitialAssignments']:             # FunctionDefinitions do not support Annotations but still need new ids and meta ids
                        ae = semanticSBML.annotate.ElementAnnotations(e,True)
                        annotstore = ae.getAnnotations()
                        ae.unsetAnnotations()#will not work since the cvterm is not delted
                        e.unsetAnnotation()
                    newMetaId = '%s_%s'%(e.getMetaId(),mid)
                    e.setMetaId(newMetaId)
                    if not etype in ["FunctionDefinitions", 'InitialAssignments']:             # FunctionDefinitions do not support Annotations but still need new ids and meta ids
                        for an in annotstore:
                            ae.addAnnotation(an)
                    #======================
                    newId = '%s_%s'%(e.getId(),mid)
                    oldId_newId[mid].append((e.getId(),newId))
                    e.setId(newId)
        #set all pointers to new ids
        elements_w_id=['Species', 'Compartments', 'Reactions', 'Parameters', 'FunctionDefinitions', 'InitialAssignments', 'Rules']
        for mid,document in enumerate(self.libsbml_documents):
            model = document.getModel()
            for etype in elements_w_id:
                for e in list(getattr(model,'getListOf'+etype)()):
                    for oldId,newId in oldId_newId[mid]:
                        self._granny.updateId(e,oldId,newId)

            #open("testpremerge.xml",'w').write(model.toSBML())#DEBUG
            #raw_input('pause')#DEBUG       
        #---------------------------------------------------------------
        self.merge_models=[]
        error_str=''
        for mid,document in enumerate(self.libsbml_documents):
            model = document.getModel()
            try:
                self.merge_models.append(OriginModel(model,mid))
            except semanticSBML.annotate.AnnotationCreationError as e:
                error_str += 'The model "%s" can not be used for merging because: %s'%(model.getName() or model.getId(), e)
                continue
        if error_str: raise MergeError(error_str)

        self._getMergeTuples()
        self.libsbm_elem2merge_elem={}

        for m in self.merge_models:
            for e in m.getElements():
                self.libsbm_elem2merge_elem[e.libsbml_element]=e

    def copyAttributes(self,e1,e2):
        '''
        '''
        try:
            e2.setCompartment(e1.getCompartment())
        except:
            pass

        getterList = ['get'+x for x in ['Id','Name']]
        for getter in getterList:
            setter = 's'+getter[1:]
            if setter in dir(e1):
                getattr(e2,setter)(getattr(e1,getter)())


    def makeCopyDoc(self,doc):
        '''
        generate brute force copy of the document (since .clone() and .deepcopy() do not work and .copy() is too shallow)
        '''
        model  = doc.getModel()
        self.__newDoc = libsbml.SBMLDocument()
        self.newModel = self.__newDoc.createModel()
        #self.newModel.setName(model.getName())
        #self.newModel.setId(model.getId())
        blacklist = ['getListOf'+x for x in ['Rules']]

        for getList in dir(model):
            if getList.startswith('getList') and getList not in blacklist:
                elementType = getList[9:]
                if elementType != 'Species':
                    elementType = elementType[:-1]
                for element in getattr(model,getList)():
                    newElement = getattr(self.newModel,'create'+elementType)()
                    self.copyAttributes(element,newElement)

        self.newModel.setId(model.getId())
        self.newModel.setName(model.getName())

        return self.__newDoc


    def _getMergeTuples(self):
        model_pop_list=copy.copy(self.merge_models)
        xnt=len(self.merge_models)                    #model_pop_list)
        nummod=0
        while model_pop_list:
            m=model_pop_list.pop(0)#get a model
            nummod+=1
            numel=len(model_pop_list)
            for xi,om in enumerate(model_pop_list):#get a 2nd model
                if not xi+1 == nummod+1:
                    setStatusMessage('Merge Core: Comparing Model %s with Model %s'%(xi+1,nummod+1))
                    for i,el in enumerate(m.getElements()):#go through elements
                        for oel in om.getElementsByType(el.type):#go through elements of 2nd model
                            score = self._granny.compareElements(el,oel)
                            if score>0: #add pair to candidate list
                                try:
                                    self._addElement2Tuple(el,oel)
                                except MergeTupleError as e:
                                    sys.stderr.wirte('Error occured while adding to tuple:',e)
        #TODO we need more than a 2nd run, as much runs as long as elements get kicked out of tuples? rethink if this is true, I don not want to get stuck in an endless loop
        #alternative is a global matrix containing all match scores where bes matches are selected
        model_pop_list=copy.copy(self.merge_models)
        nummod=0
        while model_pop_list:
            m=model_pop_list.pop(0)#get a model
            nummod+=1
            numel=len(model_pop_list)
            for xi,om in enumerate(model_pop_list):#get a 2nd model
                setStatusMessage('Merge Core: 2nd Run Comparing Model %s with Model %s'%(xi+1,nummod+1))
                for i,el in enumerate(m.getElements()):#go through elements
                    for oel in om.getElementsByType(el.type):#go through elements of 2nd model
                        if not el.getTuple() and not oel.getTuple() and (el.getTupleCanidates() or oel.getTupleCanidates()):
                            score = self._granny.compareElements(el,oel)
                            if score>0:
                                try:
                                    self._addElement2Tuple(el,oel)
                                except MergeTupleError as e:
                                    sys.stderr.write('Error occured while adding to tuple:',e)

        tlen=len(self.getTuples())

        #--------------------------------------------
        #create new model
        if not '__hidden_new_sbml_document' in dir(self):
            self.__hidden_new_sbml_document = []
            self.__hidden_new_sbml_document.append(libsbml.SBMLDocument())
            new_libsbml_model=self.__hidden_new_sbml_document[-1].createModel()
            self.merged_model =  MergedModel(new_libsbml_model, self.merge_models)   #self.libsbml_documents)  


        oldid2newid=[]
        for tuple in self.getTuples():#create base elements
            oldid2newid.append(tuple.setMergedElement(self.merged_model,False))
        for tuple in self.getTuples():
            self.merged_model.addElement(tuple.getMergedElement())
        #--------------------------------------------
        #prepare merged elemtns - find conflicts ...
        for i,tuple in enumerate(self.getTuples()):
            setStatusMessage('Merge Core: Merging Tuple %s/%s'%(i+1,tlen))
            tuple.getMergedElement().prepare()  

    def resolveUsingModel(self,merge_models):
        '''
        ignore all conflicts and construct merged elements using only one element in a tuple
        first model has the highest priority
        @type merge_models: [L{OriginModel}]
        '''
        for t in self.getTuples():
            t.setMergedElement(self.merged_model,True,merge_models)

    def createSBMLmodel(self,merged_model_name=''):
        '''
        create a new SBML model from the semanticSBML abstraction
            - go tough tuples and and non tuple elements
                - collect units (in a set to prevent doubling)
                - get rules and add them to the new model
                - get initial assignments and add them to the new model
            - add units
            - add mergable elements
        @rtype: libsbml.Document
        '''
        #create a blank model
        document = libsbml.SBMLDocument()
        document.setLevelAndVersion(2,4,False)
        model    = document.createModel()

        #add all mergable elements
        all_merge_elements = []
        math_used          = [] 

        #--------------------------------
        #collect elements and check for unresolved conflicts in tuples
        error_str=''
        for tuple in self.getTuples():
            me = tuple.getMergedElement()
            if (me.bioEntity_inConflict and not me.bioEntity_conflictResolved) or (me.bioQuantity_inConflict and not me.bioQuantity_conflictResolved) or (me.modelStatement_inConflict and not me.modelStatement_conflictResolved):
                error_str+='Element %s has conflicting values. '%me.repr
            else:
                all_merge_elements.append(me)

        if error_str: raise MergeError(error_str[:-1])
        for e in self.getNonTupleElements():
            all_merge_elements.append(e)



        #--------------------------------
        #merge units by collecting all available units in all models and merging them non redundant
        #TODO this still contains the bug of not detecting different unit definitions with same id!!!!
        def getLibsbmlUnitDefinition(unit_def,no_list=False):
            if not unit_def:    return []
            if unit_def.type=='libsbml':#recursion anchor
                if no_list:     return unit_def
                else:       return [unit_def]
            elif unit_def.type=='composed':
                return [getLibsbmlUnitDefinition(ud,True) for ud in unit_def.unit_defs]
            else: 
                if no_list:     return None
                else:       return []
        units_used = []
        for e in self.getNonTupleElements()+[tuple.getMergedElement() for tuple in self.getTuples()]:
            units_used+=getLibsbmlUnitDefinition(e.getUnit())
        #--------------------------------
        # check in all reactions for internal parameters referring to units
        unit_ids = set()
        all_elements = [t.getMergedElement() for t in self.getTuples()]+self.getNonTupleElements()
        for reaction_me in [e for e in all_elements if e.type == 'reaction']:
            try:
                for parameter in list(reaction_me.libsbml_element.getKineticLaw().getListOfParameters()):
                    unit_ids.add(parameter.getUnits())
            except:
                pass
        #now go into each model and check if there is this unit we need TODO this is a hack
        for xmodel in self.merge_models:
            for uid in unit_ids:
                try:
                    units_used += getLibsbmlUnitDefinition(xmodel.global_unit_defs[uid])
                except KeyError:
                    pass
        #clear None entries
        units_used = [ud for ud in units_used if not ud==None]    

        #removing duplicates
        units_used = list(dict(zip([str(ud) for ud in units_used],[ud.libsbml_ud for ud in units_used])).values())
        for u in units_used:
            model.addUnitDefinition(u)
        #--------------------------------
        #add attached elements (rules,initialAssignments) to model
        #in tuple all formulas have been translated
        #def XgetFunctionDefinition(object):#remove when libsbml fixed
        #   try:
        #       for fd in object.functionDefinitions:
        #           model.addFunctionDefinition(fd)
        #   except AttributeError:
        #       pass
        #--------------------------------
        fd_math_str2id={}
        def getFunctionDefinition(merge_element):
            if merge_element.getInitialAssignment():
                object= merge_element.getInitialAssignment()
            elif merge_element.getStatement():
                object= merge_element.getStatement()
            else:
                return
            libsbml_element=merge_element.libsbml_element
            lofd=[]
            try:
                lofd=object.functionDefinitions
            except AttributeError:
                pass
            else:
                for fd in lofd:
                    if libsbml.formulaToString(fd.getMath()) in list(fd_math_str2id.keys()):
                        #translate linking math to already existing
                        self._granny.updateId(libsbml_element,fd.getId(),fd_math_str2id[libsbml.formulaToString(fd.getMath())])
                    else:
                        model.addFunctionDefinition(fd)
                        fd_math_str2id[libsbml.formulaToString(fd.getMath())]=fd.getId()

        for e in [t.getMergedElement() for t in self.getTuples()]+self.getNonTupleElements():
            getFunctionDefinition(e)

        for t in self.getTuples():
            e=t.getMergedElement()
            if e.getInitialAssignment():
                model.addInitialAssignment(e.getInitialAssignment())
            if e.getStatement(): 
                if not e.type=='reaction':
                    model.addRule(e.getStatement())#insert rule into the new model

        #in non tuple element formulas need to be translated
        for e in self.getNonTupleElements():
            if e.getInitialAssignment():
                self._granny.updateMathElementIDs(e.getInitialAssignment(),e.merge_model,False)
                e.getInitialAssignment().setId(e.getId())
                model.addInitialAssignment(e.getInitialAssignment())
            if e.getStatement(): 
                self._granny.updateMathElementIDs(e.getStatement(),e.merge_model,False)
                if not e.type == 'reaction':
                    e.getStatement().setId(e.getId())
                    model.addRule(e.getStatement())#insert rule into the new model
                #else:
                #   e.setStatement(e.getStatement())

        #------------------------------------------------------
        #translate compartment outsides
        for e in all_merge_elements:
            if e.type == 'compartment':
                if e.getLocation():
                    if e.getLocation().getTuple():
                        e.libsbml_element.setOutside(e.getLocation().getTuple().getMergedElement().getId())
                    else:
                        e.libsbml_element.setOutside(e.getLocation().getId())
            elif e.type in ['species']:
                if e.getLocation().getTuple():
                    e.libsbml_element.setCompartment(e.getLocation().getTuple().getMergedElement().getId())
                else:
                    e.libsbml_element.setCompartment(e.getLocation().getId())

        #translate reaction participants
        for e in all_merge_elements:
            if e.type=='reaction':
                e.setBioEntities(e.getBioEntities())
        #------------------------------------------------------
        #add all mergable elements to new model
        for elmnt in all_merge_elements:
            e=elmnt.libsbml_element
            obj = getattr(model,'add'+e.getElementName().capitalize())(e)
            #obj.setId(e.getId())
        #------------------------------------------------------
        #set name
        model.setName(merged_model_name)
        #------------------------------------------------------
        #aggregate MIRIAM annotations
        merge_model_annotations = semanticSBML.annotate.ElementAnnotations(model,True)#get annotations of model element
        for ld in self.libsbml_documents:
            lm = ld.getModel()
            for an in semanticSBML.annotate.ElementAnnotations(lm,True).getAnnotations():
                #check for double is annotation to same database
                if an.qualifier=='is':  #TODO not very suffisticated yet since qualifier type is ignored
                    an.setQualifier('bio','hasPart')
                merge_model_annotations.addAnnotation(an)

        #------------------------------------------------------
        #TODO add normal annotations
        #------------------------------------------------------
        return document

    def getTuples(self):
        '''
        get all L{MergeTuple}s in a list sorted by their type
        @rtype: [L{MergeTuple}]
        '''
        #TODO sort tuples fist
        ret=[]
        elements=[]
        for model in self.merge_models:
            elements+=model.getElements()
        tuples=set([e.getTuple() for e in elements])
        tuples.remove(None)
        for type in ['compartment', 'species', 'reaction', 'parameter']:
            ret+=[x for x in tuples if x.type==type]
        return ret


    def getMergeElement(self,element):
        '''
        @param element: libsbml.Species or libsbml.Compartment  or ...
        @rtype: L{OriginElement}
        '''
        for t in [libsbml.Species,libsbml.Compartment,libsbml.Reaction,libsbml.Compartment,libsbml.AssignmentRule,libsbml.RateRule,libsbml.AlgebraicRule,libsbml.Event,libsbml.Model,libsbml.Parameter]:
            if isinstance(element, t):
                try:
                    return self.libsbm_elem2merge_elem[element]
                except KeyError:
                    return None

    def getNonTupleElements(self,merge_model=None):
        '''
        get all L{OriginElement}s that do not belong to a L{MergeTuple}
        @rtype: [L{OriginElement}]
        '''
        if merge_model:
            return [x for x in merge_model.getElements() if not x.getTuple()]
        else:
            elements=[]
        for m in self.merge_models:
            elements += m.getElements()
        return [e for e in elements if not e.getTuple()]

    def getMergeElementsOtherModels(self,model_id,type):
        other_merge_elements=[]
        for model in self.merge_models:
            if model.mid != model_id:
                other_merge_elements += model.getElementsByType(type)
        return other_merge_elements

    def remElementFromTuple(self,element):
        '''
        remove an element from a tuple
        @param element: element that should be removed from its tuple (element must have a tuple)
        @type element: L{OriginElement}
        '''
        self._granny.elements_update_queue=set()
        #-------------------------------------------

        tuple = element.getTuple()
        TupleDestroyed=None
        try:
            tuple.remElement(element)
        except MergeTupleError as e:
            TupleDestroyed=e
        else:#if tuple was not destroyed update its merged element  
            tuple.setMergedElement(self.merged_model)

        #-------------------------------------------
        #update depending elements
        #print 'rem updating: ',self._granny.elements_update_queue
        for t in [e.getTuple() for e in self._granny.elements_update_queue if e.getTuple()]:
            t.setMergedElement(self.merged_model)
            setStatusMessage('Updating Merged Element %s'%t.getMergedElement().repr)
        if TupleDestroyed:
            #self.merge_model.rem#TODO TODO TODO clean up merged model
            raise MergeTupleError(e)


    def addElementToTuple(self,this_element,to_element):
        '''
        force an add of this_element to the tuple of to_element
        @type this_element: L{OriginElement}
        @type to_element: L{OriginElement}
        '''
        self._granny.elements_update_queue=set()
        #-------------------------------------------
        try:
            if this_element.getTuple():
                this_element.getTuple().remElement(this_element)
        except MergeTupleError as e:#only one element in tuple left old tuple deleted
            pass
        self._addElement2Tuple(this_element,to_element,True)
        #print 'MERGE: updating md e ',this_element.mid#DEBUG
        this_element.getTuple().setMergedElement(self.merged_model)

        #-------------------------------------------
        #updating depending elements
        #print 'add updating: ',self._granny.elements_update_queue
        for t in [e.getTuple() for e in self._granny.elements_update_queue if e.getTuple()]:
            t.setMergedElement(self.merged_model)
            setStatusMessage('Updating Merged Element %s'%t.getMergedElement().repr)
        return this_element.getTuple()

    def _addElement2Tuple(self,el1,el2,force=False):
        '''add an element to a merge tuple or create a new tuple if it does not exist yet'''
        if el1.getTuple() and el2.getTuple():#if both elements are already in tuples
            #print 'both in tuple',el1.mid,el2.mid#DEBUG
            if not el2.getTuple() == el1.getTuple():
                (el1,el2)=BioRelations().chooseBestTuple(el1,el2)#reorder elements so the first element is the one where the adding of the second element would fit better than to the adding of the first element to the second
                #now try adding the elements to the tuple of the other, if it fails, no problem there was a better element, that we do not want to replace 
                try:
                    old_tuple=el2.getTuple()
                    el1.getTuple().addElement(el1,el2,force)
                    try:
                        old_tuple.remElement(el2)
                    except MergeTupleError:#only one element left in tuple
                        pass
                except MergeTupleError:#adding was not successful try to add other element the other tuple
                    try:
                        old_tuple=el1.getTuple()
                        el2.getTuple().addElement(el2,el1,force)#should never get here with force=True
                        try:
                            old_tuple.remElement(el1)
                        except MergeTupleError:#only one element left in tuple
                            pass
                    except MergeTupleError:
                        pass
            else:
                pass
        elif el1.getTuple():#if the first element is already in a tuple 
            #print '1st in toupel',el1.repr,el2.repr
            try:
                if not force:#during buildup,just add the new element
                    el1.getTuple().addElement(el1,el2,force)
                else:#adding to tuple was forced, this means that el1 must be added to el2  
                    #print 'forced adding, but first remove'#DEBUG
                    #remove this_element/el1 from its old tuple
                    try:
                        #print 'remove forced in el2tup'#DEBUG
                        el1.getTuple().remElement(el1)
                    except MergeTupleError:#only one element left in tuple
                        pass
                        #print 'only one left after rearrangement of tuples'
                    #add it to its new tuple
                    self._addElement2Tuple(el1,el2,force)

            except MergeTupleError:#adding was not successful, not so bad
                #print 'no success'
                pass    
        elif el2.getTuple():
            #print '2nd in toupel',el2.repr,el1.repr
            try:
                el2.getTuple().addElement(el2,el1,force)
            except MergeTupleError:
                pass    
        else: #if non of the elements is in a tuple create a new tuple
            #print 'new tuple added'
            MergeTuple(el1,el2)

    def getTuple(self,element,element2=None):
        if element2:
            #print 'el2 set'
            if element.getTuple() == element2.getTuple():
                    #print 'and found',self.element2tuple[element.mid].getElements()
                    return element.getTuple()
            return None
        return element.getTuple()

