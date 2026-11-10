#python libraries
import re
import copy
from . import xmlobject
import os
import sys
#external libraries
import libsbml
#semanticSBML libraries
from . import merge_datamodel
import libSBAnnotation.database_new
import libSBAnnotation.config

g_known_units = None
g_bio_qualifier_list=[]
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

class BioRelations(object):
    '''
    singleton class for the comparison of elements
    creation of globally unique identifiers
    update of mathematical expressions that contain identifiers
    '''
    __instance = None

    def __new__(cls, *args, **kargs): 
        if cls.__instance is None:
                cls.__instance = BioRelations.Local(*args, **kargs)
        return cls.__instance        

    @classmethod
    def reset(cls):
        '''forget the state of the last merge'''
        cls.__instance = None

    class Local:

        def __init__(self,merger=None):
            self.merger=merger
            self._ldb = libSBAnnotation.database_new.Meta_DB()
            self.element_score=self.getScoreMatrix(libSBAnnotation.config.Config().getpath('resources','element_score_matrix'))
            self.idlist=[]
            self.model2id2newid={}
            self.elements_update_queue=set()


        def setNewId(self,element,merge_model):
            '''
            set an unique id for an element,
            if id is already unique, leave it, if not set a new one
            @param element: merge(d) element or libsmbl element
            @type element: MergedElement | OriginElement | libsmbl.(Species|...)
            '''
            #TODO set unique metaid too!!!!
            #----------------------------------------------------
            try:
                if element.newidset:
                    pass
            except AttributeError:
                element.newidset=True
            #----------------------------------------------------
            i  = 'M'
            id = element.getId()
            #----------------------------------------------------
            while element.getId() in self.idlist:
                element.setId('%s_%s'%(id,i))
                if i=='M':
                    i = 0
                i+=1
            self.idlist.append(element.getId())
            try:
                self.model2id2newid[merge_model][id]=element.getId()
            except KeyError:
                self.model2id2newid[merge_model]={id:element.getId()}


        def getMathElementNames(self,mathobj):
            '''
            get a list of strings that are identifiers of sbml elements used in the inserted math obj
            this is a recursive function
            @type mathobj: libSBML.math
            @return: list of identifier strings
            @rtype: [str]
            '''
            if mathobj.getNumChildren()==0:
                if mathobj.isName():
                    return [mathobj.getName()]
                else:
                    return []
            else:
                names=[]
                #check all children, recursive
                for i in range(mathobj.getNumChildren()):
                    childnames=self.getMathElementNames(mathobj.getChild(i))
                    if childnames:names+=childnames
                #get this name/element
                if mathobj.isName() or mathobj.isFunction():
                    names.append(mathobj.getName())
                return names

        def updateMathElementIDs(self,mathobj,merge_model,create_copy=True):
            '''
            update all IDs in mathematical expressions with IDs of merged elements if the old element is merged
            '''
            if create_copy: mathobj = mathobj.clone()
            math=mathobj
            try:
                math = mathobj.getMath()
            except AttributeError:
                pass
            if math.getNumChildren()==0:
                #--------------------------------------------------------
                if math.isName() or math.isFunction():
                    try:
                        if merge_model.id2merge_element[math.getName()].getTuple():#check if tuple is available
                            math.setName(merge_model.id2merge_element[math.getName()].getTuple().getMergedElement().getId())
                        elif self.model2id2newid[merge_model][math.getName()]!=math.getName():#check if id was updated
                            math.setName(self.model2id2newid[merge_model][math.getName()])   # the original assigned to math.setName instead of calling it
                    except KeyError:#no!? internal parameter!
                        #print 'HUH internal param ',math.getName(),merge_model.repr
                        pass
                #elif math.isFunction():
                #   if self.model2id2newid[merge_model][math.getFunction()]!=math.getFunction():
                #       math.setFunction = self.model2id2newid[merge_model][math.getFunction()]

                #--------------------------------------------------------
            else:
                #check all children, recursive
                for i in range(math.getNumChildren()):
                    childnames=self.updateMathElementIDs(math.getChild(i),merge_model,False)
            #--------------------------------------------------------
            #check for fkt defs
            try:
                newfds=[]
                for ofd in mathobj.functionDefinitions:
                    xfd = ofd
                    if create_copy: xfd = ofd.clone()
                    newfds.append(self.updateMathElementIDs(xfd,merge_model,False))
                if create_copy: mathobj.functionDefinitions=newfds
            except AttributeError:# no fkt def attached
                pass
            #--------------------------------------------------------
            return mathobj


        def chooseBestTuple(self,el1,el2):
            '''
            reorder elements so that adding the second element to the tuple of the first gives a better result than adding the first element to the tuple of the second
            @type el1: L{OriginElement}
            @parameter el1: L{OriginElement} that belongs to a L{MergeTuple}
            @rtype: (OriginElement,OriginElement)
            @return: reorder list of the input
            '''
            t1=el1.getTuple()
            e1_combined_pair_score=0
            e1_combined_location_score=0
            for e in t1.getElements():
                if not e.mid==el1.mid:
                    e1_combined_pair_score+=t1.pair_score[frozenset([el1.mid,e.mid])]
                    e1_combined_location_score+=t1.pair_location_score[frozenset([el1.mid,e.mid])]
            t2=el2.getTuple()
            e2_combined_pair_score=0
            e2_combined_location_score=0
            for e in t2.getElements():
                if not e.mid==el2.mid:
                    e2_combined_pair_score+=t2.pair_score[frozenset([el2.mid,e.mid])]
                    e2_combined_location_score+=t2.pair_location_score[frozenset([el2.mid,e.mid])]

            if e1_combined_location_score==e2_combined_location_score:
                if e1_combined_pair_score>=e2_combined_pair_score:
                    return (el1,el2)
                else:
                    return (el2,el1)
            elif e1_combined_location_score>=e2_combined_location_score:
                return (el1,el2)
            else:
                return (el2,el1)
        def compareElements(self,el1,el2):
            '''
            compare the annotations of two elements and return a similarity score
            @param el1: element containing annotations
            @param el2: element containing annotations
            @type el1: libSBML base element or L{OriginElement}
            @type el2: libSBML base element or L{OriginElement}
            '''
            if isinstance(el1, merge_datamodel.OriginElement):
                el1annotations=el1.getMiriamAnnotations()
            else:
                el1annotations=ElementAnnotations(el1).getAnnotations()
            if isinstance(el2, merge_datamodel.OriginElement):
                el2annotations=el2.getMiriamAnnotations()
            else:
                el2annotations=ElementAnnotations(el2).getAnnotations()
            score=0
            for an1 in el1annotations:#FIXME must use pop list to no repeat
                for an2 in el2annotations:
                    if self.isAnnotationEqual(an1,an2,el1.type):
                        score += self.element_score[an1.qualifier][an2.qualifier]
            #if el1.id.lower() in ['glci','glcx'] or el2.id.lower() in ['glci','glcx']:
            #   if el1.id.lower() in ['glci','glcx'] and el2.id.lower() in ['glci','glcx']:
            #       print 'yepp',el1.id,el2.id, score
            return score

        def compareLocation(self,el1, el2):
            #if not el1.bioQuantity or not el2.bioQuantity:#uncomparable due to missing bioquantities
            #   return 0
            def compare(em1,em2):
                '''internal function for comparing the location, this is needed due to different element types that have different number of locations'''
                score=0
                if em1 and em2:
                    tuple=self.merger.getTuple(em1,em2)
                    if tuple:
                        #print 'found for ',em1,em1.repr,em1.merge_model.repr,em2,em2.repr,em2.merge_model.repr,tuple.pair_score
                        try:
                            return tuple.pair_score[frozenset([em1.mid,em2.mid])]
                        except KeyError:
                            print('should not happen')
                            asdf
                    else:
                        return score
                return score

            if el1.type=='reaction':#reactands and products may be in differnent comapartments
                score=0
                #print 'reactands and products'#DEBUG
                for n in range(2):
                    for x in el1.getLocation()[n]:
                        for y in el2.getLocation()[n]:
                            score += compare(x,y) 
                return score
            else:
                try:
                    return compare(el1.getLocation(),el2.getLocation())
                except NotImplementedError:
                    return 0

        def isAnnotationEqual(self,an1,an2,type):
            '''
            return if 2 annotations can be considered equal
            @type an1: Annotation
            @type an2: Annotation
            @rtype: bool
            '''
            if str(an1)==str(an2):
                #print 'exact match',an1.getName(),type#DEBUG
                return True
            #elif type=='reaction':
                #print 'reaction, comparing directly'
            #   return str(an1)==str(an2)#TODO do a database powered comparison
                #TODO now do something more sophisticated like checking if the reactands are equal
            #elif type == 'compartment':
                #print 'compartment, comparing directly'
            #   return str(an1)==str(an2)
            elif an1.getName() and an2.getName():
                # The original looked up the name of an1 with searchEntries(), which answered
                # from a cache that at this point always held an1 itself, and then compared
                # an1.db+an1.id with an1.db+an2.id. The same decision, written out:
                return an1.id == an2.id
            #elif an2.getName():#this should make not difference since if one is not known we can only do string matching
            #   print an2.getName(),self._ldb.get_id_exact(an2.getName(),type)
            else:#FALLBACK ids not in database comparing directly
                #print 'Id is not in database, comparing directly'
                return str(an1)==str(an2)
        def getScoreMatrix(self,filename):
            '''
            TODO
            triangular to square matrix 
            '''
            try:
                raw_cmp_scr = open(filename,'r').read()
            except IOError:
                raise merge_datamodel.MergeError('missing score matrix file %s' % filename)
            raw_cmp_scr = raw_cmp_scr.split()
            raw_cmp_scr = [int(x) for x in raw_cmp_scr if re.match('-?\\d+',x)]
            #print raw_cmp_scr#DEBUG
            vs=g_bio_qualifier_list
            #vs=[' ','is','hasPart', 'isPartOf', 'isVersionOf', 'hasVersion', 'isHomologTo', 'isDescribedBy','encodes','isEncodedBy', 'unknown']
            #vs=[' ','is','hasPart', 'isPartOf', 'isVersionOf', 'hasVersion', 'isHomologTo', 'isDescribedBy', 'unknown']
            score_matrix={}
            n=0
            try:
                for v in vs:
                    if v not in score_matrix:score_matrix[v]={}
                    for iv in vs[n:]:
                        if iv not in score_matrix:score_matrix[iv]={}
                        score_matrix[v][iv]=raw_cmp_scr.pop(0)
                        score_matrix[iv][v]=score_matrix[v][iv]
                    n+=1
            except IndexError as e:
                raise merge_datamodel.MergeError('Faulty score matrix file %s' % filename)
            return score_matrix

        def updateId(self,libsbml_element,oldId,newId):

            if isinstance(libsbml_element,libsbml.Model):
                elements_w_id=['Species', 'Compartments', 'Reactions', 'Parameters', 'FunctionDefinitions', 'InitialAssignments', 'Rules']
                for etype in elements_w_id:
                    for e in list(getattr(self._merged_element.libsbml_model,'getListOf'+etype)()):
                        self.updateId(e,oldId,newId)

            #-------
            #reaction participants: Reactant, Modifier, Product
            if libsbml_element.getElementName()=='reaction':
                self._updateReactionPointers(libsbml_element,oldId,newId)
            #-------
            #location: compartment, outside
            try:
                self._updateLocationPointers(libsbml_element,oldId,newId)
            except AttributeError:
                pass
            #-------
            #math statement: kineticLaw, initialAssignment,
            try:
                try:
                    libsbml_element=libsbml_element.getKineticLaw()
                except AttributeError:
                    pass
                math = libsbml_element.getMath()
                self._updateMathPointers(math,oldId,newId)
            except AttributeError:
                pass
            #-------

        #--------------------------------------------------------
        def _updateMathPointers(self,math,oldId,newId):

            if math.isFunction():
                if math.getName()==oldId:
                    math.setName(newId)
                for i in range(math.getNumChildren()):
                    childnames=self._updateMathPointers(math.getChild(i),oldId,newId)
            elif math.getNumChildren()==0:
                if math.isName():
                    if math.getName()==oldId:
                        math.setName(newId)
            else:
                #check all children, recursive
                for i in range(math.getNumChildren()):
                    childnames=self._updateMathPointers(math.getChild(i),oldId,newId)

        #--------------------------------------------------------
        def _updateLocationPointers(self,libsbml_element,oldId,newId):
            try:
                if libsbml_element.getOutside()==oldId: libsbml_element.setOutside(newId)
            except AttributeError:
                if libsbml_element.getCompartment()==oldId: libsbml_element.setCompartment(newId)


        #--------------------------------------------------------
        def _updateReactionPointers(self,libsbml_element,oldId,newId):
            types=['Reactant','Modifier','Product']
            for type in types:
                listof=getattr(libsbml_element,'getListOf'+type+'s')()
                for i,elem in enumerate(list(listof)):
                    if elem.getSpecies()==oldId:
                        elem.setSpecies(newId)
        #--------------------------------------------------------

def resetScoreMatrix(filename):
    '''
    create an empty score matrix
    the score matrix should be filled to contain values that represent the relationship between MIRIAM annotations with different or same qualifier
    the score matrix is in human readable format
    '''
    FH=open(filename,'w')
    #vs=(' ','is','hasPart', 'isPartOf', 'isVersionOf', 'hasVersion', 'isHomologTo', 'isDescribedBy','encodes','isEncodedBy', 'unknown')
    vs=['']+g_bio_qualifier_list

    tmpstr = ''.join(['%-14s' for i in range(len(vs))])+'\n'
    FH.write(tmpstr % tuple(vs)) 
    n=len(vs)
    for v in vs:
        if n<len(vs):
            v = '%-14s'%v
            for x in range(len(vs)-1-n):
                v+='%-14s' % ' '
            for x in range(n):
                v+='%-14d' % 1
            FH.write(v+'\n')
        n-=1
    FH.close()

