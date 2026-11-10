import libsbml, os
#from sbmlcompartment import sbmlcompartment
from . import annotate
#from config import Config
from  xml.sax import saxutils

import libSBAnnotation.database_new



class CheckResult:
    def __init__(self, name):
        self.name = name
        self.errors = 0
        self.warnings = 0
        self.infos = 0
        self.warning_messages = []
        self.error_messages = []
        self.info_messages = []

    def add_warning(self, message):
        #print 'current message: ', message
        #print 'all messages: ', self.warning_messages

        self.warnings += 1
        self.warning_messages.append(message)

    def add_error(self, message):
        self.errors += 1
        self.error_messages.append(message)

    def add_info(self, message):
        self.infos += 1
        self.info_messages.append(message)
    
    def add_result(self, result_object):
        for lst, adder in (("error_messages", "add_error"),
                           ("warning_messages", "add_warning"),
                           ("info_messages", "add_info")):
            for msg in getattr(result_object, lst):
                getattr(self, adder)(msg)


class Checker:

    def __init__ (self, model = None):
        self.output = None
        if model:
            self.set_model(model)
        self.database = libSBAnnotation.database_new.Meta_DB()

            
    def set_model(self, model):
        self.document = model.get_libsbml_document()
        self.output   = model.get_callback()
        self.model    = self.document.getModel()
        

    def run(self):
        result_list = []
        #result_semantic_check = self.semantic_check()
        #result_list.append(result_semantic_check)
        
        # if file could be read
        #if result_semantic_check.errors == 0:
        for check in [self.conservation_constraint_check,# self.overlap_check,
                      self.semantic_dependency_check, self.annotation_check]:
                      #self.physical_value_check, self.rules_check]:
                      result = check()
                      if result:
                          result_list.append(result)
                    

        return result_list


    def semantic_check( self ):
        """Check for semantic errors"""
    
        self.message('Performing semantic check')
        result = CheckResult('Semantic')

        self.document.checkConsistency()

        #for err_type in ('Error'):
        #    for err_idx in xrange(getattr(self.document, 'getNum' + err_type + 's')()):
        #        result.add_error(getattr(self.document, 'get' + err_type)(err_idx).getMessage())

        return result


    def annotation_check( self ):
        '''this function checks if the model is annotated and writes warnigs if this is not the case'''

        self.message('Checking the model for annotations')
        result = CheckResult('Annotations')
        warnings = 0
        
        try:
                maes=annotate.ModelElementsAnnotations(self.model)
        except Exception as e:
                #result.add_error(str(e))
                maes=annotate.ModelElementsAnnotations(self.model,True)
        
        if maes.getNumNotAnnotatedElements():
            for elem in maes.getElementAnnotations():
                if not elem.isAnnotated():
                    result.add_warning('%s %s: no annotation recognized ' % (elem.type, elem.name or elem.id))
                re=elem.getReadErrors()
                if re: result.add_error(re)
        if result.warnings > 0:
            result.add_info('Some element annotations are missing or could '\
                            + 'not be recognized.\nPlease annotate the model first.'\
                            + ' Further model checks depend on correct element annotations\n')
        else:
            self.message('Finished annotation check: all elements are annotated\n')

        return result

            

    def semantic_dependency_check( self ):#TODO take care of different qualifiers (right now qualifiers have to identical for a match)
        """
        Search for elements that are annotated with the same MIRIAM annotation
        this will be a lot more complicated with qualifiers, must rewrite completely!
        """

        self.message('Performing semantic dependence check')
        result = CheckResult('Duplicate elements')
        
        items_not_checked = {}
        for elt_type, extra_check in (("compartments", None),
                                      ("reactions", None),
                                      ("species", lambda i1, i2: i1.getCompartment() != i2.getCompartment())):
            items_not_checked[elt_type] = 0
            items = getattr(self.model, "getListOf" + elt_type.capitalize())()

            for pos, item1 in enumerate(items):

                try:
                        ae1 = annotate.ElementAnnotations(item1,True)
                except AnnotationCreationError as e:
                            continue
                annotations1 = ae1.getAnnotations()
                if len(annotations1) == 0:
                    items_not_checked[elt_type] += 1
                    continue

                for item2 in items[pos+1:]:
                    ae2 = annotate.ElementAnnotations(item2,True)
                    annotations2 = ae2.getAnnotations()
                    if len(annotations2) == 0:
                        continue
                    if extra_check and extra_check(item1,item2):
                        continue
                    for a in annotations1:
                        if not a in annotations2:
                            break
                        else:
                            #print a
                            #print item1.getName(), ' ', item2.getName()
                            current_warning = ('The %s %s and %s contain identical annotations'
                                         % (elt_type,
                                            item1.getName() or item1.getId(),
                                            item2.getName() or item2.getId()))
                            if not current_warning in result.warning_messages:
                                result.add_warning(current_warning)
                             
        for name in ['compartment', \
                     ['species', 'species'], \
                     'reaction']:
            counter = items_not_checked[type(name) == type("") and (name + "s") or name[1]]
            if counter > 0:
                result.add_warning('%s could not be checked due to missing annotations'
                                   % self.format_counter(counter, name))
                
            self.message('Semantic dependence check finished - %s found\n' %self.format_errors(result.errors))
        
        return result
    
                    
    def overlap_check( self ):#not working, sbmlcomp uses sbmlcompartment to perform the check
        """ check whether one of the compartments is part of another """
        if self.model.getNumCompartments()<2:
            return
        self.message( 'Performing check if compartments are overlapping' )
        result = CheckResult('Overlapping compartments')

        compartmentsNotChecked = 0
        sbmlcomp = sbmlcompartment()
        compartments = self.model.getListOfCompartments()

        for pos, comp1 in enumerate(compartments):
            try:
                ae1=annotate.ElementAnnotations(comp1,True)
            except AnnotationCreationError as e:
                    continue
            annotations1=ae1.getAnnotations()
            go1 = None
            for a in annotations1:
                if a.db=='Gene Ontology':
                    go1=a.id

            if not go1:
                compartmentsNotChecked += 1
                continue
            
            for comp2 in compartments[pos+1:]:

                ae2=annotate.ElementAnnotations(comp2,True)
                annotations2=ae2.getAnnotations()
                go2=None
                for a in annotations2:
                    if a.db=='Gene Ontology':
                        go2=a.id

                if not go2:
                    continue
                        
                if sbmlcomp.isPartOf( go1, go2 ):
                    result.add_warning('The compartments ' + comp1.getId() + ' and ' + comp2.getId() + ' are overlapping')
                 
        if compartmentsNotChecked > 0:
            result.add_warning('%s could not be checked for overlap due to missing Gene Ontology annotation' % self.format_counter(compartmentsNotChecked, 'compartment'))
        self.message('Overlap check finished - %s found\n' % self.format_errors(result.warnings))
        return result


    def conservation_constraint_check( self ):
        """Check for balance of atom numbers in reactions"""
        
        self.message('Performing mass conservation constraint check')
        result = CheckResult('Atom number balance')

        #wsdl = 'http://soap.genome.jp/KEGG.wsdl'
        #serv = WSDL.Proxy( wsdl )

        reactionsNotChecked=0
        reactions = self.model.getListOfReactions()

        for r in reactions:
            formulas = {}

            #self.output.write( 'Starting kegg database lookup' )

            everyElementHasKeggId = True

            # collect empirical formulas of reatants
            for part_type in ['reactants', 'products']:
                formulas[part_type] = []
                for part in getattr(r, 'getListOf' + part_type.capitalize())():
                    species = self.model.getSpecies(part.getSpecies())
                    if species == None:
                        continue
                    kegglink=None
                    ae=annotate.ElementAnnotations(species,True)
                    for a in ae.getAnnotations():
                        if a.db == 'KEGG' or a.db == 'KEGG Compound':
                            kegglink=a.id
                    if not kegglink:
                        everyElementHasKeggId = False
                        continue
                    try:
                        elem_formula = self.database.database_instance.getAnnotationObjectByIdDict({"kegg_compound":kegglink}).dictionary["kegg_formula"]
                    except:
                        everyElementHasKeggId = False
                        continue
                        
                    formulas[part_type].append((part.getSpecies(),
                                                part.getStoichiometry(),
                                                elem_formula))

            if not everyElementHasKeggId:
                reactionsNotChecked += 1
                continue


            atoms = {}
            for part_type in ['reactants', 'products']:
                atoms[part_type] = {}

                for part in formulas[part_type]:
                    multiplier = part[1]
                    formula = part[2]

                    for i in range(len(formula)):

                        if formula[i].isalpha() and i < len(formula)-1:
                            substance = formula[i]
                            indive = 1

                            if formula[i+1].isalpha() and formula[i+1].islower():
                                substance += formula[i+1]
                                i += 1

                            if i < len(formula)-1 and formula[i+1].isdigit():
                                indive = formula[i+1]

                            if indive != 1 and i < len(formula)-2 and formula[i+2].isdigit():
                                indive += formula[i+2]

                            indive = int(indive)

                            if substance in atoms[part_type]:
                                atoms[part_type][substance] += multiplier * indive
                            else:
                                atoms[part_type][substance] = multiplier * indive
                        else:
                            if formula[i].isalpha():
                                substance = formula[i]
                                if substance in atoms[part_type]:
                                    atoms[part_type][substance] += multiplier
                                else:
                                    atoms[part_type][substance] = multiplier


            if atoms['reactants'] != atoms['products']:
                atomsReactants = ''
                atomsProducts = ''
                
                for at in atoms['reactants']:
                    #atomsReactantsInt[at] = int(atoms['reactants'][at])
                    atomsReactants = atomsReactants + str(at) + ' atoms: ' + str(int(atoms['reactants'][at])) + ', '
                    
                for at in atoms['products']:
                    #atomsProductsInt[at] = int(atoms['products'][at])
                    atomsProducts =  atomsProducts + str(at) + ' atoms: ' + str(int(atoms['products'][at])) + ', '

                atomsReactants.rstrip(', ')
                atomsProducts.rstrip(', ')                

                warningStr = 'Reaction [%s] violates atom number conservation\nSubstrate side: %s \nProduct side: %s'%(str(r.getId()), atomsReactants, atomsProducts)

                result.add_warning(warningStr)
        
        if reactionsNotChecked > 0:
            result.add_warning('%s could not be checked for conservation relations due to missing annotations' % self.format_counter(reactionsNotChecked, 'reaction'))

        self.message('Conservation constraint check finished - %s found\n' % self.format_errors(result.errors))
        return result

                        
    def physical_value_check( self ):

        self.message('Performing physical value check')
        result = CheckResult('Physical value')

        for species in self.model.getListOfSpecies():
            if species.getInitialConcentration() < 0:
                result.add_error('species ' + species.getId() + ' is set to a negative concentration')
        for comp in self.model.getListOfCompartments():
            if comp.getSize() < 0:
                result.add_error('compartment ' + comp.getId() + ' is set to a negative size')

        self.message('Physical value check finished - %s found.\n' % self.format_errors(result.errors))
        return result


    def rules_check( self ):

        self.message('Performing check on rules')
        result = CheckResult('Rules')

        # check for algebraic loops
        parameters = []
        rule_variables = []
        for rule in self.model.getListOfRules():

            type = rule.getTypeCode()
            type_string = libsbml.SBMLTypeCode_toString(type)

            if type in [libsbml.SBML_ASSIGNMENT_RULE,
                        libsbml.SBML_RATE_RULE,
                        libsbml.SBML_SPECIES_CONCENTRATION_RULE,
                        libsbml.SBML_COMPARTMENT_VOLUME_RULE,
                        libsbml.SBML_PARAMETER_RULE]:

                # check if the defined entity is being used as a parameter before
                if rule.getVariable() in parameters:
                    result.add_error(('%s variable %s is modified in a rule, ' 
                                      + 'although it is a parameter of a preceding rule')
                                     % (type_string, rule.getVariable()))

                # check if two ruled define the same entity
                if rule.getVariable() in rule_variables:
                    result.add_error(('Identifier %s is set in defined by a rule, '
                                      + 'although is has already been defined by a previous rule')
                                     % rule.getVariable())

                elements = self.getElementsFromMathML(rule.getMath())
                # check if identifier of rule contained in parameters of assignment rule
                if type == libsbml.SBML_ASSIGNMENT_RULE \
                       and rule.getVariable() in elements:
                    result.add_error('%s variable %s is a parameter of the rule by which it is being defined' \
                                     % (type_string, rule.getVariable()))

                parameters.extend(elements)
                rule_variables.extend(rule.getVariable())
        
                # check if rule modifies a species which is already modified in a reaction
                for reaction in self.model.getListOfReactions():
                    for param_getter in [reaction.getListOfReactants, reaction.getListOfProducts]:
                        for param in param_getter():
                            species = param.getSpecies()
                            if rule.getVariable() == species \
                                   and not self.model.getSpecies(species).getBoundaryCondition():
                                result.add_error(('Species %s is being modified by a %s '
                                                  + 'although it is a substrate of  reaction %s without boundaryCondition')
                                                 % (rule.getVariable(), type_string, reaction.getId()))

        self.message('Check on rules finished - %s found\n' % self.format_errors(result.errors))
        return result
                

    def getElementsFromMathML(self, ASTNode):
        """parse <ci> elements out of MathML expression"""
    
        elements = []
        if not ASTNode.isOperator() and ASTNode.getType() == libsbml.AST_NAME:
            elements.append(ASTNode.getName())
        for i in range(ASTNode.getNumChildren()):
            elements.extend(self.getElementsFromMathML(ASTNode.getChild(i)))
        return elements

    # output helper methods
    
    def format_errors(self, errors):
        return self.format_counter(errors, 'error')

    def format_counter(self, counter, names):
        if type(names) == type(str()):
            names = [names, names + 's']
        if counter == 0: return "no " + names[1]
        if counter == 1: return "1 " + names[0]
        return str(counter) + " " + names[1]

    def message(self, string):
        if self.output:
            self.output.write(string)
        
