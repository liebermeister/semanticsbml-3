import libsbml 
import re 
import sys 
import pickle 
from . import sbmlmisc 
#import database#FIXME must delete this old shit
from . import annotate 
#import nonconform_annotation#FIXME
#from annotationparser import annotationparser
from .annotate import ElementAnnotations, Annotation, database2regex
from SOAPpy import WSDL
from .Model import Model
from libSBAnnotation.config import Config
import libSBAnnotation.database_new, libSBAnnotation.annotate

'''
B{Creation of SBML files}
This module helps to create SBML files by insertion of reaction identifiers and reaction strings
'''

def get_reaction_id_by_substrates( libsbml_reaction, model, database_new):
    return
    # find the kegg reaction id for a reaction with annotated reactants
    possible_reactions=None
    for direction in ['Reactants', 'Products']:
        for species in [ model.getSpecies(ref.getSpecies()) for ref in getattr( libsbml_reaction, 'getListOf'+direction )()]:
            ae = annotate.ElementAnnotations(species)
            try:
                kegg_id = ae.getAnnotations(['KEGG', 'KEGG Compound'])[0].id
            except: 
                return None
            

            reacs = set(libSBAnnotation.database.kegg_get_reaction_from_substrate(kegg_id))
            if not possible_reactions:
                possible_reactions = reacs
            else:
                possible_reactions = possible_reactions.intersection( reacs )
            if len(possible_reactions)==0:
                return None


    if len(possible_reactions)==1:
        return possible_reactions.pop()
    else:
        return None



class KEGG2sbmlError(Exception):
    pass

class kegg2sbml:
    sbo_type2id = { 'activator':           21,\
                    'inhibitor':           20,\
                    'enzyme':              14,\
                    'mrna':                278,\
                    'transcription':       183,\
                    'translation':         184\
                  }
    
    

    """Class to create an SBML file from a list of KEGG reaction IDs """
    def __init__(self,\
                 reactions                      ,\
                 reaction2compartment      =  {},\
                 species2compartment       =  {},\
                 transport_compounds       =  {},\
                 create_modifier_enzymes   =  False,\
                 compartment_diffusion     =  False,\
                 default_metabol_reg       =  False,\
                 default_genreg            =  False,\
                 default_transl_reg        =  False,\
                 ncbi_id                   =  None,\
                 mode                      =  'offline',\
                 model_id                  =  'model',\
                 reaction2gene             =  {},\
                 reaction2modifiers        =  {},\
                 reaction2metabolic_input_functions = {},\
                 reaction2genetic_input_functions ={},\
                 reaction2use              =  {},\
                 reaction2transcr_reg      =  {},\
                 abbrevations              =  {}):


        # not really nice TODO: move all reactionmatching to single class that returns a config object (which can be used by all interfaces)
        kegg_reaction_re = re.compile("(r|R)\\d{5}")
        kegg_pathway_re  = re.compile("[a-zA-Z]{3}\\d{5}")

        for r in reactions:
            if kegg_reaction_re.match(r):
                continue
            elif kegg_pathway_re.match(r):
                k = kegg('online')
                pathway_reacs = k.get_reactions_for_pathway( r )
                reactions.remove(r)
                reactions = reactions + pathway_reacs
                continue
            else:
                continue
                reactions.remove(r)


        self._multiple_compartments = ( reaction2compartment!={} or species2compartment!={} )  # do we have to create more than one compartment ?
        self._cme = create_modifier_enzymes
        self._abbrev=abbrevations
        self._species = []
        self._nameless_reaction_counter=0
        self._reaction2transcr_reg=reaction2transcr_reg
        self._species2compartment = species2compartment
        for sp in self._species2compartment: # remove the disturbing characters
            self._species2compartment[sp] = sbmlmisc.make_valid_xml_id( self._species2compartment[sp] )
        self._reaction2compartment = reaction2compartment
        #self._database = database.Meta_DB()
        self._database_new = libSBAnnotation.database_new.Meta_DB()
        self._errors=[]
        


        # do the logic ...
        if reaction2modifiers!={}:
            default_metabol_reg = True
        if reaction2transcr_reg!={}:
            self._cme=True
        if reaction2transcr_reg!={} or reaction2genetic_input_functions!={} or reaction2transcr_reg!={}:
            default_genreg = True
            default_transl_reg = True
                

        if not mode:
            mode='offline'

        self.kegg = kegg(mode)
        self.re_ec = re.compile("([1-6](\\.\\d+){3})")
        
        if reactions==[]:
            raise KEGG2sbmlError( 'No reactions to create' )

        # doebug
        self._cme = default_genreg = default_transl_reg = default_metabol_reg = False
        # /doebug
        self._cme = True


        
        #if reaction2compartment!={}:
        #    self.standard_compartment = reaction2compartment.values()[0]
        #elif species2compartment!={}:
        #    self.standard_compartment = species2compartment.values()[0]
        #else:
        self.standard_compartment = 'c'  # name of the standard compartment
        self.model, self.abstract_model, self._compartment2name = self._initializeSBMLDoc( reaction2compartment, species2compartment,ncbi_id,model_id)
        
        for r in reactions:
            try: 
                if reaction2use[r]==0: continue # we do not want this stinky reaction
            except: pass

            try:
                comp = self._compartment2name[ self._reaction2compartment[r] ]
            except:  comp = self.standard_compartment
            try:
                gene=reaction2gene[r]
            except:  gene=None

            try:
                [enzyme,reac] = self._make_reaction(r, comp, gene )
            except Exception as e:
                raise
                sys.stderr.write('\n%s\n' %str(e))
                continue
            if not reac:
                continue
            if r in reaction2modifiers and default_metabol_reg:
                inp_func=None
                try:
                    inp_func=reaction2metabolic_input_functions[r]
                except:  pass                    
                self._add_modifiers( reac, reaction2modifiers[r], comp, inp_func )
            if default_transl_reg:  
                [production_reac, decay_reac] = self._add_translational_regulation( reac, enzyme )
            try:
                genetic_input_function = reaction2genetic_input_functions[r]
            except:  genetic_input_function = None
            if default_transl_reg:
                [mrna, mrna_prod, mrna_degr] = self._add_transcriptional_regulation( reac, comp, production_reac, genetic_input_function )

        if compartment_diffusion and transport_compounds=={}:
            self._create_interchange_compartments()

        self._create_transport_compartments( transport_compounds )


    def _write_error(self, text):
        self._errors.append(text+'\n')
        sys.stderr.write(text+'\n')

    
    def _initializeSBMLDoc(self, reaction2compartment, species2compartment, ncbi_id, model_id='id2sbml'):
        abstract_model = Model(None)
        model = abstract_model.get_libsbml_model()
        model_id_new = sbmlmisc.make_valid_xml_id(model_id)
        model.setId(model_id_new)
        model.setName(model_id_new)
        
        if ncbi_id:
            na = nonconform_annotation.nonconform_annotation( model, abstract_model.get_libsbml_document() )
            na.write( 'ncbi_id', ncbi_id )
            #    model.setMetaId( 'Model created by id2SBML, ncbi taxonomy id:%s' %ncbi_id )
        compartment2name = {}
        created_compartments=[]
        
        if self._multiple_compartments:
            for type in [reaction2compartment, species2compartment]:
                for i, id in enumerate( type ):
                    if type[id] in created_compartments:
                        continue
                    comp_id = sbmlmisc.make_valid_xml_id(type[id])
                    if not comp_id:
                        raise KEGG2sbmlError('Invalid compartment ID.')
                    self._create_compartment(model, comp_id, 1, type[id] )
                    compartment2name[ type[id] ] = comp_id
                    created_compartments.append( comp_id )
            self.standard_compartment = model.getCompartment(0).getId()
        else:
            self._create_compartment(model,self.standard_compartment,1,None)
            
        #else:  # create standard compartment
        #    self._create_compartment(model,self.standard_compartment,1,None)

        return (model,abstract_model,compartment2name)
            
            
    def _create_compartment(self, model, compartment, vol, go_annotation):
        sbml_id = sbmlmisc.make_valid_xml_id(compartment)
        c=None
        try:
            c=model.getCompartment(sbml_id)
        except:  pass
        if c:
            return c
        
        comp=model.createCompartment()
        comp.setId( sbml_id )
        comp.setName( sbml_id )
        comp.setVolume(vol)
        try:
            if go_annotation:
                ae=ElementAnnotations(comp)
                ae.addAnnotation( Annotation('Gene Ontology',go_annotation,qualifier='is') )
                comp.setName( self._database_new.getName( db='Gene Ontology', id=go_annotation ) )
        except:
            self._write_error('Error annotating compartment %s' % go_annotation)
        return comp

    
    def _make_reaction( self, id, comp, gene ):
        """ Function that looks up the KEGG database for the given reaction ID and creates this reaction in the SBML model """
        sys.stderr.write('Creating reaction %s' % id)
        no_kegg=False
        id = sbmlmisc.unify_reaction_scheme_delimiter(id)
        lib_id=id
        if '<=>' in id:
            no_kegg=True
            parsed_reaction_info={'equation': id, 'definition': id, 'name': id, 'enzyme': ''}
            lib_id='reaction_'+str(self._nameless_reaction_counter)
            self._nameless_reaction_counter+=1
        else:
            parsed_reaction_info = self.kegg.get_reaction(id)
            if not parsed_reaction_info:
                self._write_error( 'Reaction %s could not be created' %id )
                return [None,None]
        
        #parse the reaction string
        [educts,stoich_educts,products,stoich_products] = sbmlmisc.parse_reaction_scheme(parsed_reaction_info['equation'])
        [name_educts,se,name_products,sp] = sbmlmisc.parse_reaction_scheme(parsed_reaction_info['definition'])    
        reac=self.model.createReaction()
        reac.setId( sbmlmisc.make_valid_xml_id(lib_id) )
        reac.setName( sbmlmisc.make_valid_xml_id(parsed_reaction_info['name']))
        
        if not no_kegg:
            ae=ElementAnnotations(reac)
            ae.addAnnotation( Annotation('KEGG Reaction',lib_id,qualifier='is') )
        
        for pos, e in enumerate(educts):
            spec=self._create_species(e, comp, name_educts[pos] )
            educ = reac.createReactant()
            educ.setSpecies(spec.getId())
            educ.setStoichiometry(stoich_educts[pos])
           
        for pos, p in enumerate(products):
            spec=self._create_species(p, comp, name_products[pos] )
            prod = reac.createProduct()
            prod.setSpecies(spec.getId())
            prod.setStoichiometry(stoich_products[pos])

        anno = get_reaction_id_by_substrates( reac, self.model, self._database_new )#FIXME
        if anno:
            ae = ElementAnnotations(reac)
            ae.addAnnotation( Annotation( 'KEGG Reaction', anno, qualifier='is') )


        enzyme=None
        if self._cme:
            # add enzyme as modifier to the reaction if self._cme flag is set
            try:
                ec = self.re_ec.findall( parsed_reaction_info['enzyme'] )[-1][0]
            except: 
                ec = 'enzyme_'+reac.getId()
            ec_xml = sbmlmisc.make_valid_xml_id(ec)
            #if 'enzyme_'+ec_xml not in self._species:
                
            enzyme = self._create_species(ec, comp, ec_xml, sbo_type='enzyme' )
            mod = reac.createModifier()
            mod.setSpecies( enzyme.getId() )

        if gene:
            na=nonconform_annotation.nonconform_annotation( reac ,self.abstract_model.get_libsbml_document() )
            na.write( 'gene', gene )

                
        sys.stderr.write(" - OK\n")
        return [enzyme,reac]

    def _create_species(self, id, compartment, name, sbo_type=None ):
        
        id_new, anno = self._get_annotation_and_name( id )
        try:
            compartment = self._species2compartment[id]
        except:
            pass
        
        id_new=id_new+'_'+compartment
        annos=None
        if id_new in self._species:
            # see if species already exists by id
            return self.model.getSpecies(id_new)
        
        if len(id)>1:
            if anno:
                sp=sbmlmisc.get_element_by_annotation(self.model,anno)
                if sp and sp.getCompartment()==compartment:
                    return sp

            # see if species already exists by annotations
            new_name = self._remove_compartments_from_name(id)
            if new_name.startswith("C") and len(new_name) == 6:
                annos = {"KEGG Compound":new_name}
            else:
                annos = self._database_new.searchEntries(new_name , 1. )
                annos = annos[0] if len(annos)>0 else {}
            for db in annos:
                try:
                    a=annotate.Annotation( db, annos[db], qualifier='is')
                except:
                    continue
                sp=sbmlmisc.get_element_by_annotation(self.model,a)
                if sp and sp.getCompartment()==compartment:
                    return sp

        spec=self.model.createSpecies()
        spec.setId( sbmlmisc.make_valid_xml_id(id_new) )
        spec.setName(name)
        self._species.append(id_new)
        spec.setCompartment(compartment)
        if sbo_type:
            spec.setSBOTerm( self.sbo_type2id[sbo_type] )
        ae=ElementAnnotations(spec)


        if anno: # if the name is a KEGG id
            ae.addAnnotation( anno )
        if annos:
            for db in annos:
                a = annotate.Annotation( db, annos[db], qualifier='is', suppress_exception = True)
                if db not in libSBAnnotation.annotate.name2uri:
                    continue
                ae.addAnnotation(a, suppress_errors=True)
                #try:
                #    
                #    a = annotate.Annotation( db, annos[db], qualifier='is', suppress_exception = True)
                #    ae.addAnnotation(a)
                #except Exception,e: #pass
                #    print 'error', e
        return spec



    def _add_modifiers( self, reac, modifiers, comp, input_func ):
        for mod,sign in modifiers:
            # at the moment they have to be in kegg compound format
            # TODO: bawawa
            try:
                a=Annotation( 'KEGG Compound', mod, qualifier='is')
                lib_mod=sbmlmisc.get_element_by_annotation( self.model, a )
            except:
                lib_mod=sbmlmisc.get_element_by_id( self.model, sbmlmisc.make_valid_xml_id(mod) )
            if not lib_mod:
                lib_mod = self._create_species(sbmlmisc.make_valid_xml_id(mod), comp, sbmlmisc.make_valid_xml_id(mod) )
            lib_mod_ref = reac.createModifier()
            lib_mod_ref.setSpecies(lib_mod.getId())
            if sign=='+':
                lib_mod_ref.setSBOTerm( self.sbo_type2id['activator'] )  # activator
            else:
                lib_mod_ref.setSBOTerm( self.sbo_type2id['inhibitor'])  # inhibitor
            
        # write input function in annotation tag
        if input_func:
            na=nonconform_annotation.nonconform_annotation( reac ,self.abstract_model.get_libsbml_document() )
            na.write( 'metabolic_input_function', input_func )


    def _add_translational_regulation( self, reac, enzyme ):
        ec_xml=enzyme.getId()
        production_reac=self.model.createReaction()
        enzyme_ref=production_reac.createProduct()
        enzyme_ref.setSpecies(enzyme.getId())
        production_reac.setId('production_'+ec_xml)
        production_reac.setSBOTerm( self.sbo_type2id['translation'] )

        decay_reac=self.model.createReaction()
        decay_reac.createReactant().setSpecies(enzyme.getId())
        decay_reac.setId('degradation_'+ec_xml)

        return [production_reac, decay_reac]
    
    def _add_transcriptional_regulation( self, reac, comp, enzyme_prod_reac, input_func ):
        # create mrna
        mrna = self.model.createSpecies()
        mrna.setId( enzyme_prod_reac.getId()+'_mrna' )
        mrna.setCompartment( comp )
        mrna.setSBOTerm( self.sbo_type2id['mrna'] )
        
        # connect mrna to enzyme production reaction
        mrna_mod_ref = enzyme_prod_reac.createModifier()
        mrna_mod_ref.setSpecies(mrna.getId())
        mrna_mod_ref.setSBOTerm( self.sbo_type2id['activator'] )

        
        # create transcription + degradation of mrna
        mrna_prod = self.model.createReaction()
        mrna_prod.setId('transcription_' + mrna.getId() )
        mrna_prod.setSBOTerm( self.sbo_type2id['transcription'] )
        mrna_ref = mrna_prod.createProduct()
        mrna_ref.setSpecies( mrna.getId() )

        mrna_degr = self.model.createReaction()
        mrna_degr.setId('degradation_' + mrna.getId() )
        mrna_degr.createProduct().setSpecies( mrna.getId() )

        if input_func:
            na=nonconform_annotation.nonconform_annotation( reac ,self.abstract_model.get_libsbml_document() )
            na.write( 'genetic_input_function', input_func )

            tfs = self._parse_genetic_input_function( input_func )
            for tf in tfs:
                id = sbmlmisc.make_valid_xml_id(tf)
                species_tf = self.model.getSpecies( id )
                if not species_tf:
                    species_tf = self.model.createSpecies()
                    species_tf.setCompartment( comp )
                    species_tf.setId(id)
                species_tf_ref = mrna_prod.createModifier()
                species_tf_ref.setSpecies( species_tf.getId() )
            
        return [mrna, mrna_prod, mrna_degr]



    def _parse_genetic_input_function( self, func ):
        delims = ['[', ']', '(', ')', '*', 'AND', 'OR', 'NOT', ',', 'modI', 'modA', '{', '}']
        for term in delims:
            func = func.replace(term, ' ')
        ret=[]
        for elem in func.split():
            if elem:
                ret.append(elem)
        return ret


    # create interchange between compartments (indicated by startflag)
    def _create_interchange_compartments(self):
        species = self.model.getListOfSpecies()
        kegg2comp={}
        kegg_comp2id={}
        for sp in species:
            comp=sp.getCompartment()
            try:
                kegg_id=None
                for db,id in [(a.db,a.id) for a in ElementAnnotations(sp).getAnnotations()]:
                    if db=='KEGG Compound':
                        kegg_id=id
                #kegg_id=annotationparser( sp, True ).getLinks()['KEGG']
                kegg_comp2id[kegg_id+comp]=sp.getId()
            except:
                self._write_error( 'Species %s not annotated' %sp.getId() )
            try:
                if comp not in kegg2comp[kegg_id]:
                    kegg2comp[kegg_id].append(comp)
            except:
                kegg2comp[kegg_id]=[comp]

        for sp in kegg2comp:
            for i, c1 in enumerate(kegg2comp[sp]):
                for c2 in kegg2comp[sp][i+1:]:
                    self._create_interchange_reaction( kegg_comp2id[sp+c1], kegg_comp2id[sp+c2] )



    # create interchange between compartments that has been defined in the input file
    def _create_transport_compartments(self, transport_compounds):
        # get the compartment annotations
        cmp_a={}
        for cmp in self.model.getListOfCompartments():
            try:
                for db,id in [(a.db,a.id) for a in ElementAnnotations(cmp).getAnnotations()]:
                    if db=='GO':
                        cmp_a[cmp.getId()]=id
                #cmp_a[ cmp.getId() ] = annotationparser( cmp, True ).getLinks()['GO']
            except:
                pass

        for comp in  transport_compounds:
            # get the two species
            species1=None
            species2=None
            for sp in self.model.getListOfSpecies():
                if ('KEGG Compound',comp) in[(a.db,a.id) for a in ElementAnnotations(sp).getAnnotations()]:
                #if annotationparser(sp, True).getLinks()['KEGG'] == comp:
                    if cmp_a[sp.getCompartment()] == transport_compounds[comp][0]:
                        species1 = sp
                    elif cmp_a[sp.getCompartment()] == transport_compounds[comp][1]:
                        species2 = sp
            if not species1 or not species2:
                raise KEGG2sbmlError('Species %s not found' %(comp) )

            self._create_interchange_reaction( species1.getId(), species2.getId() )
            
    # create a reaction converting species1 to species2
    def _create_interchange_reaction(self,species1,species2):
        reac=self.model.createReaction()
        reac.setId( sbmlmisc.make_vald_xml_id('Transport_%s_to_%s' %(species1, species2)) )
        for type, id in [('Reactant',species1), ('Product',species2)]:
            sr = getattr(reac, 'create'+type)()
            sr.setSpecies(id)
            sr.setStoichiometry(1)
            

    def __depr_xml_stringer(self,string):
        return sbmlmisc.make_valid_xml_id(string)
        #return string.replace(' ','_').replace('.','_').replace('-','_')


    def writeSBML(self, outfile = None):
        writer = libsbml.SBMLWriter()
        if outfile and outfile != '-':
            sys.stderr.write('Writing SBML to file %s' % outfile)
            #writer.write(self.doc outfile)
            self.abstract_model.save_as( outfile )
        else:
            sys.stderr.write('Writing SBML to stdout')
            print(writer.writeToString( self.abstract_model.get_libsbml_document() ))
        del writer
        sys.stderr.write(' - OK\n')

    def getModel(self):
        return self.abstract_model

    def getErrors(self):
        print('geterr')
        return '\n'.join( self._errors )

    
    def _get_annotation_and_name( self, name ):
        id=sbmlmisc.make_valid_xml_id(name)
        if name in self._abbrev:
            id=self._abbrev[name]
        return ( id, sbmlmisc.get_annotation_for_id(id) )

    def _remove_compartments_from_name( self, name ):
        re_species_comp_replace   =  re.compile( '\\[\\S+\\]' ) 
        return re_species_comp_replace.sub( '', name )
        







class kegg:
    wsdl_path = 'http://soap.genome.jp/KEGG.wsdl'
    def __init__(self, mode):

        if mode not in ['online', 'offline']:
            raise KEGG2sbmlError('Undefined mode %s' %mode)
        self.mode=mode
        

        if mode=='offline':
            self.database = libSBAnnotation.database_new.Meta_DB()
        else:
            try:
                self.kegg_wsdl = WSDL.Proxy( self.wsdl_path )
            except:
                raise KEGG2sbmlError( 'Could not connect to %s . Aborting.' %self.wsdl_path)


    def get_reaction(self,id):
        return eval("self.get_%s('%s')" %(self.mode,id) )


    def get_offline(self,id):        
        dic          = {'entry': '', 'name': '', 'definition': '', 'equation': '', 'enzyme': '', 'other': ''}
        formula      = self.database.getReactionDescription("kegg_reaction",id)
        if not formula:
            sys.stderr.write(' - Could not find reaction %s\n' % id)
            return 
        name         = formula[formula.rfind('[')+2:formula.rfind(']')-1]
        formula      = formula[:formula.rfind('[')-1]            
        dic['equation']   = formula
        dic['definition'] = formula
        dic['name']     = name
        
        return dic
        
    def get_online(self,id):
        reaction_info = self.kegg_wsdl.bget( 'REACTION ' + id )
        if reaction_info=='':
            sys.stderr.write(' - Could not find reaction %s\n' % id)
            return
        
        return self.__parseKeggReaction(reaction_info)

    def get_reactions_for_pathway(self,path):
        reacs = self.kegg_wsdl.get_reactions_by_pathway('path:'+path)
        reacs = [ x[3:] for x in reacs ]
        return reacs
    

    def __parseKeggReaction(self,string):
        """ Sorts the important information from a KEGG request object into a dictionary """        
        dic = {'entry': '', 'name': '', 'definition': '', 'equation': '', 'enzyme': '', 'other': ''}
        last = ''
        for line in string.split('\n'):
            if line.startswith('ENTRY'):
                dic['entry'] = line[12:]
                last = 'entry'
            elif line.startswith('NAME'):
                dic['name'] = line[12:]
                last = 'name'
            elif line.startswith('DEFINITION'):
                dic['definition'] = line[12:]
                last = 'definition'
            elif line.startswith('EQUATION'):
                dic['equation'] = line[12:]
                last = 'equation'
            elif line.startswith('ENZYME'):
                dic['enzyme'] = line[12:]
                last = 'enzyme'
            elif line.startswith('    '):
                dic[last] += line[11:]
            else:
                last='other'
        return dic



class kegg2sbml_parser:

    line_reaction_compartment_re = re.compile("kegg:(r|R)\\d{5}\\s+(G|g)(O|o):\\d{7}")
    line_transport_re = re.compile("Transport\\[kegg:(c|C)\\d{5}\\]\\s+((G|g)(O|o):\\d{7}\\s+){2}")
    line_transport_all_re = re.compile("TransportAll\\s+((G|g)(O|o):\\d{7}\\s+){2}")
    kegg_reaction_re = re.compile("(r|R)\\d{5}")
    kegg_compound_re = re.compile("(c|C)\\d{5}")
    go_id_re = re.compile("([gG][oO]:\\d{7})")
    

    def parse_file_old(self, f ):
        ids = []
        reaction2compartment = {}
        transport_compounds = {}
    

        for line_no, line in enumerate(f):
            line=line.strip()
            if line.startswith('#') or line == "":
                continue
            
            if self.line_reaction_compartment_re.match(line):
                try:
                    reaction = self.kegg_reaction_re.search(line).group()
                except:
                    pass
                try:
                    compartment = self.go_id_re.search(line).group()
                except:
                    pass

            elif self.line_transport_re.match(line):
                try:
                    compound = self.kegg_compound_re.search(line).group()
                    compartments = self.go_id_re.findall(line)
                except:
                    pass
                if len(compartments)==2:
                    transport_compounds.update( {compound: compartments} )
                continue
            else:
                reaction = line


            if not self.kegg_reaction_re.match(reaction):
                raise ParseError("Input file is not in the right format. Parsing error in line %d:\n%s" % (line_no + 1, line))
            ids.append(reaction)

            try:
                reaction2compartment.update({reaction: compartment})
            except:
                pass

        return ( ids, reaction2compartment, transport_compounds )


    def parse_file(self, f):

        ret = { 'reactions':                              [],\
                'reaction2compartment':                   {},\
                'species2compartment':                    {},\
                'transport_compounds':                    {},\
                'reaction2use':                           {},\
                'abbrevations':                           {},\
                'reaction2modifiers':                     {},\
                'reaction2metabolic_input_functions':     {},\
                'reaction2genetic_input_functions':       {},\
                'reaction2gene':                          {},\
                'reaction2transcr_reg':                   {},\
                'dummy':                                  {}\
                
            }

        head2key = { 'reaction':              'reactions',\
                     'compartment':           'reaction2compartment',\
                     'buildreaction':         'reaction2use',\
                     'buildexpression':       'reaction2transcr_reg',\
                     'modifiers':             'reaction2modifiers',\
                     'metabolicregulation':   'reaction2metabolic_input_functions',\
                     'transcriptionalregulation' :'reaction2genetic_input_functions',\
                     'gene':                  'reaction2gene',\
                     'translationalregulation': 'dummy'\
                     }
        header=[]

        for line_no, line in enumerate(f):
            # compile the regexes
            re_species_comp_replace   =  re.compile( '\\[\\S+\\]' ) # part of the species compartment syntax that is not kegg conform and needs to be replaced
            #re_species_comp           =  re.compile( '^|(?<=\+)\S+\[\S+\]' )  # regex for species compartmen
            re_species_comp           = re.compile( '(?:(?<=\\+)|^|(?<=\\>)).*?\\[\\S+\\]' )
            
            line=line.strip()
            if line=='' or line.startswith('#'):
                continue

            if header==[]:
                header = [ x.lower() for x in line.split() ]
                continue


            if len(line.split('\t'))==1:
                if self.kegg_reaction_re.match(line.strip()):
                    # line is a kegg reaction id
                    ret['reactions'].append(line.strip())
                    continue

            
            for i, elem in enumerate(line.split('\t')):
                if elem=='' or elem=='-':
                    continue
                if header[i]=='reaction':
                     for match in re_species_comp.findall(elem):
                         match=match.strip()
                         sep          =  match.find('[')
                         #species      =  match[:sep]
                         species      =  match # we take it all :)
                         compartment  =  match[sep+1:-1]
                         ret['species2compartment'].update({species: compartment})

                     #this does not make sensereaction = re_species_comp_replace.sub( '', elem)
                     ret['reactions'].append( elem )
                     r=elem
                     continue
        

                elif header[i]=='modifiers':
                    mods=[]
                    for term in elem.split(' '):
                        if (term=='') or (not term[0] in  ['+', '-']):
                            continue
                        mods.append( (term[1:], term[0]) )
                    ret['reaction2modifiers'][r]=mods
                    continue
                            

                elif header[i]=='generegulation':
                    ret['reaction2input_functions'][r]=[]
                    for m in elem.split('&'):
                        ret['reaction2input_functions'][r].append(m)
                    continue
                else:
                    try:
                        ret [head2key[header[i]] ].update({r: elem})
                    except:
                        raise Exception('Error parsing file at line %s' %(line_no))
        return ret



            
            
