"""Sbml2Dot converts a SBMLDocument into a .dot file.

functions:
    write_to_dotfile(SBMLDocument, file):
        provides a method to write a dot file from the given SBMLDocument
requires:
    graphviz from www.graphviz.org for displaying the dot file
    """

import libsbml
import os
import re
import subprocess
import time
from . import annotate
from . import  sbo
import libSBAnnotation.config


# alternative coloring for dot file
leftcolor   = '#ffe7be'
rightcolor  = '#c8ddff'
middlecolor = '#b3b0ff'
colormap    = {0:'#ff2200',1:'#bbbbff',2:'#bbffbb',3:'#ffbbbb',4:'#ffff88',5:'#ff88ff',6:'#88ffff',7:'#e5e5ba',8:'#e5bae5',9:'#bae5e5'}

def __get_color_by_number(number):
    return colormap[number % 10]

def overwrite_colors(colored):
    global speciescolor, speciesstyle, modifiercolor, substratearrowcolor, productarrowcolor
    global modifierarrowcolor, reactioncolor, reactionstyle, reactionlabelblank, compartmentcolor
    global irrevsubstratearrowhead, revsubstratearrowhead, revproductarrowtail
    global revsubstratearrowtail
    
    if colored:
        speciescolor = '#a7ebff'
        speciesstyle = 'filled'
        modifiercolor = '#ffda80'
        substratearrowcolor = '#ff8942'
        productarrowcolor = '#ff0342'
        modifierarrowcolor = 'blue'
        reactioncolor = '#d0d0d0'
        reactionstyle = 'filled'
        irrevsubstratearrowhead = 'normal'
        revsubstratearrowhead = 'normal'
        revproductarrowtail = 'normal'
        revsubstratearrowtail = 'normal'
        reactionlabelblank = False
        compartmentcolor           = '#ffda80'
    else:
        speciescolor               = 'black'
        speciesstyle               = 'solid'
        modifiercolor              = 'black'
        substratearrowcolor        = 'black'
        productarrowcolor          = 'black'
        modifierarrowcolor         = 'black'
        reactioncolor              = 'black'
        reactionstyle              = 'solid'
        irrevsubstratearrowhead    = 'none'
        revsubstratearrowhead      = 'none'
        revproductarrowtail        = 'none'
        revsubstratearrowtail      = 'none'
        reactionlabelblank         = True
        compartmentcolor           = 'black'

irreversiblearrowstyle     = 'bold'
speciesshape               = 'ellipse'
reactionshape              = 'box'
modifierarrowstyle         = 'dashed'
revproductarrowhead        = 'normal'
irrevsubstratearrowtail    = 'none'
irrevproductarrowhead      = 'normal'
irrevproductarrowtail      = 'none'
modifierarrowhead          = 'odiamond'
modifierarrowtail          = 'none'

compartmentstyle           = 'rounded'
fontsize                   = '14'
fontname                   = 'Helvetica'
max_number_species         = 50
first_start                = True

def __is_enzyme(anno):
            return 'EC' in [an.db for an in anno.getAnnotations()]
    #return 'EC' in anno.getLinks()

def __write_node_species(file,species,onlymodifier,sourcenumber):
    # Writes a single node to the dot-file and returns the file with a further line
    # Input is a file and a species defined by SBML
    # Output is the name of the node written
    # sourcenumber is the model, the species came from (-1 if the source should be ignored)
    # try getting the label from the annotation
    species_id = species.getId()
    anno = annotate.ElementAnnotations(species,1)

    label = species_id
    keggidstring = ''

    check_viewNames = libSBAnnotation.config.Config().getint('sbml2dot','view_names')

    if check_viewNames == 0:
        if annotate.ElementAnnotations(species,True).getAnnotations():
            label = annotate.ElementAnnotations(species,True).getAnnotations()[0].getName() or species.getId()
            keggidstring =  annotate.ElementAnnotations(species,True).getAnnotations()[0].getURIAction() or ''
            if keggidstring:
                keggidstring='href="%s"' %(keggidstring)
        else:
            label = str(species.getId())
    else:
        label = str(species.getId())

    # getting possible links
    #links = anno.getLinks()

    #if links != None and 'KEGG' in links:
    #    keggidstring = ' href="http://www.genome.jp/dbget-bin/www_bget?cpd:%s"' % links['KEGG']
    color = ''
    # for alternative output
    if sourcenumber != -1:
        color = __get_color_by_number(sourcenumber)
    
    if species_id not in onlymodifier:
        if color == '':
            color = speciescolor
        file.write('\t%s [label=\"%s\" shape=\"%s\" style=\"%s\" fontname=%s, fontsize=\"%s\" color=\"%s\"%s];\n' \
                   % (species_id, str(label), speciesshape, speciesstyle, fontname, fontsize, color, keggidstring))
        
    else:
        # the modifierspecies will be handled later, if it is an enzyme
        if not __is_enzyme(anno):
            if color == '':
                color = modifiercolor
            file.write('\t%s [label=\"%s\" shape=\"%s\" style=\"%s\" fontname=%s, fontsize=\"%s\" color=\"%s\"%s];\n' \
                       % (species_id, str(label), speciesshape, speciesstyle, fontname, fontsize, color, keggidstring))

    return species_id


def __write_node_reaction(file,reaction,onlymodifier,sourcenumber,document,species_blacklist=[]):
    # writes a node for the reaction and binds them to the species
    name = reaction.getName()
    reaction_id = reaction.getId() or name
    reversible = reaction.getReversible()
    # write node for the reaction
    # label is set to the enzyme that modifies it, if the enzyme is recognized as one and is only a modifier to all reactions
    enzymemodifiers = []
    label = ''
    for modifier in reaction.getListOfModifiers():
        # modifier is enzyme and only a modifier
        # get the modifier as a speciesPtr
        modifierbaseelement = None
        for thismodifierbaseelement in document.getModel().getListOfSpecies():
            if thismodifierbaseelement.getId() == modifier.getSpecies():
                modifierbaseelement = thismodifierbaseelement
                break
        
        # get the name from the annotation
        if annotate.ElementAnnotations(modifierbaseelement,True).getAnnotations():
                species = annotate.ElementAnnotations(modifierbaseelement,1).getAnnotations()[0].getName()

        #anno = annotate.ElementAnnotations(modifierbaseelement)
        ##anno = annotationparser(modifierbaseelement, True)
        #species = anno.getNameFromDatabaseLinksOrReturnId()
        #if (modifier.getSpecies() in onlymodifier) and __is_enzyme(anno):
        #    # the modifier will be the label for the reaction
        #    if label == '':
        #        label = species
        #    else:
        #        label += ' or\\n'+species
        #    enzymemodifiers.append(modifier.getSpecies())
    nodeshape = reactionshape
    if label == '' and name != '':
        # no modifying enzyme found
        label = name
    else:
        label = reaction_id

    if reactionlabelblank:
        label = ""
        additional_options = ", width=.15, height=.15"
    else:
        additional_options = ""

    # get the color
    colorstring = ''
    if sourcenumber != -1:
        colorstring += ' color=\"' + __get_color_by_number(sourcenumber) + '\"'
    
    file.write('\t%s [label=\"%s\" shape=\"%s\" style=\"%s\" fontname=%s, fontsize=\"%s\" %s %s];\n' \
        % (reaction_id, label, nodeshape, reactionstyle, fontname, fontsize, colorstring, additional_options) )

    # write bindings to substrates
    for substrate in reaction.getListOfReactants():
        sid = substrate.getSpecies()
        stoichlabel = " taillabel=\""+str(int(substrate.getStoichiometry()))+"\"" if substrate.getStoichiometry() != 1. else ""
        if sid in species_blacklist:
            continue
        
        if not reversible:
            file.write('\t'+sid+' -> '+reaction_id+' [tailport=s headport=n color=\"'+substratearrowcolor+'\" arrowtail="'+irrevsubstratearrowtail+'" arrowhead="'+irrevsubstratearrowhead+'" style="'+irreversiblearrowstyle+'"'+stoichlabel+'];\n')
        else:
            file.write('\t'+sid+' -> '+reaction_id+' [arrowtail="'+revsubstratearrowtail+'" arrowhead="'+revsubstratearrowhead+'" tailport=s headport=n color=\"'+substratearrowcolor+'\"'+stoichlabel+'];\n')
            
    for product in reaction.getListOfProducts():
        sid = product.getSpecies()
        stoichlabel = " headlabel=\""+str(int(product.getStoichiometry()))+"\"" if product.getStoichiometry() != 1. else ""
        if sid in species_blacklist:
            continue

        if not reversible:
            file.write('\t'+reaction_id+' -> '+sid+' [tailport=s headport=n color=\"'+productarrowcolor+'\" arrowtail="'+irrevproductarrowtail+'" arrowhead="'+irrevproductarrowhead+'" style="'+irreversiblearrowstyle+'"'+stoichlabel+'];\n')
        else:
            file.write('\t'+reaction_id+' -> '+sid+' [arrowtail="'+revproductarrowtail+'" arrowhead="'+revproductarrowhead+'" tailport=s headport=n color=\"'+productarrowcolor+'\"'+stoichlabel+'];\n')

    for modifier in reaction.getListOfModifiers():
        sid = modifier.getSpecies()
        if sid in enzymemodifiers and not reactionlabelblank:
            # the species is an enzyme and will be drawn into the reaction node
            pass
        else:
            file.write('\t'+sid+' -> '+reaction_id+' [arrowhead='+modifierarrowhead+' arrowtail='+modifierarrowtail+' style='+modifierarrowstyle+' color=\"'+modifierarrowcolor+'\"];\n')


    
def __get_list_of_models(id2origin):
    origin = []
    for name in id2origin:
        if id2origin[name] not in origin and id2origin[name] != 0:
            origin.append(id2origin[name])
    return origin

def write_legend(file,id2origin):
    """ Making the subgraph that is the legend.
    
    Depending on the output that is actually written,
    the function should decide which legend to write out."""
    
    file.write('\tsubgraph cluster0 {\n')
    file.write('\t\tcolor=\"black\";\n')
    
    if id2origin != None:
        lastspeciesnumber = 0
        lastnode = 'species'
        file.write('\t\tspecies%s [shape=\"%s\" style=\"%s\" color=\"%s\" fontname=%s, fontsize=\"%s\" label=\"species from more\\nthan one model\"];\n' \
                    % (str(lastspeciesnumber), speciesshape, reactionstyle, fontname, fontsize, __get_color_by_number(lastspeciesnumber)) )
                   
        for originnumber in __get_list_of_models(id2origin):
            file.write('\t\tNODE [shape=\"%s\" style=\"%s\" fontname=%s, fontsize=\"%s\" color=\"%s\" label=\"species from\\n\"];\n' \
                       % (speciesshape, speciesstyle, fontname, fontsize, __get_color_by_number(originnumber, originlogger.getOriginlogger().getFilenameByNumber(originnumber)) ) )
            
            file.write('\t\tspecies%s->species%s [color=\"white\"];\n' % (str(lastspeciesnumber), str(lastspeciesnumber+1) ) )
            lastspeciesnumber += 1
            lastnode = 'species' + str(lastspeciesnumber)
    else:
        file.write('\t\tspecies [shape=\"%s\" style=\"%s\" fontname=%s, fontsize=\"%s\" color=\"%s\", label=\"Species\"];\n' \
                   % (speciesshape, speciesstyle, fontname, fontsize, speciescolor) )
        lastnode = 'species'
    values = globals().copy()
    values.update({'lastnode': lastnode})
    values.update({'fontsize': fontsize})
    values.update({'fontname': fontname})
    if reactionlabelblank:
        values.update({'label': ''})
        values.update({'additional_options': ', width=.15, height=.15'})
    else:
        values.update({'label': 'Reaction'})
        values.update({'additional_options': ''})
    file.write("""\t\tNODE [shape=%(reactionshape)s, style=%(reactionstyle)s, fontname=%(fontname)s, fontsize=%(fontsize)s, color="%(reactioncolor)s", label=%(label)s %(additional_options)s];
\t\t%(lastnode)s->reaction [color="white"];\n
\t\t{rank=same; hidden1 [shape=point color="white"];
\t\thidden2 [shape=point color="white"];}
\t\t{rank=same; hidden3 [shape=point color="white"];
\t\thidden4 [shape=point color="white"];}
\t\t{rank=same; hidden5 [shape=point color="white"];
\t\thidden6 [shape=point color="white"];}
\t\t{rank=same; hidden7 [shape=point color="white"];
\t\thidden8 [shape=point color="white"];}
\t\t{rank=same; hidden9 [shape=point color="white"];
\t\thidden10 [shape=point color="white"];}
\t\thidden1->hidden2 [arrowtail=%(revsubstratearrowtail)s arrowhead=%(revsubstratearrowhead)s tailport=e headport=w color="%(substratearrowcolor)s" length=.5 fontname=%(fontname)s headlabel="Substrate" labeldistance=3 labelfontsize=%(fontsize)s];
\t\thidden3->hidden4 [arrowtail=%(revproductarrowtail)s arrowhead=%(revproductarrowhead)s tailport=e headport=w color="%(productarrowcolor)s" length=.5 fontname=%(fontname)s headlabel="Product" labeldistance=3 labelfontsize=%(fontsize)s];
\t\thidden5->hidden6 [tailport=e headport=w color="%(modifierarrowcolor)s" length=.5 fontname=%(fontname)s headlabel="Modifier" labeldistance=3 labelfontsize=%(fontsize)s dir=forward arrowhead=%(modifierarrowhead)s arrowtail=%(modifierarrowtail)s style="%(modifierarrowstyle)s"];
\t\thidden7->hidden8 [arrowtail=none arrowhead=none tailport=e headport=w color="black" length=.5 fontname=%(fontname)s headlabel="Reversible" labeldistance=3 labelfontsize=%(fontsize)s dir=forward];
\t\thidden9->hidden10 [arrowtail=none arrowhead=none tailport=e headport=w color="black" length=.5 fontname=%(fontname)s headlabel="Irreversible" labeldistance=3 labelfontsize=%(fontsize)s dir=forward style="%(irreversiblearrowstyle)s"];
\t\treaction->hidden1 [arrowtail="normal" tailport=s headport=n color="white" length=.05];
\t\thidden1->hidden3 [arrowtail="normal" tailport=s headport=n color="white" length=.05];
\t\thidden3->hidden5 [arrowtail="normal" tailport=s headport=n color="white" length=.05];
\t\thidden5->hidden7 [arrowtail="normal" tailport=s headport=n color="white" length=.05];
\t\thidden7->hidden9 [arrowtail="normal" tailport=s headport=n color="white" length=.05];
\t}\n\n""" % values)

def get_compartments_for_nodes(document):
    """Return a list of compartments and all ids of the nodes which live in it.
    
    Key None says, that the node should not be taken into any cluster."""
    
    compartment2node = {}
    species2compartment = {}
    compartment2node[ None ] = []
    compartment_names = []
    # append all compartments to the list
    for compartment in document.getModel().getListOfCompartments():
        #print compartment.getId()
        compartment2node[ compartment.getId() ] = []
        compartment_names.append(compartment.getName())
    # append all species
    for species in document.getModel().getListOfSpecies():
        #print species.getCompartment()
        compartment2node[ species.getCompartment() ].append( species.getId() )
        species2compartment[ species.getId() ] = species.getCompartment()
    # appending all reactions
    for reaction in document.getModel().getListOfReactions():
        # get all participating species
        participating_compartments = []
        for p in reaction.getListOfReactants():
            compartment = species2compartment[ p.getSpecies() ]
            if compartment not in participating_compartments:
                participating_compartments.append( compartment )
        for p in reaction.getListOfProducts():
            compartment = species2compartment[ p.getSpecies() ]
            if compartment not in participating_compartments:
                participating_compartments.append( compartment )
        for p in reaction.getListOfModifiers():
            compartment = species2compartment[ p.getSpecies() ]
            if compartment not in participating_compartments:
                participating_compartments.append( compartment )
        # if all species are in one compartment
        # we want the reaction to be in this compartment too
        if len(participating_compartments) == 1:
            compartment2node[ participating_compartments[0] ].append( reaction.getId() )
        else:
            # otherwise we are putting it into the none compartment
            compartment2node[ None ].append( reaction.getId() )
    # return the map
    return compartment2node, compartment_names



        
# this function is not used at the moment?
def get_cluster_for_nodes(document,id2origin):
    """Return a list of clusters and all ids of the nodes which live in it.
    
    Key None says, that the node should not be taken into any cluster.
    This will also might modify the originlogger."""
    import sbmldivide
    divider = sbmldivide.sbmldivide(document,id2origin)
    myid2origin = divider.getOriginlogger()
    compartment2node = {}
    species2compartment = {}
    compartment2node[ None ] = []
    # append all compartments to the list
    for compartment in document.getModel().getListOfCompartments():
        compartment2node[ compartment.getId() ] = []
    # append all species
    for species in document.getModel().getListOfSpecies():
        compartment2node[ species.getCompartment() ].append( species.getId() )
        species2compartment[ species.getId() ] = species.getCompartment()
    # appending all reactions
    for reaction in document.getModel().getListOfReactions():
        # get all participating species
        participating_compartments = []
        for p in reaction.getListOfReactants():
            compartment = species2compartment[ p.getSpecies() ]
            if compartment not in participating_compartments:
                participating_compartments.append( compartment )
        for p in reaction.getListOfProducts():
            compartment = species2compartment[ p.getSpecies() ]
            if compartment not in participating_compartments:
                participating_compartments.append( compartment )
        for p in reaction.getListOfModifiers():
            compartment = species2compartment[ p.getSpecies() ]
            if compartment not in participating_compartments:
                participating_compartments.append( compartment )
        # if all species are in one compartment
        # we want the reaction to be in this compartment too
        if len(participating_compartments) == 1:
            compartment2node[ participating_compartments[0] ].append( reaction.getId() )
        else:
            # otherwise we are putting it into the none compartment
            compartment2node[ None ].append( reaction.getId() )
    # return the map
    return (compartment2node,myid2origin)
    
    

def write_to_dotfile(document,filename,id2origin = None,max_reaction_participations = 1000):
    '''
    Writes species and reactions into a .dot file.
    The model has to be present in the document.
    Id2Origin contains a map from which model the elements came at first. id -> 1/2
    '''
    with open(filename,'w') as file:
        file.write(sbml2dot_string(document,max_reaction_participations,id2origin))

def write_to_stream(document,file,max_reaction_participations,id2origin = None):

    max_number_species = libSBAnnotation.config.Config().getint('sbml2dot','max_number')

    if document.getModel().getNumSpecies() > max_number_species:
        raise Exception('Document contains %s species, maximum allowed are %s.'%(document.getModel().getNumSpecies(),max_number_species))

    file.write('digraph SBMLreactions {\n')
    write_legend(file, id2origin )
    
    # get a list of the species, which are only modifiers in all reactions
    onlymodifier = []
    # get all modifiers
    for reaction in document.getModel().getListOfReactions():
        for modifier in reaction.getListOfModifiers():
            if modifier.getSpecies() not in onlymodifier:
                onlymodifier.append(modifier.getSpecies())
    
    # delete the species, that occur as reactants or products
    for reaction in document.getModel().getListOfReactions():
        for modifier in reaction.getListOfReactants():
            if modifier.getSpecies() in onlymodifier:
                onlymodifier.remove(modifier.getSpecies())
        for modifier in reaction.getListOfProducts():
            if modifier.getSpecies() in onlymodifier:
                onlymodifier.remove(modifier.getSpecies())

    # get a list of the species which participate in too many reactions
    species2participations = dict([[s.getId(),0] for s in document.getModel().getListOfSpecies()])
    for reaction in document.getModel().getListOfReactions():
        for reference_list in [reaction.getListOfReactants(),reaction.getListOfProducts(),reaction.getListOfModifiers()]:
            for reference in [r.getSpecies() for r in reference_list]:
                species2participations[reference] = species2participations[reference] + 1
                
    species_blacklist = [s.getId() for s in document.getModel().getListOfSpecies() if species2participations[s.getId()]>max_reaction_participations]

    
    # get a list in which compartment a certain species or reaction lives
    compartment2node, compartment_names = get_compartments_for_nodes(document)
    #(compartment2node,id2origin) = get_cluster_for_nodes(document,id2origin)
    
 
    # writing the code for the compartimentalisation
    compartmentcounter = 0
    for compartment, nodelist in compartment2node.items():
        # writing the code for the subgraph
        if compartment:
            # we have a key None, whose node will not be clustered
            compartmentcounter = compartmentcounter + 1
            file.write('\tsubgraph cluster'+str(compartmentcounter)+' {\n')
            file.write('\t\tfontname="'+fontname+'";\n')
            file.write('\t\tlabel="'+compartment_names[compartmentcounter-1]+'";\n')
            file.write('\t\tstyle=\"'+compartmentstyle+'\";\n')
            file.write('\t\tcolor=\"'+compartmentcolor+'\";\n')

        # first, the nodes have to be written
        # they represent the species
        file.write('NODE [shape='+speciesshape+', color=\"'+speciescolor+'\", style='+speciesstyle+']\n')
        nodes = []
        
        for species in document.getModel().getListOfSpecies():
            # eventually skip this species
            if species.getId() not in nodelist:
                continue
            # leave out species participating in too many reactions
            if species.getId() in species_blacklist:
                continue
            # from which model the species came from
            if id2origin != None:
                sourcenumber = id2origin[species.getId()]
            else:
                sourcenumber = -1
            nodes.append(__write_node_species(file,species,onlymodifier,sourcenumber))
        
        # second, the reaction nodes and their bindings to the species have to be written
        if reactionlabelblank:
            additional_options = ", label=\"\", width=.15, height=.15"
        else:
            additional_options = ""
        file.write('NODE [shape='+reactionshape+',style='+reactionstyle+',color=\"'+reactioncolor+'\"'+additional_options+']\n')
        
        for reaction in document.getModel().getListOfReactions():
            # eventually skip this reaction
            if reaction.getId() not in nodelist:
                continue
            # from which model the reaction came from
            if id2origin != None:
                sourcenumber = id2origin[reaction.getId()]
            else:
                sourcenumber = -1

            __write_node_reaction(file,reaction,onlymodifier,sourcenumber,document,species_blacklist=species_blacklist)
        
        # eventually close the subgraph
        if compartment:
            file.write('}')

    file.write('}')
    first_start == False

    
def write_rules_and_alternatives(ids,model1,model2,filename,sbmlmodelmerger):
    '''
    Function to write ids to nodes, their dependencies to edges and also their alternatives.
    Returns a string containing the SVG Image.
    '''
  
    p = subprocess.Popen([libSBAnnotation.config.Config().getpath('resources','graphvizdot'),'-Tsvg','-Gdpi=72'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE, text=True)
    (file,dot_stdout,dot_stderr) = (p.stdin, p.stdout, p.stderr)
    # at first collect all rules
    allrules = []
    for rule in model1.getListOfRules():
        allrules.append(rule)
    for rule in model2.getListOfRules():
        allrules.append(rule)
    
    #file = open(filename,'w')
    file.write('digraph SBMLrules {\n')
    
    # first, the nodes in ids have to be written
    # they represent the species which cannot be computed due to circles in rule definitions
    file.write('\tNODE [shape=box,style=filled,color=\"red\"]\n')
    for id in ids:
        file.write('\t'+id+' [label=\"'+id+'\"];\n')
    # all other nodes which are already computed should be written in green
    file.write('\tNODE [shape=box,style=filled,color=\"green\"]\n')
    
    # now the selected rules will be displayed
    for dependendid in ids:
        formula = ''
        for rule in allrules:
            try:
                if rule.getVariable() == dependendid and not sbmldelete.isDeleted(rule):
                    if rule.isSetMath():
                        rule.setFormulaFromMath()
                    if isinstance(rule, libsbml.AssignmentRulePtr):
                        formula = rule.getVariable()+' = '+rule.getFormula()
                    elif isinstance(rule, libsbml.RateRulePtr):
                        formula = 'd '+rule.getVariable()+' / dt = '+rule.getFormula()
                    rule.setMathFromFormula()
                    break
            except AttributeError:
                # for algebraic rules
                pass
        # write the reactionnode
        file.write('\tRule2'+dependendid+' [label="'+formula+'" style="normal" color=\"yellow\"];\n')
        file.write('\tRule2'+dependendid+' -> '+dependendid+' [style="normal" tailport="n" headport="s"];\n')
        for targetid in ids[dependendid]:
            file.write('\t'+targetid+' -> Rule2'+dependendid+' [style="normal" tailport="n" headport="s"];\n')
            
    # possible ids of compartments, species and parameters will be collected
    referto = []
    for elt_type in ("Compartments", "Species", "Parameters"):
        for model in model1, model2:
            for elt in getattr(model, 'getListOf' + elt_type)():
                if elt.getId() and elt.getId() not in referto:
                    referto.append(elt.getId())
            
    # now the alternative rules will be displayed
    for dependendid in ids:
        for rule in allrules:
            try:
                if rule.getVariable() == dependendid and sbmldelete.isDeleted(rule):
                    # an alternative rule is found
                    # get the formula
                    formula = ''
                    if rule.isSetMath():
                        rule.setFormulaFromMath()
                    if isinstance(rule, libsbml.AssignmentRulePtr):
                        formula = rule.getVariable()+' = '+rule.getFormula()
                    elif isinstance(rule, libsbml.RateRulePtr):
                        formula = 'd '+rule.getVariable()+' / dt = '+rule.getFormula()
                    rule.setMathFromFormula()
                    # write the rulenode
                    file.write('\tAlternative2'+dependendid+' [label="'+formula+'" style="dotted" color=\"yellow\"];\n')
                    file.write('\tAlternative2'+dependendid+' -> '+dependendid+' [style="dotted" tailport=n headport=s];\n')
                    # search for references in the mathml object
                    for referid in referto:
                        if sbmlmodelmerger.__idoccursinmathml__(rule.getMath(),referid):
                            file.write('\t'+referid+' -> Alternative2'+dependendid+' [style="dotted" tailport=n headport=s"];\n')
            except AttributeError:
                # rule has no variable
                pass
    

    file.write('}')
    file.close()
    outfilename = '.'.join( filename.split('.')[:-1] )+'.gif'
    sleeping_time = 0.001    
    max_sleeping_time = 10
    total_time = 0
    while p.poll() == None and total_time < max_sleeping_time:
        time.sleep(sleeping_time)
        total_time = total_time + sleeping_time     
    svg = dot_stdout.read()
    print(dot_stderr.read())
    dot_stdout.close()
    dot_stderr.close()
    return re.sub('font-size:\\d+.\\d\\d', 'font-size:11px', svg)

def sbml2dot_string(sbml_document,max_reaction_participations = 1000,id2origin = None):
    '''
    @return: the graph of the model in the graphviz dot language
    @rtype: str
    '''
    use_colors = libSBAnnotation.config.Config().getint('sbml2dot','colors')
    overwrite_colors(use_colors == 2)

    class datastore:
        def __init__(self):
            self.store = []

        def write(self,inp):
            self.store.append(inp)

    ds = datastore()
    write_to_stream(sbml_document,ds,max_reaction_participations = max_reaction_participations,id2origin = id2origin)
    return ''.join(ds.store)

def sbml2image(sbml_document,format = 'svg',max_reaction_participations = 1000,id2origin = None):
    '''
    run graphviz on the graph of the model
    @param format: any output format of graphviz (svg, png, pdf ...)
    @return: image data, str for svg, bytes otherwise
    '''
    dot_code = sbml2dot_string(sbml_document,max_reaction_participations,id2origin)
    p = subprocess.Popen([libSBAnnotation.config.Config().getpath('resources','graphvizdot'),'-T'+format,'-Gdpi=72'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,stderr = p.communicate(dot_code.encode('utf-8'))
    if p.returncode != 0:
        raise Exception('Graphviz failed: %s' % stderr.decode('utf-8','replace'))
    if format == 'svg':
        return out.decode('utf-8')
    return out

def sbml2svg(sbml_document,max_reaction_participations = 1000):
    return sbml2image(sbml_document,'svg',max_reaction_participations)

if __name__ == "__main__":
    import libsbml
    d = libsbml.SBMLDocument()
    m = d.createModel()
    m.createCompartment().setId('c')
    for species in 'substrate product modifier'.split():
        s = m.createSpecies()
        s.setId(species)
        s.setCompartment('c')
        
    r = m.createReaction()
    r.setId('reaction1')
    rs = r.createReactant()
    rs.setSpecies('substrate')
    rs.setStoichiometry(2)
    r.createProduct().setSpecies('product')
    r.createModifier().setSpecies('modifier')
    r.setReversible(False)

    r2 = m.createReaction()
    r2.setId('reaction2')
    r2.createReactant().setSpecies('substrate')
    r2.createProduct().setSpecies('product')
    r2.setReversible(True)
   
    svg_code = sbml2svg(d)
    f=open('test.svg','w')
    f.write(svg_code)
    f.close()
