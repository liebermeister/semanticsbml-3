#!/usr/bin/env python
biomodeldir = "curated/"
import semanticSBML.annotate, libsbml, os
from functools import cmp_to_key

def count_resources_old():
    import popen2
    greps = ["","obo.go","obo.chebi",
             "kegg.compound","kegg.reaction","kegg.pathway",
             "biomodels.db","taxonomy","ec-code",
             "reactome","uniprot","pubmed","obo.fma",
             "icd","omim","pubchem.substance","interpro",
             "pirsf","doi","ensembl","cas",
             "3dmet","intact:EBI","ipi","bind"]
    grep2count = []
    for grep in greps:
        command = "grep urn:miriam:"+grep+" "+biomodeldir+"*.xml | wc -l"
        (stdout,stdin) = popen2.popen2(command)
        count = int(stdout.readlines()[0].strip("\n"))
        grep2count.append([grep,count])


    grep2count.sort(key=cmp_to_key(lambda x,y: y[1]-x[1]))
    print("""\\begin{table}[h]
  \\begin{tabular}{lr}
    Web resource & Number of appearances \\\\\\hline""")
    for [resource,count] in grep2count:
        print("   ",resource,"&",count,"\\\\")
    print(""" \\end{tabular}
  \\caption{Number of appearances of web resources in annotations of compartments, species, and reactions in the BioModels Database.}
  \\label{table:resource_count}
\\end{table}""")

def count_resources(only_compartments_species_reactions=False):
    resourcecount = {}

    for filename in os.listdir(biomodeldir):
        if not filename.endswith(".xml"):
            continue
        d = libsbml.readSBML(biomodeldir+filename)
        ma = semanticSBML.annotate.ModelElementsAnnotations(d.getModel(),True)
        for ea in ma.getElementAnnotations():
            if only_compartments_species_reactions and ea.type not in ['compartment','species','reaction']:
                continue
            for a in ea.getAnnotations():
                if a.db not in resourcecount:
                    resourcecount[a.db] = 0
                resourcecount[a.db] = resourcecount[a.db] + 1

    resourcecount = [list(x) for x in resourcecount.items()]
    resourcecount.sort(key=cmp_to_key(lambda x,y: y[1]-x[1]))
    print("""\\begin{table}[h]
  \\begin{tabular}{lr}
    Web resource & Number of appearances \\\\\\hline""")
    for [resource,count] in resourcecount:
        print("   ",resource,"&",count,"\\\\")
    print(""" \\end{tabular}
  \\caption{Number of appearances of web resources in annotations of compartments, species, and reactions in the BioModels Database.}
  \\label{table:resource_count}
\\end{table}""")

def count_qualifiers(only_compartments_species_reactions=False):
    qualifiercount = {}
    for filename in os.listdir(biomodeldir):
        if not filename.endswith(".xml"):
            continue
        d = libsbml.readSBML(biomodeldir+filename)
        ma = semanticSBML.annotate.ModelElementsAnnotations(d.getModel(),True)
        for ea in ma.getElementAnnotations():
            if only_compartments_species_reactions and ea.type not in ['compartment','species','reaction']:
                continue
            for a in ea.getAnnotations():
                if a.qualifier not in qualifiercount:
                    qualifiercount[a.qualifier] = 0
                qualifiercount[a.qualifier] = qualifiercount[a.qualifier] + 1
    qualifiercount = [[x,y] for x,y in qualifiercount.items()]
    qualifiercount.sort(key=cmp_to_key(lambda x,y: y[1]-x[1]))
    print("""\\begin{table}[h]
     \\begin{tabular}{lr}
       Web resource & Number of apprearances \\\\\\hline""")
    for [qualifier,count] in qualifiercount:
        print("   ",qualifier,"&",count,"\\\\")
    print("""  \\end{tabular}
  \\caption{Number of appearances of qualifiers in annotations in the BioModels Database.}
  \\label{table:qualifier_count}
\\end{table}""")
        

def count_degree_of_annotation(only_compartments_species_reactions=False):
    import semanticSBML.annotate, libsbml, os
    element_count = 0
    number_annotated_elements = 0
    annotation_count = 0
    for filename in os.listdir(biomodeldir):
        if not filename.endswith(".xml"):
            continue
        d = libsbml.readSBML(biomodeldir+filename)
        ma = semanticSBML.annotate.ModelElementsAnnotations(d.getModel(),True)
        for ea in ma.getElementAnnotations():
            if only_compartments_species_reactions and ea.type not in ['compartment','species','reaction']:
                continue
            element_count = element_count + 1
            annotations = ea.getAnnotations()
            if len(annotations) > 0:
                number_annotated_elements = number_annotated_elements + 1
            annotation_count = annotation_count + len(annotations)
    print("Number elements:",element_count)
    print("Annotated elements:",number_annotated_elements)
    print("Number annotations:",annotation_count)
    print("Degree of annotated elements:",float(number_annotated_elements)/float(element_count))
    print("Average number of annotations:",float(annotation_count)/float(element_count))
    print("Average number of annotations per annotated element:",float(annotation_count)/float(number_annotated_elements))

if __name__ == "__main__":
    count_qualifiers(True)
    count_resources(True)
    #count_degree_of_annotation()
    count_degree_of_annotation(True)
