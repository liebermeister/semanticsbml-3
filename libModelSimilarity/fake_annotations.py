#!/usr/bin/env python

class ModelElementsAnnotations:
    '''
    Represents all L{ElementAnnotations}s (each contains a list of L{Annotation}s)
    functions that operate on multiple L{ElementAnnotations}s
    '''
    def __init__(self,list_of_annotations):
        '''
        Read the SBML file elements that could be annoated.
        Elements are stored as L{ElementAnnotations}s
        @param suppress_errors: 
        @type suppress_errors: bool 
        '''
        self._annotationElements =  [ElementAnnotations(x) for x in list_of_annotations]
        self.libsbml_model = None
        self.libsbml_document = None

    def getElementAnnotations(self):
        '''
        return a list of L{ElementAnnotations} that can contain MIRIAM annotations
        @return: return all annotatable elements of the model
        @rtype: [L{ElementAnnotations}]
        '''
        return self._annotationElements

class ElementAnnotations:
    def __init__(self, annotations):
        self._annotations = annotations
        self.id = annotations[0].id
        self.name = self.id
        self.libsbml_element = FakeLibsbmlElement(self.id)
        self.type = 'species'
    def getAnnotations(self):
        return self._annotations

class FakeLibsbmlElement:
    def __init__(self,id,name=""):
        self.id   = id
        self.name = name or id
    def getId(self):
        return self.id
    def getName(self):
        return self.name
    def getElementName(self):
        return 'species'

class Annotation:
    def __init__(self, db, id, qualifier):
        self.db = db
        self.id = id
        self.qualifier = qualifier

if __name__ == "__main__":
    from . import cluster
    cluster.similarity_measure_version = 4
    cluster.apply_caching = False
    da=cluster.DatabaseAnnotation()
    da.load_from_path("curated/")
    ma = cluster.ModelAnnotation(ModelElementsAnnotations([[Annotation("ChEBI","CHEBI:17234","is")],
                                                           [Annotation("ChEBI","CHEBI:4170","is")],
                                                           [Annotation("ChEBI","CHEBI:28013","is")],
                                                           [Annotation("ChEBI","CHEBI:16108","is")],
                                                           [Annotation("KEGG Compound","C00197","is")],
                                                           [Annotation("KEGG Compound","C00031","is")]
                                                           ]),"current query")
    ma = cluster.ModelAnnotation(ModelElementsAnnotations([[Annotation("KEGG Compound","C00002","is")
                                                           ]]),"current query")
    print(da.compare_model_annotation_object(ma))
    #print da.compare_to_file("bmtest/BIOMD0000000211.xml")
