"""
Simple functions for comparing whole models by their annotations
(added in the Python 3 port; the computation itself is in cluster.py).

    names, matrix = similarity_matrix(['a.xml', 'b.xml', 'c.xml'])
"""

import os

import libsbml

import semanticSBML.annotate
from . import cluster


def reset():
    """forget the models loaded so far (the similarity code keeps them in global variables)"""
    cluster.height_cache.clear()
    cluster.all_tops.clear()
    cluster.hash_count_dict.clear()
    cluster.pseudo_hash_count = 1
    cluster.DatabaseAnnotation._shared_state.clear()


def model_name(path):
    return os.path.splitext(os.path.basename(path))[0]


def load_models(paths, names=None):
    """
    read SBML files for comparison
    @return: cluster.DatabaseAnnotation, its attribute model_annotations holds one entry per file
    """
    reset()
    names = names or [model_name(p) for p in paths]
    if len(set(names)) != len(names):
        names = list(paths)
    collection = cluster.DatabaseAnnotation()
    for path, name in zip(paths, names):
        document = libsbml.readSBML(path)
        if document.getModel() is None:
            raise ValueError('Could not read %s as SBML' % path)
        annotations = semanticSBML.annotate.ModelElementsAnnotations(document.getModel(), True)
        collection.model_annotations.append(cluster.ModelAnnotation(annotations, name))
    collection.update_information_content()
    return collection


def similarity_matrix(paths, names=None):
    """
    pairwise similarities (between 0 and 1) of annotated SBML models
    @return: (list of model names, list of rows of the symmetric matrix)
    """
    collection = load_models(paths, names)
    models = collection.model_annotations
    n = len(models)
    matrix = [[0.] * n for i in range(n)]
    for i in range(n):
        for j in range(i, n):
            sim = models[i].similarity(models[j])
            matrix[i][j] = matrix[j][i] = sim
    return [m.model_name for m in models], matrix


def cluster_order(collection):
    """order of the models after hierarchical (single linkage) clustering"""
    if len(collection.model_annotations) < 2:
        return [m.model_name for m in collection.model_annotations]
    return [m.model_name for m in collection.cluster().elements]
