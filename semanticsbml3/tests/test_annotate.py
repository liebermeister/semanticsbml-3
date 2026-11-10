import libsbml
import pytest

import libSBAnnotation.annotate as lib_annotate
import semanticSBML.annotate as annotate
from semanticSBML import cli


@pytest.mark.parametrize('resource, database, identifier', [
    ('urn:miriam:obo.chebi:CHEBI%3A17234', 'ChEBI', 'CHEBI:17234'),
    ('urn:miriam:kegg.compound:C00031', 'KEGG Compound', 'C00031'),
    ('http://identifiers.org/chebi/CHEBI:17234', 'ChEBI', 'CHEBI:17234'),
    ('https://identifiers.org/chebi/CHEBI:17234', 'ChEBI', 'CHEBI:17234'),
    ('https://identifiers.org/CHEBI:17234', 'ChEBI', 'CHEBI:17234'),
    ('http://identifiers.org/obo.chebi/CHEBI%3A17234', 'ChEBI', 'CHEBI:17234'),
    ('http://identifiers.org/kegg.compound/C00031', 'KEGG Compound', 'C00031'),
    ('http://identifiers.org/go/GO:0005829', 'Gene Ontology', 'GO:0005829'),
    ('https://identifiers.org/GO:0005829', 'Gene Ontology', 'GO:0005829'),
    ('http://identifiers.org/ec-code/2.7.1.1', 'Enzyme Nomenclature', '2.7.1.1'),
    ('http://identifiers.org/taxonomy/4932', 'Taxonomy', '4932'),
    ('https://identifiers.org/uniprot:P12345', 'UniProt', 'P12345'),
])
def test_resource_spellings_are_understood(resource, database, identifier):
    a = lib_annotate.Annotation(urn=resource, qualifier='is')
    assert (a.db, a.id) == (database, identifier)


def test_all_spellings_denote_the_same_annotation():
    spellings = ['urn:miriam:obo.chebi:CHEBI%3A17234', 'http://identifiers.org/chebi/CHEBI:17234',
                 'https://identifiers.org/CHEBI:17234']
    written = set(str(lib_annotate.Annotation(urn=s, qualifier='is')) for s in spellings)
    assert written == {'http://identifiers.org/chebi/CHEBI:17234'}


def test_unknown_namespace_is_kept():
    a = lib_annotate.Annotation(urn='http://identifiers.org/mamo/MAMO_0000046', qualifier='is')
    assert a.id == 'MAMO_0000046'
    assert str(a) == 'http://identifiers.org/mamo/MAMO_0000046'


def test_names_from_database():
    assert lib_annotate.Annotation('ChEBI', 'CHEBI:17234', qualifier='is').getName() == 'glucose'
    assert lib_annotate.Annotation('Taxonomy', '4932', qualifier='is').getName() == 'Saccharomyces cerevisiae'


def test_malformed_resource_is_an_error():
    with pytest.raises(lib_annotate.AnnotationCreationError):
        lib_annotate.Annotation(urn='this is not an annotation', qualifier='is')


@pytest.mark.parametrize('fixture', ['teusink', 'teusink_2011'])
def test_read_annotations_old_and_new_file(fixture, request):
    document = libsbml.readSBML(request.getfixturevalue(fixture))
    species = document.getModel().getSpecies('GLCi')
    found = set((a.db, a.id, a.qualifier) for a in annotate.ElementAnnotations(species).getAnnotations())
    assert ('ChEBI', 'CHEBI:17234', 'is') in found
    assert ('KEGG Compound', 'C00293', 'is') in found


def test_add_and_remove_annotation(teusink, tmp_path):
    document = libsbml.readSBML(teusink)
    element = annotate.ElementAnnotations(document.getModel().getSpecies('GLCi'))
    before = len(element.getAnnotations())
    element.addAnnotation(annotate.Annotation('ChEBI', 'CHEBI:4167', qualifier='isVersionOf'))
    assert len(element.getAnnotations()) == before + 1

    # survives writing and reading
    out = str(tmp_path / 'out.xml')
    libsbml.writeSBMLToFile(document, out)
    document = libsbml.readSBML(out)     # keep the document alive: libSBML frees the model with it
    element = annotate.ElementAnnotations(document.getModel().getSpecies('GLCi'))
    assert ('ChEBI', 'CHEBI:4167', 'isVersionOf') in [(a.db, a.id, a.qualifier) for a in element.getAnnotations()]

    for a in list(element.getAnnotations()):
        element.remAnnotation(a)
    assert element.getAnnotations() == []
    assert not element.isAnnotated()


def test_removing_the_last_annotation_leaves_a_writable_model(teusink):
    # current libSBML crashes when writing a qualifier group without entries
    document = libsbml.readSBML(teusink)
    for element in cli.element_annotations(document):
        element.unsetAnnotations()
    text = libsbml.writeSBMLToString(document)
    assert 'identifiers.org' not in text
    reread = libsbml.readSBMLFromString(text)
    assert reread.getModel().getNumSpecies() == 26


def test_automatic_annotation_by_name():
    document = libsbml.SBMLDocument(2, 4)
    model = document.createModel()
    model.setId('m')
    compartment = model.createCompartment()
    compartment.setId('c')
    species = model.createSpecies()
    species.setId('s1')
    species.setName('Ethanol')
    species.setCompartment('c')
    element = annotate.ElementAnnotations(species)
    element.addAnnotationAutomatic()
    assert ('ChEBI', 'CHEBI:16236') in [(a.db, a.id) for a in element.getAnnotations()]

    # "pyruvate" has two entries in the database (the anion and the acid): no automatic choice
    species.unsetCVTerms()
    species.setName('Pyruvate')
    element = annotate.ElementAnnotations(species)
    with pytest.raises(annotate.AnnotationError):
        element.addAnnotationAutomatic()


def test_search_suggestions(teusink):
    document = libsbml.readSBML(teusink)
    element = annotate.ElementAnnotations(document.getModel().getSpecies('GLCi'))
    names = [a.getName() for a in element.getSuggestions('glucose')]
    assert 'glucose' in names


def test_cli_annotate_roundtrip(teusink, tmp_path, capsys):
    out = str(tmp_path / 'annotated.xml')
    assert cli.main(['annotate', teusink, '--clear', 'GLCi', '--add', 'GLCi', 'ChEBI', 'CHEBI:17234', '-o', out]) == 0
    capsys.readouterr()
    assert cli.main(['annotations', out, '--tsv']) == 0
    lines = [l.split('\t') for l in capsys.readouterr().out.splitlines()]
    glci = [l for l in lines if l[1] == 'GLCi']
    assert [(l[3], l[4], l[5], l[6]) for l in glci] == [('is', 'ChEBI', 'CHEBI:17234', 'glucose')]
