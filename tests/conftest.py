import os
import pytest

os.environ.setdefault('SEMANTICSBML_QUIET', '1')

MODELS = os.path.join(os.path.dirname(__file__), 'models')


def model_path(name):
    return os.path.join(MODELS, name)


@pytest.fixture
def hynne_2011():
    """BIOMD0000000061 (Hynne 2001, yeast glycolysis), file as of 2011 with urn:miriam annotations"""
    return model_path('BIOMD0000000061.xml')


@pytest.fixture
def teusink_2011():
    """BIOMD0000000064 (Teusink 2000, yeast glycolysis), file as of 2011 with urn:miriam annotations"""
    return model_path('BIOMD0000000064.xml')


@pytest.fixture
def hynne():
    """BIOMD0000000061, current file from BioModels with identifiers.org annotations"""
    return model_path(os.path.join('current', 'BIOMD0000000061.xml'))


@pytest.fixture
def teusink():
    """BIOMD0000000064, current file from BioModels with identifiers.org annotations"""
    return model_path(os.path.join('current', 'BIOMD0000000064.xml'))
