import os
import pytest

os.environ.setdefault('SEMANTICSBML_QUIET', '1')
# The tests use the annotation database in data/ (built with 'semanticsbml build-database'),
# or the data folder given by the environment variable SEMANTICSBML_DATA.
# tests/test_against_original.py needs the database of February 2011 (see that file).
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.abspath(os.environ.get('SEMANTICSBML_DATA') or os.path.join(ROOT, 'data'))
os.environ['SEMANTICSBML_DATA'] = DATA
HAVE_DATABASE = os.path.exists(os.path.join(DATA, 'db', 'internal_db.main'))
# the database of the original program has no description file
DATABASE_2011 = HAVE_DATABASE and not os.path.exists(os.path.join(DATA, 'db', 'database_sources.txt')) \
    and not os.path.exists(os.path.join(DATA, 'database_sources.txt'))


def pytest_configure(config):
    config.addinivalue_line('markers', 'nodb: test does not need the annotation database')


def pytest_collection_modifyitems(config, items):
    if HAVE_DATABASE:
        return
    skip = pytest.mark.skip(reason='no annotation database in %s; build it with: semanticsbml build-database' % DATA)
    for item in items:
        if 'nodb' not in item.keywords:
            item.add_marker(skip)

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
