#standard python modules
import re
import os
import sys
import urllib.request, urllib.parse, urllib.error
#own modules
from . import database_new, config, xmlobject

#////////////////////////////////////////////////////
#Error Classes
class AnnotationError(Exception):
    'this is the base class of the errors that can occure in the annotation module'
    def __init__(self,message):
        self.errorMessage = message
    def __repr__(self):
        return self.errorMessage
    def __str__(self):
        return self.errorMessage
class AnnotationDatabaseError(AnnotationError):
    def __init__(self,message):
        self.errorMessage = message
class AnnotationCreationError(AnnotationError):
    def __init__(self,message):
        self.errorMessage = 'Error in Annotation: %s'% message
#////////////////////////////////////////////////////
ldb=None
#////////////////////////////////////////////////////
def redownload(message,option):
    # the original tool tried to download missing resource files here;
    # the download locations no longer exist
    raise AnnotationDatabaseError('Resource file for "%s" is missing or unreadable. %s' % (option, message))
def checkAndGetXml(section,option):
    '''
    check if xml file exists and can be loaded with xmlobject
    inputs are same as get for from config.Config
    @param section: the section in the config
    @type section: str
    @param option: the config entry
    @type option: str
    '''
    #--------------------------------
    #--------------------------------
    filepath= config.Config().getpath(section,option)
    if not os.path.exists(filepath):
        msg = "Missing %s \nTrying to get it ...\n"%(filepath)
        sys.stderr.write(msg)
        #trigger download
        redownload(msg,option)
    try:
        return xmlobject.XMLFile(path=filepath)
    except Exception as e:
        msg = "Error reading %s because: %s \nTrying to get it ...\n"%(filepath,e)
        sys.stderr.write(msg)
        #trigger download
        redownload(msg,option)
#//////////////////////////////////////////////////////
#build lookup dictionaries
#////////////////////////////////////////////////////
#get known BioModels.net qualifiers stuff
#////////////////////////////////////////////////////
bio_qualifiers        = []
model_qualifiers      = []
x = checkAndGetXml('resources','libsbml_resources')
for bqr in x.root.BioQualifier[0].qualifier:
    bio_qualifiers.append(bqr.name)
for mqr in x.root.ModelQualifier[0].qualifier:
    model_qualifiers.append(mqr.name)
#//////////////////////////////////////////////////////
uri2name       = {}#maps ever known (also deprecated) to the db name
name2uri       = {}#this will map the db name to the correct uri
database2elements  = {}
elements2databases = {}#maps libsbml element types to db
name2altURIs       = {}#maps db name to deprecated uris
database2regex     = {}
name2www       = {}#maps db name to weblink in fto

#////////////////////////////////////////////////////
#OLD CODE, Keep to stay compatible, kill once MIRIAM DB is updated
#////////////////////////////////////////////////////
x= checkAndGetXml('resources','listofresources')
for resource in x.root.resource:
    try:
        uri2name['urn:miriam:'+resource.urn] = resource.name
        name2uri[resource.name] = 'urn:miriam:'+resource.urn
    except AttributeError as e:
        pass
        #open('annotate.err.log','a').write("%s\n"%(e))
    uri2name[resource.uri] = resource.name
    if not resource.name in name2uri: name2uri[resource.name] = resource.uri
    name2altURIs[resource.name]=[]
    try:
        database2regex[resource.name]=resource.idPattern
    except AttributeError:
        database2regex[resource.name]=None
    try:
        name2www[resource.name]=resource.action+'%s'
    except AttributeError:
        name2www[resource.name]=None
    try:
        for alt in resource.alternateUris.split(','):
            uri2name[alt] = resource.name
        name2altURIs[resource.name]+=alt.split(',')
    except AttributeError:
        pass
    try:
        database2elements[resource.name]=resource.elements.split()
        for elem in resource.elements.split():
            try:
                elements2databases[elem].add(resource.name)
            except KeyError:
                elements2databases[elem]=set([resource.name])
    except AttributeError:
        pass

#////////////////////////////////////////////////////
#NEW CODE for resources
#////////////////////////////////////////////////////
x = checkAndGetXml('resources','MIRIAM_resources')
for datatype in x.root.datatype:
    try:
        if datatype.deprecated=='true':
            pass
    except AttributeError:
        pass
    name = datatype.name[0]._text

    name2altURIs[name]=[]

    #////////////////////////////////////////////////////
    for xuri in datatype.uris[0].uri:
        try:
            if xuri.deprecated=="true":
                name2altURIs[name].append(xuri._text)
        except AttributeError:
            pass
            #open('annotate.err.log','a').write("%s %s\n"%(name,xuri))
            name2uri[name] = xuri._text
        uri2name[xuri._text] = name 
    #////////////////////////////////////////////////////
    name2www[name]= re.sub('\\$id','%s',str(datatype.resources[0].resource[0].dataEntry[0]))

    #////////////////////////////////////////////////////
    elems=[]
    try: datatype.annotation
    except AttributeError: pass
    else:
        for fo in datatype.annotation[0].format:
            if fo.name=="SBML":
                for elm in fo.elements[0].element:
                    if elm._text=='rule':
                        elems+=['rateRule','assignmentRule','algebraicRule']
                        for e in elems:
                            try:
                                elements2databases[e].add(name)
                            except KeyError:
                                elements2databases[e]=set([name])
                    else:
                        elems.append(elm._text)
                    try:
                        elements2databases[elm._text].add(name)
                    except KeyError:
                        elements2databases[elm._text]=set([name])
            database2elements[name]=elems
    #////////////////////////////////////////////////////
    database2regex[name] = datatype.pattern


#////////////////////////////////////////////////////
# identifiers.org support (added in the Python 3 port)
#
# In 2011 annotations were written as  urn:miriam:<namespace>:<id>.
# Today the standard form is  http(s)://identifiers.org/<namespace>/<id>
# (or the compact form https://identifiers.org/<PREFIX>:<id>).
# Internally all of them are mapped onto the 'urn:miriam:<namespace>' key
# that the resource tables (uri2name ...) of this module use.
#////////////////////////////////////////////////////
# identifiers.org namespace -> namespace used in the 2011 MIRIAM export
namespace_aliases = {
    'chebi': 'obo.chebi', 'go': 'obo.go', 'sbo': 'biomodels.sbo',
    'eco': 'obo.eco', 'pato': 'obo.pato', 'fma': 'obo.fma', 'bto': 'obo.bto',
    'psimod': 'obo.psi-mod', 'mod': 'obo.psi-mod', 'mi': 'obo.mi', 'so': 'obo.so',
    'ncbigene': 'entrez.gene', 'kegg': 'kegg.compound', 'ec': 'ec-code',
    'uniprotkb': 'uniprot', 'ncbitaxon': 'taxonomy', 'biomodels': 'biomodels.db',
    'pubchem.cid': 'pubchem.compound', 'pubchem.sid': 'pubchem.substance',
}
_known_namespaces = set(u[len('urn:miriam:'):] for u in uri2name if u.startswith('urn:miriam:'))

def _read_identifiers_org_namespaces():
    """
    namespaces of the identifiers.org registry that the 2011 MIRIAM table does not have
    (file written by 'semanticsbml update-database'); they get name, identifier pattern and link
    """
    path = os.path.join(os.path.dirname(config.Config().getpath('resources', 'MIRIAM_resources')),
                        'identifiers_org_namespaces.tsv')
    if not os.path.exists(path):
        return
    with open(path, encoding='utf-8') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 3:
                continue
            prefix, name, pattern = fields[0], fields[1], fields[2]
            if prefix in _known_namespaces or prefix in namespace_aliases or not prefix:
                continue
            if not name or name in name2uri:
                name = prefix
                if name in name2uri:
                    continue
            uri = 'urn:miriam:' + prefix
            uri2name[uri] = name
            name2uri[name] = uri
            name2altURIs[name] = []
            database2regex[name] = None     # identifiers are not checked against the pattern: nothing is rejected
            name2www[name] = 'https://identifiers.org/' + prefix + ':%s'
            _known_namespaces.add(prefix)

_read_identifiers_org_namespaces()
namespace_aliases = dict((k, v) for k, v in namespace_aliases.items()
                         if k not in _known_namespaces and v in _known_namespaces)
# 2011 namespace -> preferred identifiers.org namespace (for writing)
_modern_namespace = {}
for _new in ('chebi', 'go', 'sbo', 'eco', 'pato', 'fma', 'bto', 'mod', 'mi', 'so', 'ncbigene'):
    if _new in namespace_aliases:
        _modern_namespace[namespace_aliases[_new]] = _new
# namespaces whose identifiers carry their own prefix (CHEBI:1234, GO:0001234)
_id_prefix = {}
for _name, _regex in list(database2regex.items()):
    _m = re.match(r'^\^([A-Za-z]+):', _regex or '')
    if _m and _name in name2uri and name2uri[_name].startswith('urn:miriam:'):
        _id_prefix[name2uri[_name][len('urn:miriam:'):]] = _m.group(1)

_identifiers_re = re.compile(r'^https?://(?:www\.)?identifiers\.org/(.+)$')

def _normalise_namespace(ns):
    ns = ns.lower()
    if ns not in _known_namespaces and ns in namespace_aliases:
        ns = namespace_aliases[ns]
    return ns

def split_resource(resource):
    """Split an annotation resource string into (database uri, identifier).

    Understands  urn:miriam:<ns>:<id>,  http(s)://identifiers.org/<ns>/<id>,
    http(s)://identifiers.org/<PREFIX>:<id>  and the very old  <uri>#<id>.
    Raises ValueError if the string follows none of these schemes."""
    resource = resource.strip()
    if resource.startswith('urn:miriam:'):
        match = re.match('urn:miriam:([^:]+):(.*)', resource)
        if not match:
            raise ValueError(resource)
        return 'urn:miriam:' + _normalise_namespace(match.group(1)), urllib.parse.unquote(match.group(2))
    match = _identifiers_re.match(resource)
    if match:
        rest = urllib.parse.unquote(match.group(1))
        if '/' in rest and (':' not in rest or rest.index('/') < rest.index(':')):
            ns, id = rest.split('/', 1)                 # .../chebi/CHEBI:17234
        elif ':' in rest:
            ns, id = rest.split(':', 1)                 # .../CHEBI:17234 , .../uniprot:P12345
        else:
            raise ValueError(resource)
        ns = _normalise_namespace(ns)
        prefix = _id_prefix.get(ns)
        if prefix and not id.upper().startswith(prefix.upper() + ':'):
            id = prefix + ':' + id
        if not ns or not id:
            raise ValueError(resource)
        return 'urn:miriam:' + ns, id
    (db, id) = resource.split('#')                      # may raise ValueError
    return db, id

def resource_strings(db_uri, id):
    """All spellings under which the annotation (db_uri,id) may occur in a file."""
    if not db_uri.startswith('urn:miriam:'):
        return ['%s#%s' % (db_uri, urllib.parse.quote(id))]
    ns = db_uri[len('urn:miriam:'):]
    spaces = [ns] + [new for new, old in namespace_aliases.items() if old == ns]
    out = []
    for n in spaces:
        for i in (id, urllib.parse.quote(id)):
            out.append('urn:miriam:%s:%s' % (n, i))
            for scheme in ('http', 'https'):
                out.append('%s://identifiers.org/%s/%s' % (scheme, n, i))
    prefix = _id_prefix.get(ns)
    if prefix and id.upper().startswith(prefix.upper() + ':'):
        for scheme in ('http', 'https'):
            out.append('%s://identifiers.org/%s' % (scheme, id))
    seen = set()
    return [x for x in out if not (x in seen or seen.add(x))]

def format_resource(db_uri, id):
    """Resource string in the configured style ('identifiers.org' or 'urn')."""
    if not db_uri.startswith('urn:miriam:'):
        return '%s#%s' % (db_uri, urllib.parse.quote(id))
    try:
        style = config.Config().get('annotate', 'uri_style')
    except Exception:
        style = 'identifiers.org'
    ns = db_uri[len('urn:miriam:'):]
    if style == 'urn':
        return '%s:%s' % (db_uri, urllib.parse.quote(id))
    return 'http://identifiers.org/%s/%s' % (_modern_namespace.get(ns, ns), id)

#////////////////////////////////////////////////////
#convenience functions
def get_db_from_uri(db):
    try:
        return uri2name[db]
    except KeyError:
        # a namespace that the (2011) MIRIAM table does not know, e.g. from
        # identifiers.org: keep the annotation, using the namespace as database name
        if db.startswith('urn:miriam:') and re.match(r'^[A-Za-z0-9._-]+$', db[len('urn:miriam:'):]):
            name = db[len('urn:miriam:'):]
            if name not in name2uri:
                uri2name[db] = name
                name2uri[name] = db
                name2altURIs[name] = []
                database2regex[name] = None
                name2www[name] = 'https://identifiers.org/' + name + '/%s'
                return name
        raise AnnotationCreationError("Database '%s' Not Supported, please modify %s to quickly add support. If you feel this database should be generally supported suggest it to the <a href='http://www.ebi.ac.uk/miriam/main/mdb?section=submit'>MIRIAM database</a>.<br/><br/>" % (db,config.Config().getpath('resources','listofresources')))#impossible 
#////////////////////////////////////////////////////


class Annotation:

    def __init__(self, db = None, id = None, qualifier_type = 'bio', qualifier = None, suppress_exception = False, urn = None):
        '''
        Annotation class represents an abstract annotation that is a tuple of an identifier that belongs to a database and a qualifier that describes the relation of the anotated element to that database identifier pair.
        the database identifier touple is called link in this context

        @param db: database uri string or name identifier
        @type db: int OR str
        @param id: identifier
        @type id: str
        @param qualifier: qualifier libsbml integer or name string 
        @type qualifier: int OR str
        @param urn: MIRIAM urn e.g. urn:miriam:foo.bar:1234 or http://foo.bar#1234
        @type urn: str
        '''
        if urn:
            try:
                db,id = split_resource(urn)
            except ValueError as e:
                raise AnnotationCreationError('Annotation resource \'%s\' does not follow a known scheme (identifiers.org URL, urn:miriam:DB:ID or URI#ID). ' % urn)
        global ldb
        if not ldb:
            ldb=database_new.Meta_DB()

        self.warning='Not Supported: '
        id=re.sub('%3A',':',id)
        self.setLink(db,id,suppress_exception)
        self.setQualifier(qualifier_type,qualifier)

    def getName(self,try_direct_name=False):
        '''
        @rtype: str
        @return: sting representation of the annotation if in internal db or empty string
        '''
        return ldb.getName(self.db,self.id,try_direct_name) or ''

    def getURIAction(self):
        '''
        @rtype: str
        @return: uri that can be intergrated int a hyperlink to find the refrenced element in the external database
        '''
        try:
            if name2www[self.db]:   return name2www[self.db]%self.id
            else:               return None
        except KeyError:
            return None

    def __eq__(self,other):
        '''
        equality operator
        '''
        eq_elems=['db','id','qualifier','qualifier_type']
        if self.db.startswith('KEGG') and other.db.startswith('KEGG') and (self.db=='KEGG' or other.db=='KEGG'): # consider KEGG compound and KEGG as equal
            eq_elems.remove('db')
        return not False in [ getattr(self,x)==getattr(other,x) for x in eq_elems ]

    def safe_str(self):
        try:
            return str(self)
        except AnnotationCreationError:
            return "unknown"

    def __str__(self,suppress_exception=False):
        '''
        returns the annotation resource as specified by the proposed MIRIAM annotatation standard
        @return: MIRIAM annotation string "urn:miriam:database[.sub]:id"
        @rtype: str
        '''
        if self.db.find(self.warning)!=-1:
            return str(self.db[len(self.warning):]+'#'+self.id)
        elif self.db == 'SBO':
            return self.db
        elif not self.db in name2uri:#check so we can merge our old code with listofresources.xml
            if not suppress_exception:
                raise AnnotationCreationError('str Database \'%s\' yet unknown, please extend %s' % (self.db,config.Config().getpath('resources','listofresources')))
            else:
                return "unknown"
        else:
            return format_resource(name2uri[self.db], self.id)

    def getResourceStrings(self):
        '''
        all spellings (urn:miriam, identifiers.org ...) under which this annotation may be stored in an SBML file
        @rtype: [str]
        '''
        if self.db in name2uri:
            return resource_strings(name2uri[self.db], self.id)
        return []


    def getUrlAnnotation(self):
        '''should not be used anymore'''
        #if self.db.find(self.warning)!=-1:
        #   return str(self.db[len(self.warning):]+'#'+self.id)
        #elif not name2url.has_key(self.db):#check so we can merge our old code with listofresources.xml
        #   raise AnnotationCreationError('url Database %s yet unknown, please extend %s' % (self.db,config.Config().get('resources','listofresources')))
        #else:
        #   return str(name2url[self.db]+'#'+self.id)#old url style now outdated
        return 'DO NOT USE THIS FUNCTION'

    def setQualifier(self,qualifier_type,qualifier):
        '''
        set the qualifier of the annotation
        the list of qualifiers used can be found at http://biomodels.net/index.php?s=Qualifiers
        libsbml encoded them in numbers between 0 and 7
        there are 3 types in current libsbml:
            - qualifier
            - bio qualifier
            - model qualifier
        this still has to be implemented
        @param qualifier: qualifier string or libsbml integer
        @type qualifier: str OR int
        '''
        #print 'qua ',qualifier
        #print qualifier_type
        qtypes = {0:'model',1:'bio',2:'unknown'} #TODO move to xml
        name_convert = {'bqmodel':'model',
                        'bqbiol':'bio'}
        if isinstance(qualifier_type,str) or isinstance(qualifier_type, str):
            self.qualifier_type = qualifier_type.lower()
            if self.qualifier_type in name_convert:
                self.qualifier_type = name_convert[self.qualifier_type]
        elif isinstance(qualifier_type,int):
            try:
                self.qualifier_type = qtypes[qualifier_type]
            except KeyError:
                #print 
                raise AnnotationCreationError('Unkown Qualifier Type, Allowed are: Model, Bio, Unknown. Used: %s'%qualifier_type)

        if self.qualifier_type == 'model':
            if qualifier == None: qualifier = 'unknown'
            if isinstance(qualifier, int):
                raise AnnotationCreationError('Inserted an int value as model qualifier, this is not allowed.')
            if qualifier in model_qualifiers:
                self.qualifier = qualifier
            else:
                raise AnnotationCreationError('Unknown Qualifier \'%s\', go to http://biomodels.net/index.php?s=Qualifiers'%qualifier)
        elif self.qualifier_type == 'bio':
            if qualifier == None: qualifier = 'unknown'
            if isinstance(qualifier, int):
                raise AnnotationCreationError('Inserted an int value as bio qualifier, this is not allowed.')
            if qualifier in bio_qualifiers:
                self.qualifier = qualifier
            else:
                raise AnnotationCreationError('Unknown Qualifier \'%s\', go to http://biomodels.net/index.php?s=Qualifiers'%qualifier)
        elif self.qualifier_type == 'unknown':
            #sys.strerror.write('unknown type',self.qualifier_type)#DEBUG
            self.qualifier = 'unknown'
        else:
            print(qualifier_type)

        #print bio_qualifier_type2name.values()
        #print 'qualifier set'#DEBUG

    def setLink(self,db,id,suppress_exception=False):
        '''
        set the database and the identifier of the Annotation class, convert the database to the main format, url,urn and name are accepted
        @param db: database urn,url string or name identifier
        @type db: int OR str
        @param id: identifier
        @type id: str
        @param qualifier: qualifier libsbml integer or name string 
        @type qualifier: int OR str
        @param suppress_exception: suppress database not supported exception TODO see if this is still needed! 
        @type suppress_exception: bool
        '''
        if not db or not id:
            if not suppress_exception:
                raise AnnotationCreationError('No Database or ID set DB: \' %s\' ID: \'%s\' '%(db,id))
        elif db.find('http')>-1 or db.find('urn:')>-1:
            odb=get_db_from_uri(db)
            if not odb:
                if not suppress_exception:
                    raise AnnotationCreationError('Database %s Not Supported' % db)
            else:
                db=odb
        elif db in name2uri:
            pass
        else:
            if not suppress_exception:
                raise AnnotationCreationError('Database \'%s\' Not Supported' % db)

        self.checkIdPattern(db,id,suppress_exception)
        self.db = db
        self.id = id

    def checkIdPattern(self,db,id,suppress_exception=False):
        '''
        check if the database can be found and if a regular expression matches the identifier of the database 
        the regular expressions and databases are provided by the listofresources.xml wich was provided by the nice folks of BioModels
        @param db: database name
        @type db: str
        @param id: identifier
        @type id: str
        @param suppress_exception: suppress database not supported exception TODO see if this is still needed! 
        @type suppress_exception: bool
        @return: retrns true if the check was succesfull or throws L{AnnotationError} if the check failed
        @rtype: bool
        '''
        if db in database2regex:
            if database2regex[db]:
                if not re.match(database2regex[db],id):
                    if not suppress_exception:
                        raise AnnotationCreationError('%s identifier \'%s\' did not match the reqular expression scheme %s' % (db, id,database2regex[db]))
                    else:
                        return False
                else:
                    return True


