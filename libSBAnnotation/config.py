import os
from os.path import expanduser, join
from configparser import ConfigParser, NoSectionError, NoOptionError

def data_dir():
    """Directory holding the resource files and the annotation database.

    Taken from the environment variable SEMANTICSBML_DATA if set,
    otherwise the 'data' folder next to the source packages.
    (The original tool used ~/.semanticSBML)"""
    d = os.environ.get('SEMANTICSBML_DATA')
    if d:
        return os.path.abspath(expanduser(d))
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')

_OLD_HOME = join('~', '.semanticSBML')

class Config:

    __borg_dict = {'parser': None}
    
    def __init__(self):
        self.__dict__ = self.__borg_dict
        if self.parser is not None:
            return
        self.parser = ConfigParser(interpolation=None)
        self.reload_parser()
        self.parser.read([join(data_dir(), 'config')], encoding='utf-8')

    def reload_parser(self):
        for sec in list(self.parser.sections()):
            self.parser.remove_section(sec)
        self.parser.add_section('resources')
        self.parser.set('resources','graphvizdot','dot')
        self.parser.set('resources','libsbml_resources',join('~','.semanticSBML','libsbmlresources.xml'))
        self.parser.set('resources','brenda_log',join('~','.semanticSBML','brenda_log.txt'))
        self.parser.set('resources','listofresources',join('~','.semanticSBML','listofresources.xml'))
        self.parser.set('resources','element_score_matrix',join('~','.semanticSBML','elementscore.txt'))
        self.parser.set('resources','MIRIAM_resources',join('~','.semanticSBML','MIRIAM_resources.xml'))
        self.parser.set('resources','qualifiers',join('~','.semanticSBML','qualifiers.csv'))
        self.parser.set('resources','SBO_resources',join('~','.semanticSBML/db','SBO_OBO.obo'))
        self.parser.set('resources','balancing_information',join('~','.semanticSBML','information.tsv'))

        self.parser.add_section('database')
        self.parser.set('database', 'database_dir', join('~', '.semanticSBML', 'db') + os.sep)
            
        self.parser.add_section('annotate')
        self.parser.set('annotate', 'search_precision', '0.8')
        self.parser.set('annotate', 'search_precision2', '0.8')            
        # how annotations are written: 'identifiers.org' (today's standard) or 'urn' (urn:miriam:..., as in 2011)
        self.parser.set('annotate', 'uri_style', 'identifiers.org')
        for et in ['species','reaction', 'rule', 'parameter', 'event', 'model', 'compartment', 'kineticLaw']:#TODO find a more generic way to do this
            self.parser.set('annotate', 'annotate_%s_with'%et, '')

        self.parser.add_section('sbml2dot')
        self.parser.set('sbml2dot','max_number','50')
        self.parser.set('sbml2dot','view_names','0')
        self.parser.set('sbml2dot','colors','0')            

        self.parser.add_section('fill')
        self.parser.set('fill','arithmeticgeometric','0')
        self.parser.set('fill','balancing','0')

        #print 'config'
        self.parser.add_section('gui')
        self.parser.set('gui', 'firststart', '1')
        self.parser.set('gui', 'firststart1.9', '1')
        self.parser.set('gui', 'maximized', '0')
        self.parser.set('gui', 'app_x', '10')
        self.parser.set('gui', 'app_y', '10')
        self.parser.set('gui', 'app_w', '640')
        self.parser.set('gui', 'app_h', '420')
        self.parser.set('gui','I2SoutfileOpenPath', os.path.realpath('.' + os.sep))
        self.parser.set('gui','I2SinfileOpenPath', os.path.realpath('.' + os.sep))
        self.parser.set('gui','CopyAnnofileOpenPath', os.path.realpath('.' + os.sep))
        self.parser.set('gui','fileOpenPath', os.path.realpath('.' + os.sep))
        self.parser.set('gui','dirOpenPath', os.path.realpath('.' + os.sep))
        self.parser.set('gui','k2sDirOpenPath', os.path.realpath('.' + os.sep))
        self.parser.set('gui','fileOpenErrorPopup','0')

    def getpath(self, section, option):
        path = self.get(section, option)
        # paths of the form ~/.semanticSBML/... refer to the data directory
        if path.startswith(_OLD_HOME):
            path = data_dir() + path[len(_OLD_HOME):]
        return expanduser(path)

    def get(self, section, option):
        try:
            return self.parser.get(section, option).strip()
        except (NoSectionError,NoOptionError):
            # try rebuilding the config file
            self.reload_parser()
        # if it is not working now, then [section] option will not be in the config file
        return self.parser.get(section, option).strip()
    
    def getint(self, section, option):
        return int(self.get(section, option))
    
    def getfloat(self, section, option):
        return float(self.get(section, option))
    
    def getbool(self, section, option):
        return self.get(section, option) in 'True'

    def set(self, section, option, value):
        self.parser.set(section, option, value)
        
    def write(self):
        if not os.path.isdir(data_dir()):
            os.makedirs(data_dir())
        with open(join(data_dir(), 'config'), 'w', encoding='utf-8') as f:
            self.parser.write(f)


