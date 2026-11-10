"""
Bring the annotation database up to date (added in the Python 3 port).

The original program built its database from files that are no longer offered
(KEGG ftp server, MIRIAM XML export ...).  This module does not rebuild the
database; it takes the existing tables and

- replaces names and relations of ChEBI, Gene Ontology and SBO by those of the
  current ontology files and adds the entries that are new,
- adds KEGG compounds and reactions that are new, from the name lists of the
  KEGG web service (two requests; KEGG data are free for academic use only),
- writes the list of identifiers.org namespaces (successor of MIRIAM).

Links between databases that were made in 2011 (e.g. ChEBI <-> KEGG) are kept
as they are; new entries are linked only where this cannot disturb old links.

    from libSBAnnotation import update
    update.run()                # downloads into <data>/downloads and updates <data>/db
"""

import gzip
import json
import os
import re
import shutil
import time
import urllib.request

from . import config, indexedTables

SOURCES = {
    'chebi': 'https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz',
    'go': 'https://purl.obolibrary.org/obo/go/go-basic.obo',
    'sbo': 'https://raw.githubusercontent.com/EBI-BioModels/SBO/master/SBO_OBO.obo',
    'registry': 'https://registry.api.identifiers.org/resolutionApi/getResolverDataset',
    'kegg_compound': 'https://rest.kegg.jp/list/compound',
    'kegg_reaction': 'https://rest.kegg.jp/list/reaction',
}
FILENAMES = {
    'chebi': 'chebi.obo.gz', 'go': 'go-basic.obo', 'sbo': 'SBO_OBO.obo',
    'registry': 'identifiers_org_registry.json',
    'kegg_compound': 'kegg_compound_list.txt', 'kegg_reaction': 'kegg_reaction_list.txt',
}
NAMESPACE_FILE = 'identifiers_org_namespaces.tsv'
INFO_FILE = 'database_sources.txt'

# relations of ChEBI are given by identifiers of the Relation Ontology today
CHEBI_RELATIONS = {
    'BFO:0000051': 'has_part', 'RO:0000087': 'has_role',
    'RO:0018033': 'is_conjugate_base_of', 'RO:0018034': 'is_conjugate_acid_of',
    'RO:0018036': 'is_tautomer_of', 'RO:0018037': 'is_substituent_group_from',
    'RO:0018038': 'has_functional_parent', 'RO:0018039': 'is_enantiomer_of',
    'RO:0018040': 'has_parent_hydride',
}
# relation types that the similarity code knows (libModelSimilarity/qualifiers.py)
KNOWN_RELATIONS = set(['is_a', 'part_of', 'has_part', 'regulates', 'positively_regulates',
                       'negatively_regulates', 'is_tautomer_of', 'is_enantiomer_of',
                       'is_conjugate_acid_of', 'is_conjugate_base_of', 'has_role',
                       'has_functional_parent', 'is_substituent_group_from', 'has_parent_hydride'])


class UpdateError(Exception):
    pass


# ---------------------------------------------------------------------------
# reading files
# ---------------------------------------------------------------------------

def download(url, path, status=print):
    status('downloading %s' % url)
    request = urllib.request.Request(url, headers={'User-Agent': 'semanticSBML/3'})
    tmp = path + '.part'
    with urllib.request.urlopen(request, timeout=120) as response, open(tmp, 'wb') as out:
        shutil.copyfileobj(response, out)
    if os.path.getsize(tmp) == 0:
        os.remove(tmp)
        raise UpdateError('%s returned an empty file' % url)
    os.replace(tmp, path)


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


_quoted = re.compile(r'"((?:[^"\\]|\\.)*)"')


def parse_obo(path):
    """
    read an OBO file
    @return: (header dictionary, list of terms); a term is a dictionary with the keys id, name,
             namespace, synonyms, alt_ids, is_a, relationships [(type, id)], xrefs [(prefix, id)],
             formula, obsolete, replaced_by, subsets
    """
    header, terms = {}, []
    term, in_term, in_header = None, False, True

    def new_term():
        return {'id': None, 'name': None, 'namespace': None, 'synonyms': [], 'alt_ids': [], 'is_a': [],
                'relationships': [], 'xrefs': [], 'formula': None, 'obsolete': False, 'replaced_by': [],
                'subsets': []}

    with open_text(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('['):
                in_header = False
                if term is not None and term['id']:
                    terms.append(term)
                in_term = line.strip() == '[Term]'
                term = new_term() if in_term else None
                continue
            if ':' not in line:
                continue
            key, value = line.split(':', 1)
            value = value.strip()
            if in_header:
                header.setdefault(key, value)
                continue
            if not in_term:
                continue
            if key == 'id':
                term['id'] = value
            elif key == 'name':
                term['name'] = value
            elif key == 'namespace':
                term['namespace'] = value
            elif key == 'synonym':
                m = _quoted.search(value)
                if m:
                    term['synonyms'].append(m.group(1).replace('\\"', '"').replace('\t', ' '))
            elif key == 'alt_id':
                term['alt_ids'].append(value)
            elif key == 'is_a':
                term['is_a'].append(value.split('!')[0].strip())
            elif key == 'relationship':
                parts = value.split('!')[0].split()
                if len(parts) >= 2:
                    term['relationships'].append((parts[0], parts[1]))
            elif key == 'xref':
                ref = value.split(' ')[0]
                if ':' in ref:
                    prefix, ident = ref.split(':', 1)
                    term['xrefs'].append((prefix, ident))
            elif key == 'property_value':
                if value.startswith('chemrof:generalized_empirical_formula'):
                    m = _quoted.search(value)
                    if m:
                        term['formula'] = m.group(1)
            elif key == 'is_obsolete':
                term['obsolete'] = value.startswith('true')
            elif key == 'replaced_by':
                term['replaced_by'].append(value)
            elif key == 'subset':
                term['subsets'].append(value)
    if term is not None and term['id']:
        terms.append(term)
    return header, terms


# ---------------------------------------------------------------------------
# the update
# ---------------------------------------------------------------------------

class Updater:
    def __init__(self, status=print):
        from . import database_new
        self.status = status
        self.status('loading the annotation database')
        self.db = database_new.Meta_DB().database_instance
        self.report = []

    # -- helpers ------------------------------------------------------------
    def _table(self, internal_name):
        return self.db.__dict__['external_item_' + internal_name]

    def _aid(self, internal_name, ident):
        rows = self._table(internal_name).global_index['data_id'].get(ident)
        if rows:
            return self._table(internal_name).rows[rows[0]][1]
        return None

    def _has_link(self, internal_name, aid):
        return bool(self._table(internal_name).global_index['abstract_item'].get(aid))

    def _filter(self, attribute, keep, change=None):
        """rebuild a table with the rows for which keep(row) is true"""
        old = self.db.__dict__[attribute]
        new = old.__class__(old.database_name) if isinstance(old, indexedTables.ExternalItem) else old.__class__()
        seen = set()
        for row in old.rows:
            if not keep(row):
                continue
            if change is not None:
                row = change(row)
                key = tuple(row)
                if key in seen:
                    continue
                seen.add(key)
            new.insertRow(row)
        self.db.__dict__[attribute] = new
        self.db._relation_cache = None
        return len(old.rows) - len(new.rows)

    def _safe_to_link(self, source, target_db, target_id):
        """
        May a new entry of database `source` be put on the abstract item of target_id?
        Only if the target is unknown, or known and not yet linked to an entry of `source`
        (otherwise links made in 2011 would be rearranged).
        """
        rows = self._table(target_db).global_index['data_id'].get(target_id)
        if not rows:
            return True
        if len(rows) > 1:
            return False
        return not self._has_link(source, self._table(target_db).rows[rows[0]][1])

    def _add(self, database2id, names, item_type, source):
        aid = self.db.addItem(database2id, names, item_type, source)
        self.db.abstract_id_to_annotation_object.clear()
        return aid

    def _begin(self, source):
        """forget relations and the main names that came from `source`"""
        db = self.db
        self._old_names = {}
        for row in db.string_repr_from_database.rows:
            if row[2] == source:
                self._old_names[row[1]] = row[0]
        removed = self._filter('abstract_item_relation', lambda row: row[3] != source)
        self._filter('string_repr_from_database', lambda row: row[2] != source)
        old = self._old_names
        # the old main names stay in the database as ordinary names
        self._filter('string_repr', lambda row: True,
                     lambda row: [row[0], row[1], False, row[3]] if row[2] and old.get(row[1]) == row[0] else row)
        return removed

    def _end(self, source):
        """entries that are gone from the source (obsolete) keep their old main name"""
        db = self.db
        kept = 0
        for aid, name in self._old_names.items():
            rows = db.string_repr.global_index['abstract_item'].get(aid, [])
            if any(db.string_repr.rows[i][2] for i in rows):
                continue
            for i in rows:
                if db.string_repr.rows[i][0] == name:
                    db.string_repr.updateRow(i, [name, aid, True, name.lower()])
                    db.string_repr_from_database.insertRow([name, aid, source, name.lower()])
                    kept += 1
                    break
        # updateRow does not maintain the row hashes
        self._filter('string_repr', lambda row: True)
        return kept

    def _relations(self, source, terms, translate=None):
        db = self.db
        table = self._table(source)
        known = table.global_index['data_id']
        before = len(db.abstract_item_relation.rows)
        for term in terms:
            relations = [['is_a', parent] for parent in term['is_a']]
            for kind, target in term['relationships']:
                if translate is not None:
                    kind = translate.get(kind, kind)
                relations.append([kind, target])
            relations = [r for r in relations if r[0] in KNOWN_RELATIONS and known.get(r[1]) and r[1] != term['id']]
            if relations:
                db.addRelationship(source, term['id'], relations)
        return len(db.abstract_item_relation.rows) - before

    def _old_identifiers(self, source, terms, obsolete, item_type):
        """
        former identifiers (alt_id) and obsolete entries with a successor: as in the original
        program they are separate entries, linked by "current entry is_a former entry"
        """
        db = self.db
        current = set(t['id'] for t in terms)
        n = 0
        pairs = [(t['id'], alt) for t in terms for alt in t['alt_ids']]
        pairs += [(new, t['id']) for t in obsolete for new in t['replaced_by']]
        for new, old in pairs:
            if new not in current or old in current or not self._table(source).global_index['data_id'].get(new):
                continue
            if self._aid(source, old) is None:
                self._add({source: old}, [], item_type(new), source)
            if self._aid(source, old) != self._aid(source, new):
                db.addRelationship(source, new, [['is_a', old]])
                n += 1
        return n

    # -- ChEBI --------------------------------------------------------------
    def update_chebi(self, path, all_entries=False):
        self.status('reading %s' % path)
        header, terms = parse_obo(path)
        if len(terms) < 50000:
            raise UpdateError('%s does not look like the ChEBI ontology (%d entries)' % (path, len(terms)))
        obsolete = [t for t in terms if t['obsolete']]
        terms = [t for t in terms if not t['obsolete'] and t['name']]
        known = self._table('chebi').global_index['data_id']
        if not all_entries:
            # manually curated entries, entries with a KEGG identifier and entries known already
            terms = [t for t in terms if '3:STAR' in t['subsets'] or known.get(t['id'])
                     or any(p == 'kegg.compound' for p, i in t['xrefs'])]
        n_known = sum(1 for t in terms if known.get(t['id']))
        self.status('ChEBI: %d entries, %d of them new' % (len(terms), len(terms) - n_known))
        dropped = self._begin('chebi')
        formulas = self._table('kegg_formula')
        n_linked = n_formula = 0
        for n, term in enumerate(terms):
            if n % 20000 == 0 and n:
                self.status('  %d entries' % n)
            names = [term['name']] + [s for s in term['synonyms'] if s != term['name']]
            database2id = {'chebi': term['id']}
            kegg = sorted(set(i for p, i in term['xrefs'] if p == 'kegg.compound'))
            is_new = not known.get(term['id'])
            if is_new and len(kegg) == 1 and self._safe_to_link('chebi', 'kegg_compound', kegg[0]):
                database2id['kegg_compound'] = kegg[0]
                n_linked += 1
            aid = self._add(database2id, names, 'species', 'chebi')
            if not is_new and len(kegg) == 1 and not self._has_link('kegg_compound', aid) \
                    and self._aid('kegg_compound', kegg[0]) is None:
                self._table('kegg_compound').insertRow([kegg[0], aid, False])
                n_linked += 1
            if term['formula'] and not self._has_link('kegg_formula', aid):
                formulas.insertRow([term['formula'], aid, False])
                n_formula += 1
        n_old = self._old_identifiers('chebi', terms, obsolete, lambda i: 'species')
        added = self._relations('chebi', terms, CHEBI_RELATIONS)
        self._end('chebi')
        version = 'release %s (%s)' % (header.get('data-version', '?'), header.get('date', '?'))
        self.report.append('ChEBI %s: %d entries (%d new), %d relations (before: %d), %d new links to KEGG '
                           'Compound, %d formulas added, %d former identifiers'
                           % (version, len(terms), len(terms) - n_known, added, dropped, n_linked, n_formula, n_old))

    # -- Gene Ontology ------------------------------------------------------
    def update_go(self, path):
        self.status('reading %s' % path)
        header, terms = parse_obo(path)
        if len(terms) < 20000:
            raise UpdateError('%s does not look like the Gene Ontology (%d entries)' % (path, len(terms)))
        obsolete = [t for t in terms if t['obsolete']]
        terms = [t for t in terms if not t['obsolete'] and t['name'] and t['namespace']]
        known = self._table('go').global_index['data_id']
        n_known = sum(1 for t in terms if known.get(t['id']))
        self.status('Gene Ontology: %d entries, %d of them new' % (len(terms), len(terms) - n_known))
        dropped = self._begin('go')
        namespace = {}
        n_linked = 0
        for term in terms:
            namespace[term['id']] = term['namespace']
            names = [term['name']] + [s for s in term['synonyms'] if s != term['name']]
            database2id = {'go': term['id']}
            ec = sorted(set(i for p, i in term['xrefs'] if p == 'EC' and i.count('.') == 3 and '-' not in i))
            if not known.get(term['id']) and len(ec) == 1 and self._safe_to_link('go', 'ec', ec[0]):
                database2id['ec'] = ec[0]
                n_linked += 1
            self._add(database2id, names, term['namespace'], 'go')
        n_old = self._old_identifiers('go', terms, obsolete, lambda i: namespace.get(i, 'biological_process'))
        added = self._relations('go', terms)
        self._end('go')
        self.report.append('Gene Ontology %s: %d entries (%d new), %d relations (before: %d), %d new links to '
                           'EC numbers, %d former identifiers'
                           % (header.get('data-version', '?'), len(terms), len(terms) - n_known, added, dropped,
                              n_linked, n_old))

    # -- SBO ----------------------------------------------------------------
    def update_sbo(self, path):
        self.status('reading %s' % path)
        header, terms = parse_obo(path)
        if len(terms) < 300 or not all(t['id'].startswith('SBO:') for t in terms):
            raise UpdateError('%s does not look like the Systems Biology Ontology' % path)
        terms = [t for t in terms if not t['obsolete'] and t['name']]
        known = self._table('sbo').global_index['data_id']
        n_known = sum(1 for t in terms if known.get(t['id']))
        dropped = self._begin('sbo')
        for term in terms:
            self._add({'sbo': term['id']}, [term['name']], 'sbo_term', 'sbo')
        added = self._relations('sbo', terms)
        self._end('sbo')
        # the ontology file itself is used by libSBAnnotation/sbo.py
        cfg = config.Config()
        for target in (cfg.getpath('resources', 'SBO_resources'),
                       os.path.join(os.path.dirname(cfg.getpath('resources', 'MIRIAM_resources')), 'SBO_OBO.obo')):
            if os.path.abspath(target) != os.path.abspath(path):
                shutil.copyfile(path, target)
        self.report.append('SBO %s: %d entries (%d new), %d relations (before: %d)'
                           % (header.get('data-version', header.get('date', '?')), len(terms),
                              len(terms) - n_known, added, dropped))

    # -- KEGG ---------------------------------------------------------------
    _stoichiometry = re.compile(r'^(\d+|\d*n|\d*m|\((?:\d*[nm])?[+-]?\d*[nm]?\))\s+(.+)$')

    def update_kegg(self, compound_path, reaction_path):
        db = self.db
        compounds = []
        with open_text(compound_path) as f:
            for line in f:
                if '\t' in line:
                    ident, names = line.rstrip('\n').split('\t', 1)
                    ident = ident.replace('cpd:', '')
                    if re.match(r'^C\d{5}$', ident):
                        compounds.append((ident, [n.strip() for n in names.split('; ') if n.strip()]))
        if len(compounds) < 10000:
            raise UpdateError('%s does not look like the KEGG compound list' % compound_path)
        n_new = 0
        first_name, any_name = {}, {}
        for ident, names in compounds:
            if names:
                first_name.setdefault(names[0], []).append(ident)
            for name in names:
                any_name.setdefault(name, []).append(ident)
            aid = self._aid('kegg_compound', ident)
            if aid is None:
                self._add({'kegg_compound': ident}, list(names), 'species', 'kegg_compound')
                n_new += 1
            elif names:
                has_name = any(db.string_repr.rows[i][2] for i in db.string_repr.global_index['abstract_item'].get(aid, []))
                if has_name:
                    db.addName(aid, names)
                else:                       # compound known only as a link from ChEBI so far
                    self._add({'kegg_compound': ident}, list(names), 'species', 'kegg_compound')
        self.report.append('KEGG Compound: %d entries in the list, %d new (names only)' % (len(compounds), n_new))

        def compound(name):
            for index in (first_name, any_name):
                found = index.get(name, [])
                if len(found) == 1:
                    return found[0]
            return None

        def side(text):
            participants = []
            for part in text.split(' + '):
                part = part.strip()
                stoich = 1
                ident = compound(part)
                if ident is None:
                    m = self._stoichiometry.match(part)
                    if m:
                        ident = compound(m.group(2).strip())
                        try:
                            stoich = int(m.group(1))
                        except ValueError:
                            stoich = 1          # n, (n+1) ... are treated as 1, as in the original
                if ident is None:
                    return None
                participants.append([{'kegg_compound': ident}, [], stoich])
            return participants

        n_reactions = n_new = n_equation = 0
        with open_text(reaction_path) as f:
            for line in f:
                if '\t' not in line:
                    continue
                ident, text = line.rstrip('\n').split('\t', 1)
                ident = ident.replace('rn:', '')
                if not re.match(r'^R\d{5}$', ident):
                    continue
                n_reactions += 1
                if self._aid('kegg_reaction', ident) is not None:
                    continue
                parts = [p.strip() for p in text.split('; ')]
                equation = parts[-1] if '<=>' in parts[-1] else None
                names = parts[:-1] if equation else parts
                substrates = products = None
                if equation:
                    left, right = equation.split('<=>', 1)
                    substrates, products = side(left), side(right)
                if substrates is None or products is None:
                    substrates, products = [], []
                    names = names or [equation or ident]
                    if equation and equation not in names:
                        names = names + [equation]
                else:
                    n_equation += 1
                db.addReaction({'kegg_reaction': ident}, substrates, products, [],
                               [[{}, names]] if names else [], 'kegg_reaction')
                db.abstract_id_to_annotation_object.clear()
                n_new += 1
        if n_reactions < 5000:
            raise UpdateError('%s does not look like the KEGG reaction list' % reaction_path)
        self.report.append('KEGG Reaction: %d entries in the list, %d new, %d of them with all participants '
                           'recognised by name' % (n_reactions, n_new, n_equation))

    # -- identifiers.org ----------------------------------------------------
    def update_registry(self, path):
        with open_text(path) as f:
            namespaces = json.load(f)['payload']['namespaces']
        if len(namespaces) < 300:
            raise UpdateError('%s does not look like the identifiers.org registry' % path)
        target = os.path.join(os.path.dirname(config.Config().getpath('resources', 'MIRIAM_resources')), NAMESPACE_FILE)
        with open(target, 'w', encoding='utf-8') as out:
            out.write('#prefix\tname\tpattern\tdeprecated\n')
            for ns in sorted(namespaces, key=lambda n: n['prefix']):
                clean = lambda s: re.sub(r'\s+', ' ', str(s or '')).strip()
                out.write('%s\t%s\t%s\t%s\n' % (clean(ns['prefix']), clean(ns['name']), clean(ns['pattern']),
                                                'deprecated' if ns.get('deprecated') else ''))
        self.report.append('identifiers.org: %d namespaces' % len(namespaces))

    def save(self):
        self.status('writing the database')
        self.db.dump()


def run(download_dir=None, use_existing=False, kegg=True, all_chebi=False, only=None, status=print):
    """
    update the database in the data folder
    @param download_dir: where the source files are kept (default: <data>/downloads)
    @param use_existing: do not download files that are in download_dir already
    @param kegg: also fetch the KEGG name lists (KEGG is free for academic use only)
    @param only: list of parts to update (chebi, go, sbo, kegg, registry); default all
    @return: list of report lines
    """
    data = os.path.dirname(config.Config().getpath('resources', 'MIRIAM_resources'))
    download_dir = download_dir or os.path.join(data, 'downloads')
    os.makedirs(download_dir, exist_ok=True)
    parts = list(only) if only else ['chebi', 'go', 'sbo', 'kegg', 'registry']
    if not kegg and 'kegg' in parts:
        parts.remove('kegg')
    unknown = [p for p in parts if p not in ('chebi', 'go', 'sbo', 'kegg', 'registry')]
    if unknown:
        raise UpdateError('unknown part of the database: %s' % ', '.join(unknown))
    files = {}
    for part in parts:
        for key in (['kegg_compound', 'kegg_reaction'] if part == 'kegg' else [part]):
            files[key] = os.path.join(download_dir, FILENAMES[key])
            if not (use_existing and os.path.exists(files[key])):
                try:
                    download(SOURCES[key], files[key], status)
                except UpdateError:
                    raise
                except Exception as e:
                    raise UpdateError('Could not download %s (%s). Nothing was changed.' % (SOURCES[key], e))
    updater = Updater(status)
    if 'chebi' in parts:
        updater.update_chebi(files['chebi'], all_chebi)
    if 'go' in parts:
        updater.update_go(files['go'])
    if 'sbo' in parts:
        updater.update_sbo(files['sbo'])
    if 'kegg' in parts:
        updater.update_kegg(files['kegg_compound'], files['kegg_reaction'])
    if 'registry' in parts:
        updater.update_registry(files['registry'])
    updater.save()
    with open(os.path.join(data, INFO_FILE), 'a', encoding='utf-8') as out:
        out.write('Updated %s\n' % time.strftime('%Y-%m-%d'))
        for line in updater.report:
            out.write('  ' + line + '\n')
    return updater.report
