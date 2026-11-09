#!/usr/bin/env python

import semanticSBML.annotate
import libSBAnnotation.database_new
import libSBAnnotation.mydifflib

d = libSBAnnotation.database_new.Meta_DB()

class InChI:
    def __init__(self, doc):
        self.doc = doc
        self.mea = semanticSBML.annotate.ModelElementsAnnotations(self.doc.getModel())
        self.id2inchi = {}
        for ea in self.mea.getElementAnnotations():
            inchi = None
            for a in ea.getAnnotations():
                aid = d.database_instance.__getAbstractIdByIdDict__({a.db:a.id})
                if aid is not None:
                    names = d.database_instance.__getAllNamesByAbstractId__(aid)
                    inchi_names = [name[6:] for name in names if name.startswith('InChI=')]
                    if len(inchi_names)>0:
                        inchi = inchi_names[0]
            if inchi is not None:
                self.id2inchi[ea.id] = inchi
        self._sim_dict = None
    def sim_dict_append(self, keyi, keyj, sim):
        for key in [keyi, keyj]:
            if key not in self._sim_dict:
                self._sim_dict[key] = {}
        self._sim_dict[keyi][keyj] = sim
        self._sim_dict[keyj][keyi] = sim
    def _compute_sim_dict(self):
        sm = libSBAnnotation.mydifflib.StringMatcher()
        self._sim_dict = {}
        keys = list(self.id2inchi.keys())
        for i in range(len(keys)):
            keyi = keys[i]
            sm.set_seq1(self.id2inchi[keyi])
            for j in range(i+1, len(keys)):
                keyj = keys[j]
                sm.set_seq2(self.id2inchi[keyj])
                self.sim_dict_append(keyi, keyj, sm.ratio())
    def get_closest_matches(self):
        if self._sim_dict is None:
            self._compute_sim_dict()
        keys = sorted(self.id2inchi.keys())
        for key in keys:
            sims = self._sim_dict[key]
            [key2, sim] = sorted([(x, sims[x]) for x in sims], key=lambda x:x[1], reverse=True)[0]
            print(key, key2, sim)
    def get_closest_matches_for_reactions(self):
        l = []
        if self._sim_dict is None:
            self._compute_sim_dict()
        for reaction in self.doc.getModel().getListOfReactions():
            substrates = [s.getSpecies() for s in reaction.getListOfReactants() if s.getSpecies() in self._sim_dict]
            products = [p.getSpecies() for p in reaction.getListOfProducts() if p.getSpecies() in self._sim_dict]
            modifiers = [m.getSpecies() for m in reaction.getListOfModifiers() if m.getSpecies() in self._sim_dict]
            for [from_list, to_list, name] in [[substrates, products, 'substrate'], [products, substrates, 'product']]:
                for _from in from_list:
                    [_to, sim] = sorted([(x, self._sim_dict[_from][x]) for x in to_list], key=lambda x:x[1], reverse=True)[0]
                    print(reaction.getId(), name, _from, _to, sim)
                    l.append([reaction.getId(), name, _from, _to, sim])
            for modifier in modifiers:
                [best_sub, sub_sim] = sorted([(x, self._sim_dict[modifier][x]) for x in substrates], key=lambda x:x[1], reverse=True)[0]
                [best_pro, pro_sim] = sorted([(x, self._sim_dict[modifier][x]) for x in products], key=lambda x:x[1], reverse=True)[0]
                all_modifiers = [x for x in modifiers if x in self._sim_dict[modifier] and self._sim_dict[modifier][x] >= 0.6]
                print(reaction.getId(), 'modifier', modifier, [best_sub, best_pro]+all_modifiers)
                l.append([reaction.getId(), 'modifier', modifier]+[best_sub, best_pro]+all_modifiers)
                
        return l
    def matrix_to_csv(self, matrix, filename):
        import csv
        w = csv.writer(open(filename, 'w'), delimiter=',', quotechar='\"')
        for l in matrix:
            w.writerow(l)
        

def main():
    import libsbml
    doc = libsbml.readSBML('/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/BIOMD0000000172.xml')
    inchi = InChI(doc)
    matrix = inchi.get_closest_matches_for_reactions()
    inchi.matrix_to_csv(matrix, '/home/schulzma/Desktop/inchi_matches.csv')
