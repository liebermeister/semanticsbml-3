#!/usr/bin/env python

import copy

def hash_elem(elem1, elem2):
    return str(elem1)+"|"+str(elem2)

class Cluster:
    def __init__(self,element_list,clustering):
        self.element_list = copy.deepcopy(element_list)
        self.linkage = clustering.linkage
        self.similarity_hash = clustering.similarity_hash
        self._dead = False
    def die(self):
        self._dead = True
    def dead(self):
        return self._dead
    def sim(self, other):
        similarities = []
        for elem1 in self.element_list:
            for elem2 in other.element_list:
                similarities.append(self.similarity_hash[hash_elem(elem1, elem2)])
        if self.linkage == "single":
            return max(similarities)
        elif self.linkage == "average":
            return sum(similarities)/float(len(similarities))
        elif self.linkage == "complete":
            return min(similarities)

class Clustering:
    def __init__(self, linkage="single"):
        self.linkage = linkage
        self.elements = []
        self.similarity_hash = {}
    def load_from_similarity_list(self, similarity_list):
        """List of [elem1, elem2, sim] as an input."""
        for [elem1, elem2, sim] in similarity_list:
            self.similarity_hash[hash_elem(elem1, elem2)] = sim
            self.similarity_hash[hash_elem(elem2, elem1)] = sim
            for elem in [elem1, elem2]:
                if elem not in self.elements:
                    self.elements.append(elem)
    def merge(self, cluster1, cluster2):
        new_cluster = Cluster(cluster1.element_list+cluster2.element_list, self)
        cluster1.die()
        cluster2.die()
        return new_cluster
    def cluster(self):
        self.cluster_history = []
        current_clusters = [Cluster([elem],self) for elem in self.elements]
        for i in range(len(current_clusters)-1):
            self.cluster_history.append(copy.deepcopy(current_clusters))
            best_indices = [None, None]
            best_sim = -1.
            for index1 in range(len(current_clusters)):
                for index2 in range(index1+1,len(current_clusters)):
                    this_sim = current_clusters[index1].sim(current_clusters[index2])
                    if this_sim > best_sim:
                        best_sim = this_sim
                        best_indices = [index1, index2]
            [index1, index2] = best_indices
            merged_cluster = self.merge(current_clusters[index1], current_clusters[index2])
            current_clusters = [c for c in current_clusters if not c.dead()]+[merged_cluster]
        self.cluster_history.append(copy.deepcopy(current_clusters))
    def to_ordered_element_list(self):
        return self.cluster_history[-1][0].element_list
    def to_ordered_matrix(self):
        element_list = self.to_ordered_element_list()
        matrix = [[""]+element_list]
        for index1 in range(len(element_list)):
            elem1 = element_list[index1]
            line = []
            for index2 in range(len(element_list)):
                if index1 == index2:
                    line.append(1.)
                else:
                    elem2 = element_list[index2]
                    line.append(self.similarity_hash[hash_elem(elem1, elem2)])
            matrix.append([element_list[index1]]+line)
        return matrix
    def to_ordered_dendogram_matrix(self):
        element_list = self.to_ordered_element_list()
        matrix = []
        for i in range(1,len(self.cluster_history)):
            line = []
            last_cluster = self.cluster_history[i][-1]
            new_element_list = last_cluster.element_list
            for elem in element_list:
                if elem in new_element_list:
                    line.append(1)
                else:
                    line.append(0)
            matrix.append(line)
        return matrix


if __name__ == "__main__":
    import random
    size = 50
    sim_list = []
    for i in range(size):
        for j in range(i+1,size):
            sim_list.append( [i, j, random.random()] )
    c=Clustering()
    c.load_from_similarity_list(sim_list)
    c.cluster()
    print(c.to_ordered_matrix())
    print(c.to_ordered_dendogram_matrix())
