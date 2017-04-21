__author__ = 'kenneth'



import  orangeonto
import  cPickle as pickle
import  os
import  math
# from orangeonto import OBOOntology, OBOObject, OBOParser, BUILTIN_OBO_OBJECTS

import  argparse
from itertools import combinations
import time
import math
import  csv
from operator import itemgetter
import sys

import pandas as pd
import matplotlib.pyplot as plt

# Path to pickle objects for traversing the ontology
picklepath= os.path.expanduser('~/Desktop/HDO_data')

filepath = os.path.join(picklepath, 'TermIndices.cpk')
filepath_indices = os.path.join(picklepath, 'IndicesTerms.cpk')




parser = argparse.ArgumentParser()
parser.add_argument("file",  help="Filepath to the obo file to parse ")
parser.add_argument("--outfile","-o",  help="Filepath where the output is written ")
parser.add_argument("--term_index", "-t", help="File path to the pickle objects with term to index mappings")
parser.add_argument("--index_term", "-i", help="File path to the pickle objects with index to term mappings")
parser.add_argument("--similarities", "-s", help="File path to the text file where similarities are written")
parser.add_argument("--data_dir", "-d", help="Data directory where the data files are located")



args = parser.parse_args()

# Create the ontology object
obo = orangeonto.OBOOntology(file=args.file)




# Find the root term of an ontology file
#  Change should be determined automatically.
ROOT = obo.root()


class ParseOntology(object):
    def __init__(self, obofile):
        self.obofile = obofile
    def num_terms(self):
        for t in self.obofile.terms:
            print t.id
        # return self.obofile.terms
    def summarise_ontology(self):
        print
        print "="*80
        print "Getting summary for the ontology located in {}".format(args.file)
        print
        print "Name of  root: [{}], the id of the root: [{}] ".format(obo.root().name, obo.root().id)
        print "Number of terms: {}".format(obo.terms().__len__())
        print "Number of current terms: {}".format(self.current_terms().__len__())
        print "Number of obsolete terms: {}".format(self.obsolete_terms().__len__())
        print "Edge types  in the ontology are {}".format(obo.edge_types().__len__())
        print
        print "="*80
    def create_indices_term(self):
        # create a dictionary of  {0:'DOID:4', n:'DOID:xyz'}
        dict_indices = {}
        count=0
        dict_indices[count]=ROOT.id
        for t in obo.traverse_bf(ROOT):
            count += 1
            dict_indices[count]=t.id
        return dict_indices
    def create_term_indices(self):
        indices=self.create_indices_term()
        dict_terms={}
        for k, v in indices.iteritems():
            dict_terms[v]=k
        return dict_terms
    def obsolete_terms(self):
        lsObsoletes =[]
        for term in  obo.terms():
            if term.get_value('is_obsolete') != []:
                lsObsoletes.append(term.id)
        return  lsObsoletes
    def current_terms(self):
        lsCurrent =[]
        for term in  obo.terms():
            if term.get_value('is_obsolete') == []:
                lsCurrent.append(term.id)
        return  lsCurrent
    def calculate_topology(self):
        lsTerms=self.current_terms()
        dictTopology={}
        for t in lsTerms:
            dictTopology[t]=obo.topology(t)
        return dictTopology

    def lowest_common_ancestor(self, term1, term2):
        if (term1 == ROOT.id) or (term2 == ROOT.id):
            return ROOT.id
        else:
            try:
                anc_term1= [i.id for i in obo.super_terms(term1)]
                anc_term2= [i.id for i in obo.super_terms(term2)]
                lca = max([term_index[i] for i in  list(set(anc_term1).intersection(set(anc_term2)))])

                return index_term[lca]
            except ValueError as e:
                print e.message
    def convert_ids_to_indices(self, ls):
        indices= self.create_term_indices()
        try:
            return [indices[term] for term in ls]
        except ValueError as e:
            print e.message

    def calculate_similarity(self, term1, term2):
        mca = self.lowest_common_ancestor(term1, term2)
        if mca == ROOT.id:
            return  0.0
        # topology_a , topology_b= obo.topology(term1), obo.topology(term2)
        topology_a , topology_b= dict_topologies[term1], dict_topologies[term2]
        try:
            similarity = (-math.log(dict_topologies[mca]))/(-math.log(max([topology_a, topology_b])))
            return  similarity
        except ZeroDivisionError:
            return None

# Create a dictionary of indices -> terms
def createIndicesTerm(root):
    # create a dictionary of  {0:'DOID:4', n:'DOID:xyz'}
    dict_bf = {}
    count=0
    dict_bf[count]=root
    for t in obo.traverse_bf(root):
        count += 1
        dict_bf[count]=t.id
    return dict_bf

# Create a dictionary of terms -> indices
def createTermIndices(dict_terms):
    dict_bf = {}
    for k, v in dict_terms.iteritems():
        dict_bf[v] = k
    return dict_bf

def writeIndicesPickle():
    with open(filepath_indices, 'wb') as handle:
        pickle.dump(createIndicesTerm(ROOT), handle, protocol=pickle.HIGHEST_PROTOCOL)


def writeTermsPickle():
    with open(filepath, 'wb') as handle:
        pickle.dump(createTermIndices(ROOT), handle, protocol=pickle.HIGHEST_PROTOCOL)

def loadTermIndex(fpath=args.term_index):
    with open(fpath, 'rb') as handle:
        termdict=pickle.load(handle)
    return  termdict

def loadIndexTerm(fpath=args.index_term):
    with open(fpath, 'rb') as handle:
        termdict=pickle.load(handle)
    return  termdict








parse_onto = ParseOntology(obo)

term_index = parse_onto.create_term_indices()
index_term = parse_onto.create_indices_term()
# topologies = parse_onto.calculate_topology()





lsTen= parse_onto.current_terms()[:10]


def pickle_topologies(topos):
    with open('Topologies.cpk', 'wb') as handle:
        pickle.dump(topos, handle, protocol=pickle.HIGHEST_PROTOCOL)

def loadTopologies(fpath='Topologies.cpk'):
    with open(fpath, 'rb') as handle:
        dictTopology=pickle.load(handle)
    return  dictTopology


dict_topologies = loadTopologies()
print parse_onto.lowest_common_ancestor('DOID:11714', 'DOID:9352')
lsCurrent= parse_onto.current_terms()

pairwise = combinations(lsCurrent, 2)

start=time.time()
with  open("disease_similarities.csv", "wb") as csvfile:
    simwriter = csv.writer(csvfile, delimiter='\t')
    columns =["disease1", "disease2", "similarity"]
    simwriter.writerow(columns)
    for t in pairwise:
        sim=parse_onto.calculate_similarity(t[0], t[1])
        # mca=parse_onto.lowest_common_ancestor(t[0], t[1])
        simwriter.writerow([t[0], t[1], sim])

        print "{0}\t{1}".format(t, sim)
end = time.time()

elapsed_time = (end - start)
print
print "="*30
print "Took {0} seconds for 100 diseases ".format(elapsed_time)



# print len(topologies)
# pickle_topologies(topologies)


#
# for k, v in dict_topologies.iteritems():
#     print k, v

# print parse_onto.create_term_indices()['DOID:7']
# print parse_onto.create_indices_term()[0]c