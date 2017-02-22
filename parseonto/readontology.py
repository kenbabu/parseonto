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
from operator import itemgetter
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



args = parser.parse_args()

# Create the ontology object
obo = orangeonto.OBOOntology(file=args.file)




# Find the root term of an ontology file
#  Change should be determined automatically.
ROOT = obo.root().id





def numOfTermsFile():
    return  len(obo.terms())


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




#  To refactor

idx_term = createIndicesTerm(ROOT)



term_idx = createTermIndices(idx_term)

# print term_idx

# Provide a brief summary of the ontology- show root and number of terms
def summariseOntology():
    print
    print "="*80
    print "Getting summary for the ontology located in {}".format(args.file)
    print
    print "Name of  root: [{}], the id of the root: [{}] ".format(obo.root().name, obo.root().id)
    print "Number of terms: {}".format(numOfTermsFile())
    print "Number of obsolete terms: {}".format(obsoleteTerms(obo).__len__())
    print "Edge types  in the ontology are {}".format(obo.edge_types().__len__())
    # print "Typedefs: {}".format(obo.typedefs())
    # print "Number of obsolete terms: {}".format(len([i for i in obo.terms() if obo.term(i).obsolete=='true']))
    print
    print "="*80




def mca (term1, term2):
    anc_term1= [i.id for i in obo.super_terms(term1)]
    anc_term2= [i.id for i in obo.super_terms(term2)]

    anc=  max([term_idx[i] for i in  list(set(anc_term1).intersection( anc_term2))])

    return idx_term[anc]

def calcSimilarity(term1, term2):
    try:
        return float(-1* (math.log(obo.topology(mca(term1, term2)))))/max([-1*math.log(obo.topology(term1)),
                                                                 -1*math.log(obo.topology(term2))])
    except ValueError:
        # print "Error: Common Ancestors undefined {0}\t{1}\t[NAN]".format(term1, term2 )
        return None
    except TypeError:
        # print "TypeError {0}\t{1}\t[NAN]".format(term1, term2 )
        return None

#     Number of pairs in the ontology
def genPairwise(lsTerms):
    return combinations(lsTerms, 2)

# print mca('DOID:679', 'DOID:9771')

# print calcSimilarity('DOID:679', 'DOID:971')

# Calc similarity ebola and malaria
def testSimCalc(term1, term2):


    MCA = mca(term1, term2)
    print "Topology of MCA ", obo.topology(mca(term1, term2))

    print "Topology of  {0}: {1} ".format(term1, obo.topology(term1))
    print "Topology of  {0}: {1} ".format(term2, obo.topology(term2))

    ICA =  (-1* math.log(obo.topology(term1)))
    print "Information Content of  {0}: {1} ".format(term1, ICA)
    ICB =  (-1* math.log(obo.topology(term2)))

    print "Information Content of  {0}: {1} ".format(term2,ICB)
    ICMCA =  (-1* math.log(obo.topology(MCA)))

    print "Information Content of  the MCA [{0}]: {1} ".format(MCA,ICMCA)

    sim = ICMCA/max([ICA, ICB])
    print "Similarity: [{0},{1}]: {2}".format(term1, term2, sim)

    # print math.log(1.0)/max(obo.topology('DOID:0060405'),obo.topology('DOID:7788'))

# Calculate the time of computation

def calcGlobalSimilarity(term):

    lsofdicts= [{(term, i.id):calcSimilarity(term, i.id) }for i in obo.terms()]

    dictTerms = { k: v for d in lsofdicts for k, v in d.items() }
    return  sorted(dictTerms.iteritems(), key=itemgetter(1), reverse=True)
    # return sorted([calcSimilarity(term, i.id) for i in obo.terms()])

# Similarity calculation for a list of  ontology terms
def calcSimilarityCombinations():
    start = time.clock()
    with open(os.path.join(picklepath, args.similarities), "w") as simFile:
        count =0
        for i in combinations(loadTermIndex().values(), 2):
            # ignore root
            if not (i[0] == ROOT or i[1] == ROOT):
                simFile.writelines("{0}\t{1}\t{2}\n".format(calcSimilarity(i[0],i[1]), i[0], i[1]))
                count +=1
                print i[0],i[1], calcSimilarity(i[0],i[1]), count

    end = time.clock()
    print "Took {0} seconds to calculate similarities for {1} pairs ".format((end-start)/1000, count )




# print obo.root().id

# print obo.topology('DOID:4')
print obo.topology('DOID:4325')
# print obo.topology('DOID:11341')

def obsoleteTerms(onto):
    lsObsoletes =[]
    for term in  onto.terms():
        if term.get_value('is_obsolete') != []:
            lsObsoletes.append(term.id)
    return  lsObsoletes



# print obsoleteTerms(obo).__len__()
summariseOntology()
#
# testSimCalc('DOID:4325', 'DOID:11341')

# print mca('DOID:4325', 'DOID:11341')




# # DOID:9084
# DOID:9092
# DOID:9093
# DOID:9103
# DOID:9105
# DOID:9109
# DOID:9114
# DOID:9115
# DOID:9117
