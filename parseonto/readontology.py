__author__ = 'kenneth'



import  orangeonto
import  cPickle as pickle
import  os
import  math
# from orangeonto import OBOOntology, OBOObject, OBOParser, BUILTIN_OBO_OBJECTS

import  argparse

parser = argparse.ArgumentParser()
parser.add_argument("file",  help="Filepath to the obo file to parse ")
parser.add_argument("--outfile","-o",  help="Filepath where the output is written ")
parser.add_argument("--term_index", "-t", help="File path to the pickle objects with term to index mappings")
parser.add_argument("--index_term", "-i", help="File path to the pickle objects with index to term mappings")



args = parser.parse_args()

obo = orangeonto.OBOOntology(file=args.file)

picklepath= os.path.expanduser('~/Desktop/HDO_data')

filepath = os.path.join(picklepath, 'TermIndices.cpk')
filepath_indices = os.path.join(picklepath, 'IndicesTerms.cpk')



# Find the root term of an ontology file
ROOT = 'DOID:4'
# def isRoot(onto):
#     return  [[i.id for i in onto.super_terms(j)] for j in onto.terms()  if onto.super_terms(j) ]
#     # return  [i for i in onto.terms() if onto(i).super_terms == [] ]





def statisticsHDO():
    print "The Statistics of the ontology "
    print "Number of Terms = {}".format(len(obo.terms()))


# Create a dictionary of indices -> terms
def createIndicesTerm(root):
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


def mca (term1, term2):
    anc_term1= [i.id for i in obo.super_terms(term1)]
    anc_term2= [i.id for i in obo.super_terms(term2)]

    anc=  max([term_idx[i] for i in  list(set(anc_term1).intersection( anc_term2))])

    return idx_term[anc]

def calcSimilarity(term1, term2):
    return float(math.log(obo.topology(mca(term1, term2))))/max([math.log(obo.topology(term1)), math.log(obo.topology(term2))])

print mca('DOID:679', 'DOID:9771')

print calcSimilarity('DOID:679', 'DOID:971')

