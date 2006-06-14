#!/usr/bin/python

from rasmus import util, graph, fasta, stats
import getopt
import sys


#    ["l:", "filterlen=", "filterlen", "AUTO<minfrac>,<minlen>"],
#    ["O:", "filterout=", "filterout", "AUTO<output file of filtered genes>"]



if len(sys.argv) < 2:
    print "usage: genepart.py [-p partitions] [-m matches]"
    sys.exit(1)


params, rest = util.getopt(sys.argv[1:], "m:p:s:l:O:f:g:")

vertices = util.Dict(2)



# read partition
seeds = {}
parts = util.readDelim(params["-p"][0])
for i in xrange(len(parts)):
    for gene in parts[i]:
        vertices[gene][i] = 1
        vertices[i][gene] = 1
        seeds[gene] = i

if "-g" in params:
    genes = util.makeset(util.readStrings(params["-g"][-1]))
else:
    genes = None

if "-O" in params:
    outfiltered = file(params["-O"][-1], "w")
else:
    outfiltered = file("filtered.genes", "w")


# all added matches form paths leading to a single partition or no partition
# this prevents connecting existing partitions


# pull in best unidirectional matches
keep = {}
for matchfile in params["-m"]:
    print >>sys.stderr, matchfile
    for line in file(matchfile):
        (gene1, gene2, score) = line.split()[:3]
        score = float(score)
        
        # skip self-matches
        if gene1 == gene2:
            continue
        
        # do not form matches between partitions
        if gene1 in seeds and gene2 in seeds:
            continue
            
        # only use genes from acceptable list
        if genes and (gene1 not in genes or gene2 not in genes):
            continue
        
        # only find best matches for genes not in partition        
        if gene1 not in seeds and \
           (gene1 not in keep or score > keep[gene1][1]):
            keep[gene1] = (gene2, score)
        
        if gene2 not in seeds and \
           (gene2 not in keep or score > keep[gene2][1]):
            keep[gene2] = (gene1, score)


# add kept genes to vertices if other gene is in partitions
for gene, match in keep.items():
    if match[0] in seeds:
        vertices[gene][match[0]] = 1
        vertices[match[0]][gene] = 1

# find connected components
def getNeighbors(vertex):
    return vertices[vertex].keys()
comps = graph.connectedComponents(vertices.keys(), getNeighbors)



# filter components

def filterSeqs(seqs, minseqfrac, minseqlen):
    lens = map(len, seqs.values())
    mid = float(stats.mean(lens))
    lens2 = map(lambda x: x/mid, lens)
    keep = map(lambda x: x[0] >= minseqfrac or \
                              x[1] >= minseqlen, zip(lens2, lens))
    ind = util.findeq(True, keep)
    keys = util.sublist(seqs.keys(), ind)

    # record filtered genes
    ind = util.findeq(False, keep)
    rkeys = util.sublist(seqs.keys(), ind)        
    for key in rkeys:
        print >>outfiltered, key
    return keys
    

def filterComps(comps, seqs, minseqfrac, minseqlen):
    comps2 = []
    for comp in comps:
        comps2.append(filterSeqs(util.subdict(seqs, comp), 
                                 minseqfrac, minseqlen))
    return comps2


if "-f" in params:
    seqs = fasta.readFasta(params["-f"][-1])
    minseqfrac, minseqlen = map(float, params["-l"][-1].split(","))
    
    comps = filterComps(comps, seqs, minseqfrac, minseqlen)
    

# output
for comp in comps:
    for gene in comp:
        # don't print partition ids
        if type(gene) == str:
            print gene,
    print

