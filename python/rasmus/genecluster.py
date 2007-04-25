"""
    genecluster.py
    
    This module contains functions for creating gene clusters based on sequence
    
"""

import os
import sys

from rasmus import algorithms
from rasmus import blast
from rasmus import graph
from rasmus import matrix
from rasmus import util
from rasmus.cluster import *



####################################
# Hierarchical clustering
#

def mergeBuh(conf, parts1, parts2, blastfiles):
    """Merge by Best Unidirectional Hits"""
    
    lookup1 = item2part(parts1)
    lookup2 = item2part(parts2)
    
    
    best = util.Dict(1, (0, None))

    util.tic("read hits")
    for blastfile, order in blastfiles:
        util.tic("determine best hits '%s'" % os.path.basename(blastfile))
        for hit in blast.BlastReader(blastfile):
            if order:
                gene1 = blast.query(hit)
                gene2 = blast.subject(hit)
            else:
                gene2 = blast.query(hit)
                gene1 = blast.subject(hit)            
            score = blast.bitscore(hit)
            
            if blast.evalue(hit) > conf["signif"]:
                continue
            
            
            if gene1 in lookup1:
                part1 = (0, lookup1[gene1])
            else:
                parts1.append([gene1])
                lookup1[gene1] = len(parts1) - 1
                part1 = (0, len(parts1) - 1)
            
            if gene2 in lookup2:
                part2 = (1, lookup2[gene2])
            else:
                parts2.append([gene2])
                lookup2[gene2] = len(parts2) - 1
                part2 = (1, len(parts2) - 1)
            
            
            if score > best[part1][0]:
                best[part1] = (score, part2)
            if score > best[part2][0]:
                best[part2] = (score, part1)
        util.toc()
        
        
        
    util.toc()

    util.tic("determine clusters")
    sets = {}
    for gene in best:
        sets[gene] = algorithms.UnionFind([gene])
    
    for blastfile, order in blastfiles:
        util.tic("read hits '%s'" % os.path.basename(blastfile))
        for hit in blast.BlastReader(blastfile):
            if order:
                gene1 = blast.query(hit)
                gene2 = blast.subject(hit)
            else:
                gene2 = blast.query(hit)
                gene1 = blast.subject(hit)            
            score = blast.bitscore(hit)
            
            if blast.evalue(hit) > conf["signif"]:
                continue
            
            
            part1 = (0, lookup1[gene1])
            part2 = (1, lookup2[gene2])        

            if score >= best[part1][0] * conf["relcutoff"]:
                sets[part1].union(sets[part2])
            if score >= best[part2][0] * conf["relcutoff"]:
                sets[part2].union(sets[part1])
        util.toc()
    
    
    sets = util.unique([x.root() for x in sets.values()])
    
    parts = []
    joining = (parts1, parts2)
    for set in sets:
        parts.append([])
        for i, row in set.members():
            parts[-1].extend(joining[i][row])
    util.toc()

    return parts



def mergeAvg(conf, parts1, parts2, blastfiles, outblastfiles):
    lookup1 = item2part(parts1)
    lookup2 = item2part(parts2)
    
    # value is [sum, total]
    hits = util.Dict(dim=2, default = [0, 0])
    
    if "accept" in conf:
        accept = conf["accept"]
    else:
        accept = False
    
    
    util.tic("read hits")
    for blastfile, order in blastfiles:
        util.tic("determine best hits '%s'" % os.path.basename(blastfile))
        for hit in blast.BlastReader(blastfile):
            if order:
                gene1 = blast.query(hit)
                gene2 = blast.subject(hit)
            else:
                gene2 = blast.query(hit)
                gene1 = blast.subject(hit)            
            score = blast.bitscore(hit)
            
            if blast.evalue(hit) > conf["signif"]:
                continue
            
            if accept and \
               (gene1 not in accept or
                gene2 not in accept):
                 continue
            
            # create a key for a partition: (side, index)
            if gene1 in lookup1:
                part1 = (0, lookup1[gene1])
            else:
                parts1.append([gene1])
                lookup1[gene1] = len(parts1) - 1
                part1 = (0, len(parts1) - 1)
            
            if gene2 in lookup2:
                part2 = (1, lookup2[gene2])
            else:
                parts2.append([gene2])
                lookup2[gene2] = len(parts2) - 1
                part2 = (1, len(parts2) - 1)
            
            val = hits[part1][part2]
            val[0] += score
            val[1] += 1
            hits[part2][part1] = val
            
        util.toc()
    util.toc()
    
    
    util.tic("read outgroup hits")    
    outbest = util.Dict(default=[0, 0])
    for blastfile, order in outblastfiles:
        util.tic("determine best hits '%s'" % os.path.basename(blastfile))
        for hit in blast.BlastReader(blastfile):
            if order:
                genein  = blast.query(hit)
                geneout = blast.subject(hit)
            else:
                geneout = blast.query(hit)
                genein = blast.subject(hit)            
            score = blast.bitscore(hit)
            
            if blast.evalue(hit) > conf["signif"]:
                continue
            
            # create a key for a partition: (side, index)
            if genein in lookup1:
                partin = (0, lookup1[genein])
            elif gene1 in lookup2:
                partin = (1, lookup2[genein])
            else:
                continue
            
            val = outbest[partin]
            val[0] += score
            val[1] += 1
            
        util.toc()
    util.toc()
    
    assert len(parts1) == len(unionPart(parts1))
    assert len(parts2) == len(unionPart(parts2))
    

    util.tic("determine clusters")
    sets = {}
    for i in xrange(len(parts1)):
        sets[(0, i)] = algorithms.UnionFind([(0, i)])
    for i in xrange(len(parts2)):
        sets[(1, i)] = algorithms.UnionFind([(1, i)])

    
    # merge top avg hits
    for part1 in hits:
        o1 = outbest[part1]
        outavg1 = float(o1[0]) / max(o1[1], 1)
        
        top = 0
        toppart = None
        
        for part2, (tot, num) in hits[part1].iteritems():
            avg = float(tot) / num
            o2 = outbest[part2]
            outavg2 = float(o2[0]) / max(o2[1], 1)
            
            if avg > outavg1 and avg > outavg2 and avg > top:
                top = avg
                toppart = part2
                
        if toppart:
            sets[part1].union(sets[toppart])
    
    sets = util.unique([x.root() for x in sets.values()])
    
    # create partition of genes
    parts = []
    joining = (parts1, parts2)
    for set in sets:
        parts.append([])
        for i, row in set:
            parts[-1].extend(joining[i][row])
    util.toc()
    
    assert len(parts) == len(unionPart(parts))
    
    return parts


"""
# merge top avg hits
    for part1 in hits:
        top = 0
        toppart = None
        for part2, (tot, num) in hits[part1].iteritems():
            avg = float(tot) / num
            
            if avg > top:
                toppart = part2
                top = avg
        
        if toppart:
            sets[part1].union(sets[toppart])
"""


def mergeTree(conf, stree, gene2species, blastFileLookup):
    util.tic("cluster all genes")

    for node in stree.nodes.values():
        node.parts = []
    
    
    # walk up tree (post-order)
    def walk(node):
        for child in node.children:
            walk(child)

        if not node.isLeaf():
            blastfiles = []
            leaves1 = node.children[0].leafNames()
            leaves2 = node.children[1].leafNames()
            
            # determine sibling blast files
            for leaf1 in leaves1:
                for leaf2 in leaves2:
                    if leaf1 in blastFileLookup and \
                       leaf2 in blastFileLookup[leaf1]:
                        blastfiles.append(blastFileLookup[leaf1][leaf2])
            
            # determine outgroup blast files (all other files, potentially)
            # go up one level, blastfiles for leaves, and subtract sibling files
            outblastfiles = []
            if node.parent:
                inleaves = leaves1 + leaves2
                outleaves = set(node.parent.leafNames()) - set(inleaves)
                
                for leaf1 in inleaves:
                    for leaf2 in outleaves:
                        if leaf1 in blastFileLookup and \
                           leaf2 in blastFileLookup[leaf1]:
                            outblastfiles.append(blastFileLookup[leaf1][leaf2])
                        
            util.tic("merging")
            util.logger("leaves1: ", leaves1)
            util.logger("leaves2: ", leaves2)
            
            if "merge" in conf and \
               conf["merge"] == "avg":
                node.parts = mergeAvg(conf,
                                      node.children[0].parts,
                                      node.children[1].parts,
                                      blastfiles,
                                      outblastfiles)
            else:
                node.parts = mergeBuh(conf,
                                      node.children[0].parts,
                                      node.children[1].parts,
                                      blastfiles)
            
            if "output" in conf and len(node.parts) > 0:
                util.writeDelim(conf["output"] + 
                                str(node.name) + 
                                ".part", node.parts)
            
            util.logger("number of parts: ", len(node.parts))
            if len(node.parts) > 0:
                util.logger("largest part:", max(map(len, node.parts)))
            
            util.toc()

    walk(stree.root)

    util.toc()

    return stree


def makeBlastFileLookup(blastfiles):
    lookup = util.Dict(2)
    
    for f in blastfiles:
        m = util.match("(|.*/)(?P<genome1>[^/_]+)_(?P<genome2>[^/\.]+)\.[\/]*", f)
        lookup[m["genome1"]][m["genome2"]] = (f, True)
        lookup[m["genome2"]][m["genome1"]] = (f, False)

    return lookup


