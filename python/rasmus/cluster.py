"""
    cluster.py
    
    This module contains functions for creating and manipulating 
    clusterings/partitionings.
    
"""

import os
import sys

import algorithms
import blast
import graph
import matrix
import util



##############################################
# comparing clusterings/partitions
#


def item2part(parts):
    """Creates a lookup dict of items to partition id"""
    
    lookup = {}
    for i in xrange(len(parts)):
        for item in parts[i]:
            lookup[item] = i
    return lookup



def confusionMatrix(parts1, parts2):
    """Returns a confusion matrix of two different partitions of the same 
       items"""
    
    confuse = util.Dict(2, 0)
    
    lookup1 = item2part(parts1)
    lookup2 = item2part(parts2)
    
    items1 = util.makeset(util.flatten(parts1, 1))
    items2 = util.makeset(util.flatten(parts2, 1))
    
    sameset = util.intersect(items1, items2)
    diffset = util.nonintersect(items1, items2)
    
    for item in sameset:
        confuse[lookup1[item]][lookup2[item]] += 1
    
    return confuse, list(diffset)


def partCrossing(parts1, parts2):
    """How often does part1 cross the boundries defined by part2"""
    
    lookup2 = util.Dict(1, -1)    
    lookup2.update(item2part(parts2))
    splits = []
    
    for part1 in parts1:
        collect = util.Dict(1, [])
        for item in part1:
            collect[lookup2[item]].append(item)
        keys = collect.keys()
        keys.sort()
        
        splits.append(map(lambda key: collect[key], keys))
    
    return splits


def sameParts(parts1, parts2):
    """Returns partitions that are exactly the same between 
       'parts1' and 'parts2'
       
       items of 'parts1' and 'parts2' should be hashable.
    """
    
    lookup = {}
    parts3 = []
    
    for part in parts1:
        lookup[tuple(util.sort(part))] = 1
    
    for part in parts2:
        if tuple(util.sort(part)) in lookup:
            parts3.append(part)
    
    return parts3
    


def part2graph(parts):
    """Creates a graph (2d dict) that represents the partitioning
    
       note: items in 'parts' should not be integers,
             (usually strings or objects)
    """
    
    # create graph
    i = 0
    vertices = {}
    for part in parts:
        for item in part:
            vertices.setdefault(item, {})[i] = 1
            vertices.setdefault(i, {})[item] = 1
        i += 1
    return vertices


def unionPart(* partsList):
    """Finds the union of several partitions"""
    
    mat = part2graph(util.concat(* partsList))
    parts = graph.connectedComponents(mat.keys(), lambda x: mat[x].keys())
    
    # remove parition ids from partitioning
    parts = map(lambda part: filter(lambda x: type(x) != int, part), parts)
    
    return parts


def parts2partids(parts, labels):
    lookup = util.list2lookup(labels)
    partids = [0] * len(labels)
    
    for i in xrange(len(parts)):
        for item in parts[i]:
            partids[lookup[item]] = i
    
    return partids


def partids2parts(partids, labels):
    nparts = max(partids) + 1
    parts = []
    
    for i in xrange(nparts):
        parts.append([])
    
    for i, label in zip(partids, labels):
        # skip i = -1, means item is not clustered
        if i >= 0:
            parts[i].append(label)
    
    return parts


def filterOne2ones(parts, gene2species):
    def isOne2one(part, gene2species):
        counts = util.histDict(map(gene2species, part))
        return (max(counts.values()) == 1)

    # get one2ones
    ones = [x for x in parts if isOne2one(x, gene2species)]

    # find maximum one2one
    maxsize = max(len(x) for x in ones)
    ones = [util.sort(x) for x in ones if len(x) == maxsize]

    return ones


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


    util.tic("read outgroup hits")    
    outbest = util.Dict(default=[0, 0])

    for blastfile, order in outblastfiles:
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
            
            # create a key for a partition: (side, index)
            if gene1 in lookup1:
                partin = (0, lookup1[gene1])
            else:
                parts1.append([gene1])
                lookup1[gene1] = len(parts1) - 1
                partin = (0, len(parts1) - 1)
            
            if gene1 in lookup2:
                partin = (1, lookup2[gene2])
            else:
                parts2.append([gene2])
                lookup2[gene2] = len(parts2) - 1
                partin = (1, len(parts2) - 1)
            
            val = outbest[partin]
            val[0] += score
            val[1] += 1
            
        util.toc()
    util.toc()
    
    
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
    

    util.tic("determine clusters")
    sets = {}
    for part in hits:
        sets[part] = algorithms.UnionFind([part])
    
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






"""

def partTreeFiles(parts, matchfiles, col1, col2, scorecol, minscore = 500):
    # format: mat[i][j] = [sum, num]
    mat = util.Dict(2, [0, 0])

    lookup = item2part(parts)
    
    # read in matches
    for f in matchfiles:
        print >>sys.stderr, f

        for line in file(f):
            tokens = line.split("\t")
            gene1 = tokens[col1]
            gene2 = tokens[col2]
            score = float(tokens[scorecol])

            if score < minscore:
                continue

            if gene1 in lookup and gene2 in lookup:
                i = lookup[gene1]
                j = lookup[gene2]
                
                if i == j:
                    continue
                
                if i < j:
                    i, j = j, i
                
                sumnum = mat[i][j]
                sumnum[0] += score
                sumnum[1] += 1
            
    return partTree(parts, mat)

def writeMatrix(mat, inuse):
    avgs = []
    for i in mat:
        for j in mat[i]:
            if inuse[i] and inuse[j]:
                sumnum = mat[i][j]
                avg = sumnum[0] / sumnum[1]
                avgs.append([avg, i, j])
    
    avgs.sort(lambda a,b: cmp(b[0], a[0]))
    
    for avg in avgs:
        print "\t".join(map(str, avg))


def get(mat, i, j):
    if i < j:
        i, j = j, i
    if j in mat[i]:
        return mat[i][j]
    else:
        return None    

def findMaxPair(mat, inuse):
    # find max avg-link
    top = [0, -1, -1]
    for i in mat:
        for j in mat[i]:
            if inuse[i] and inuse[j]:
                sumnum = mat[i][j]
                avg = sumnum[0] / sumnum[1]
                if avg > top[0]:
                    top = [avg, i, j]
    if top[1] != -1:
        return top
    else:
        return None


def merge(mat, inuse, i, j):
    # merge clusters i and j
    l = len(inuse)

    inuse.append(1)
    for k in xrange(l):
        if k == i or k == j:
            continue
        
        mat1 = get(mat, i, k)
        mat2 = get(mat, j, k)
        
        if mat1 == None:
            if mat2 == None:
                continue        
            mat1 = [0, 0]
        else:
            if mat2 == None:
                mat2 = [0, 0]
        
       
        
        a, b = l, k
        if a < b:
            a, b = b, a
        
        mat[a][b] = [mat1[0] + mat2[0], mat1[1] + mat2[1]]
    
    inuse[i] = 0
    inuse[j] = 0
    
    return l



def partTree(parts, mat):
    # mark parts that are in use
    inuse = [1] * len(parts)
    roots = {}.fromkeys(range(len(parts)), 1)
    
    tree = algorithms.Tree()
    for i in xrange(len(parts)):
        tree.add(algorithms.TreeNode(i))
    
    while True:
        util.log("number of parts left %d" % len(roots))
        #writeMatrix(mat, inuse)
        #print

        top = findMaxPair(mat, inuse)
    
        if top != None:
            #util.log("top is " + top)
    
            newnode = merge(mat, inuse, top[1], top[2])
            roots[newnode] = 1
            
            tree.add(algorithms.TreeNode(newnode))
            tree.addChild(tree.nodes[newnode], tree.nodes[top[1]])
            tree.addChild(tree.nodes[newnode], tree.nodes[top[2]])
            del roots[top[1]]
            del roots[top[2]]
            
        else:
            break

    # create roots for tree
    tree.makeRoot()
    for root in roots:
        tree.addChild(tree.root, tree.nodes[root])

    return tree

"""

