"""
    cluster.py
    
    This module contains functions for creating and manipulating 
    clusterings/partitionings.
    
"""

import os
import sys

from rasmus import algorithms
from rasmus import graph
from rasmus import matrix
from rasmus import util

from rasmus.bio import blast


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
    
    items1 = set(util.flatten(parts1, 1))
    items2 = set(util.flatten(parts2, 1))
    
    sameset = items1 & items2
    diffset = items1.symmetric_difference(items2)
    
    for item in sameset:
        confuse[lookup1[item]][lookup2[item]] += 1
    
    return confuse, list(diffset)


def partLookup(parts1, parts2):
    """For each part in part1, which parts in parts2 share the same items"""
    
    lookup2 = util.Dict(default=-1)    
    lookup2.update(item2part(parts2))
    splits = []
    
    for part1 in parts1:
        hits = set()
        for item in part1:
            hits.add(lookup2[item])
        splits.append(sorted(list(hits)))
    
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

