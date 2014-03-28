"""
    cluster.py
    
    This module contains functions for creating and manipulating 
    clusterings/partitionings.
    
"""

import os
import sys

from rasmus import graph
from rasmus import util

from compbio import blast


# TODO: update function name styles


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
    
    confuse = util.Dict(dim=2, default=0)
    
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
        counts = util.hist_dict(map(gene2species, part))
        return (max(counts.values()) == 1)

    # get one2ones
    ones = [x for x in parts if isOne2one(x, gene2species)]

    # find maximum one2one
    maxsize = max(len(x) for x in ones)
    ones = [util.sort(x) for x in ones if len(x) == maxsize]

    return ones


