#!/usr/bin/env python

import sys, os, re
from rasmus import util, phyloutil, env, blast, graph, synphyllib

options = [
    ["i:", "inhits=", "inhits", "<sorted blast file>"],
    ["x:", "crosshits=", "crosshits", "<blast hits crossing in & out>"],
    ["I:", "ingenomes=", "ingenomes", "<genome1>,<genome2>..."],
    ["O:", "outgenomes=", "outgenomes", "<genome1>,<genome2>..."],
    ["s:", "sortblast=", "sortblast", "<unsorted blast file>"],
    ["p:", "part=", "part", "<output partition>", {
        "single": True,
        "req": True}],
    ["P:",  "paths=", "paths", "<data path>"]
]


def findBlastOrders(crossBlast, ingenomes, outgenomes):
    orders = []
    
    for blastfile in crossBlast:
        try:
            groups = re.match(r"(^|.*/)(?P<g1>[^/]+)_(?P<g2>.+)\.blastp$", 
                              blastfile).groupdict()
            
            if groups["g1"] in ingenomes and \
               groups["g2"] in outgenomes:
                orders.append(True)
            elif groups["g1"] in outgenomes and \
               groups["g2"] in ingenomes:
                orders.append(False)
            else:
                raise "Blast file '%s' does not cross the in and out groups" \
                      % blastfile
        except:
            raise "Bad blast filename '%s'" % blastfile
    
    return orders
        
def cluster(inBest, outBest, sortedBlast):
    mat = util.Dict(2, 0)
    
    util.tic("cluster")
    for hit in blast.BlastListReader(sortedBlast):
        gene1 = blast.query(hit)
        gene2 = blast.subject(hit)
        score = blast.bitscore(hit)
        
        if inBest[gene1] >= score and \
           inBest[gene2] >= score:
            mat[gene1][gene2] = 1
            mat[gene2][gene1] = 1
    util.toc()
    
    comps = graph.connectedComponents(mat.keys(), lambda x: mat[x].keys())
    return comps


def main(param):
    ingenomes = param["ingenomes"][-1].split(",")
    outgenomes = param["outgenomes"][-1].split(",")

    crossBlast = param["crosshits"]
    orders = findBlastOrders(crossBlast, ingenomes, outgenomes)
    
    sortedBlast = param["inhits"]


    # cluster by outgroup
    if True:
        inBest, outBest = synphyllib.outgroupCutoff(crossBlast, orders)
        
        util.log("%d hits from IN to OUT" % len(inBest))
        
        #comps, scores = synphyllib.clusterIngroup(inBest, parts, sortedBlast)
        
        comps = cluster(inBest, outBest, sortedBlast)

        util.writeDelim(param["part"], comps)



main(util.parseOptions(sys.argv, options, quit=True))
