#!/usr/bin/env python

import os, sys
from rasmus import util, env, genomeutil, genomeio 


options = [
    ["g:", "genomes=", "genomes", "<genome1>,<genome2>,..."],
    ["c:", "orthcomp=", "orthcomp", "<ortholog components>"],
    ["s:", "synteny=", "synteny", "<synteny>"], 
    ["r:", "refgenome=", "refgenome", "<reference genome>",
        {"single": True}],
    ["S:", "smap=", "smap", "<gene2species map file>"],
    ["P:", "paths=", "paths", "<data path1>:<data path2>:...", 
        {"default": "."}]
]



# check args
param = util.parseOptions(sys.argv, options, quit=True)


def main(param):
    # setup paths
    env.addEnvPaths("DATAPATH")
    env.addPaths(param["paths"])

    # read species map
    if "smap" in param:
        gene2species = genomeutil.readGene2species(* map(env.findFile, param["smap"]))
    else:
        gene2species = genomeutil.gene2species

    # read gene components
    comps = util.readDelim(env.findFile(param["orthcomp"][-1]))
    genomes = param["genomes"][-1].split(",")
    
    # filter components such that they only refer to the given genomes
    comps2 = []
    for comp in comps:
        comp2 = filter(lambda x: gene2species(x) in genomes, comp)
        if len(comp2) > 0:
            comps2.append(comp2)
    comps = comps2
    
    matching = readMatching(param, genomes, comps, gene2species)
        
    
    calcMatching(param, genomes, matching, gene2species)
    calcCoverage(param, genomes, matching, gene2species)


def readMatching(param, genomes, comps, gene2species):
    util.tic("reading genomes")
    matching = genomeutil.Matching()
    for genome in genomes:
        util.tic("reading %s" % genome)
        matching.readGenomes(env.findFile("%s.coord" % genome), 
                             gene2species)
        util.toc()
    matching.autoconf()
    util.toc()
    
    matching.setGeneComponents(comps)
    for f in param["synteny"]:
        genomeio.readBlockDimensions(matching, env.findFile(f))
    
    return matching
    

def calcMatching(param, genomes, matching, gene2species):
    # print stats 
    stat = [["genome", "#chroms", "#genes", "#matched", "percent", "size"]]
    
    ngenes = util.mapdict(util.groupby(lambda x: gene2species(x.name), 
                                       util.flatten(matching.comps)),
                          valfunc=len)
    
    # display stats
    for genomename in genomes:
        genome = matching.genomes[genomename]
        
        stat.append([
            genome.name,
            len(genome.chroms), 
            len(genome.genes),
            ngenes[genomename], 
            100 * ngenes[genomename] / float(len(genome.genes)),
            genome.size
            ])
    
    print
    print "Matching stats"
    util.printcols(stat, spacing=2)


def genomeCompositions(names, genomes):
    counts = [0] * len(genomes)
    lookup = util.list2lookup(genomes)
    
    for name in names:
        counts[lookup[name]] += 1
    
    return tuple(counts)


def calcCoverage(param, genomes, matching, gene2species):
    cov = util.Dict(2, 0)
    nblocks = {}
    alignHist = util.Dict(1, 0)
    alignBasesHist = util.Dict(1, 0)
    alignDone = {}
    
    
    for genome1 in genomes:
        util.tic("calc coverage: %s" % genome1)
        multiblocks = genomeutil.makeGenomeMultiBlocks({}, matching, 
                                      matching.genomes[genome1])
        nblocks[genome1] = len(multiblocks)
        
        # calc coverage
        for block in multiblocks:
            seqlen = block.segments[0].end - block.segments[0].start + 1
            
            cov2 = {}
            for seg in block.segments[1:]:
                cov2[seg.genome.name] = 1
            
            for genome2 in cov2:
                cov[genome1][genome2] += seqlen
            
            key = genomeCompositions(map(lambda seg: seg.genome.name, 
                                         block.segments), genomes)
            if key not in alignDone:
                alignHist[key] += seqlen
                alignBasesHist[key] += sum(map(
                    lambda seg: seg.end - seg.start + 1, block.segments))
        
        # make keys in alignHist as done
        for key in alignHist:
            alignDone[key] = 1
            
        util.toc()
    
    # print stats 
    print
    stat = [["genome", "#blocks"]]
    for genome in genomes:
        stat.append([genome, nblocks[genome]])
    util.printcols(stat, spacing=2)
    
    
    
    print
    print "Most common alignments"
    stat = [genomes + ["length", "%len", "%bases"]]
    items = alignHist.items()
    items.sort(lambda a,b: cmp(b[1], a[1]))
    total = float(sum(util.cget(items, 1)))
    allbases = float(sum(map(lambda x: x.size, matching.genomes.values())))
    
    # calculate the coverage of +,+,+,+,+
    atleastone = 0
    atleastoneBases = 0
    for key, val in items:
        if min(key) > 0:
            atleastone += val
    for key, val in alignBasesHist.iteritems():
        if min(key) > 0:
            atleastoneBases += val
    stat.append(["+"]*len(genomes) + 
                [atleastone, 
                 100 * atleastone/total,
                 100 * atleastoneBases/allbases])
    
    for key, val in items:
        counts = []
        for count in key:
            if count == 0:
                counts.append("-")
            else:
                counts.append(str(count))
        stat.append(counts + 
                    [val, 
                     100 * val/total,
                     100 * alignBasesHist[key]/allbases])
    util.printcols(stat, spacing=2)
    
    
    # print stats 
    stat = [[" "] + genomes]
    for genome1 in genomes:
        stat.append([genome1])
        for genome2 in genomes:
            stat[-1].append("%5.2f" % (
                            100 * cov[genome1][genome2] / 
                            float(matching.genomes[genome1].size)))
        
    print
    print "Coverage stats"
    util.printcols(stat, spacing=2)


main(param)
