#!/usr/bin/env python

from rasmus import util, genomeio, descriptions, genomeutil, fasta, env, cluster
import sys



options = [
    ["c:", "orthcomp=", "orthcomp", "<ortholog components file>"],
#    ["r:", "refgenome=", "refgenome", "<reference genome>"],
#    ["s:", "splitgenome=", "splitgenome", "<split genome>"],
    ["d:", "delete=", "delete", "<acceptable deletion size>"],
    ["g:", "genomes=", "genomes", "<genome>"],
    ["S:", "smap=", "smap", "<species map>"],
    ["o:", "output=", "output", "<output file>"],
    ["P:", "paths=", "paths", "<data paths>",
        {"default": "."}]
    ]

# parse options
param = util.parseOptions(sys.argv, options, quit=True)




def countSplits(genes, parts, refgenome, splitgenome, gene2species, diff=1000):
    """Returns number of splits per ortholog group"""
                                     
    splits = {}
    lookup = cluster.item2part(parts)
    
    refgenes = filter(lambda x: gene2species(x) == refgenome, lookup.keys())
    
    # process each reference gene
    for refgene in refgenes:
        part = parts[lookup[refgene]]
        genomes = map(gene2species, part)
        
        if splitgenome in genomes:
            splits[refgene] = numSplit(genes, part, refgene, 
                                       splitgenome, gene2species)
    return splits


def numSplit(genes, part, refgene, splitgenome, gene2species, diff=1000):
    """Returns number of splits for a partition"""

    genomes = map(gene2species, part)
    assert splitgenome in genomes, "missing genomes"
    
    # get the reference gene and its length
    ref = genes[refgene]
    reflen = ref.end - ref.start + 1
    
    # splits can only exist if multiple genes from splitgenome are present
    if util.counteq(splitgenome, genomes) > 1:
        splitgenes = filter(lambda x: gene2species(x) == splitgenome, part)
        
        groups = util.groupby(lambda gene: genes[gene].chrom, splitgenes)
        
        splits = [1]
        
        for chrom, splitgenes2 in groups.items():
            start = min([genes[x].start for x in splitgenes2])
            end = max([genes[x].end for x in splitgenes2])
            
            # are orthologs closer togther than length of reference gene?
            if (reflen - diff < end - start + 1 < reflen + diff):
                # record how many splits there are
                splits.append(len(splitgenes2))
        
        # report maximum spliting for each gene
        return max(splits)
    else:
        # there are no splits
        return 1



def main(param):
    # setup paths
    env.addEnvPaths("DATAPATH")
    env.addPaths(param["paths"])

    # read species map
    gene2species = genomeutil.readGene2species(* map(env.findFile, param["smap"]))

    # read orthologs
    parts = []
    for f in param["orthcomp"]:
        parts.extend(util.readDelim(f))
    
    genomes = param["genomes"][-1].split(",")
    
    # read genomes
    matching = genomeutil.Matching()
    genomeio.readGenomes(matching, param["genomes"][-1].split(","), gene2species)
    genes = matching.getGenes()
    
    if "output" in param:
        out = util.openStream(param["output"][-1], "w")
    
    output = []
    
    output.append(["refgenome", "splitgenome", "splits", "no splits", "excluded"])
    
    for refgenome in genomes:
        for splitgenome in genomes:
            if refgenome == splitgenome:
                continue
            
            splits = countSplits(genes, parts, refgenome, 
                                 splitgenome, gene2species)
            splits2 = filter(util.neqfunc(False), splits.values())
                                 
            output.append([refgenome, splitgenome, 
                    util.countge(2, splits2),
                    util.countlt(2, splits2),
                    util.counteq(False, splits.values())])
            
            # output gene names of splits
            if "output" in param:
                for gene, split in splits.iteritems():
                    if split > 1:
                        print >>out, "%s\t%s\t%d" % (gene, splitgenome, split)
    
    print
    print
    util.printcols(output, spacing=2)

main(param)
