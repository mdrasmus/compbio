#!/usr/bin/env summon

# python libs
import os, sys

# rasmus libs
from rasmus import util, env, genomeio, genomeutil, clustalw, fasta
from rasmus.genomeutil import *

# graphics libs
from summon import *
from rasmus.vis import dotplot


options = [
    ["g:", "genomes=", "genomes", "<genome1>,<genome2>,..."],
    ["c:", "comps=", "comps", "<ortholog components file>"],
    ["s:", "synteny=", "synteny", "<synteny file>"],
    ["q", "sequence", "sequence", "", 
        {"help": "automatically read sequences"}],
    ["S:", "speciesmap=", "speciesmap", "<species map>"],
    ["P:", "paths=", "paths", "<data paths>", 
        {"default": ".",
         "single": True}]
]

param = util.parseOptions(sys.argv[1:], options, quit=True)


# globals
# m    - matching
# seqs - sequences
# vis  - visualization
# conf - visualization configuration
#

def main(param):
    env.addPaths(param["paths"])
    env.addEnvPaths("DATAPATH")



    # execute
    util.tic("read")

    # get species map
    if "speciesmap" in param:
        gene2species = readGene2species(* map(env.findFile, param["speciesmap"]))
    else:
        gene2species = genomeutil.gene2species

    # read genomes
    m = Matching()
    genomes = param["genomes"][-1].split(",")
    globals()["m"] = m
    genomeio.readGenomes(m, genomes, gene2species)
    
    
    util.tic("read syntenic orthologs")
    for orth in param["comps"]:
        genomeio.readSimpleMatches(m, env.findFile(orth))
    util.toc()

    util.tic("read synteny blocks")
    for block in param["synteny"]:
        genomeio.readBlockDimensions(m, env.findFile(block), True)
    util.toc()
        
    util.tic("read sequences")    
    seqs = fasta.FastaDict()
    globals()["seqs"] = seqs
    if "sequence" in param:
        for genome in genomes:
            try:
                seqfile = env.findFile("%s.fasta" % genome)
                util.tic("reading '%s'" % seqfile)
                seqs.read(seqfile)
                util.toc()
            except: pass
    util.toc()
    
    util.toc()
    
    
    
    # draw dotplot
    util.tic("drawing")
    set_bgcolor(1,1,1)
    conf = dotplot.initConf()
    vis = dotplot.Dotplot()
    globals()["vis"] = vis
    globals()["conf"] = conf
    vis.setSeqs(seqs)
    add_group(vis.draw(conf, m, genomes))
    home()
    util.toc()








print
dotplot.displayBindings()


def mark(markColor = color(0,0,1, .2)):
    names = []
    while True:
        line = sys.stdin.readline().rstrip()
        if line == "": break
        names.extend(line.rstrip().split())
    vis.markGenes(names, markColor)


def align(* names):
    seqs = map(lambda x: genes[x].protein(), names)
    aln = clustalw.clustalw(seqs)
    clustalw.printAlign(aln)



main(param)
    
