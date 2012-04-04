#!/usr/bin/python -i

# python libs
import os, sys

# rasmus libs
from rasmus import util, env
from rasmus.bio import genomeio, genomeutil, clustalw, fasta, gff

# graphics libs
from summon.core import *
from rasmus.vis import dotplot


options = [
    ["G:", "genomes=", "genomes", "<genome1>,<genome2>",
        {"single": True}],
    ["g:", "gff=", "gff", "<gff file>"],
    ["c:", "comps=", "comps", "<ortholog components file>"],
    ["s:", "synteny=", "synteny", "<synteny file>",
        {"default": []}],
    ["q", "sequence", "sequence", "", 
        {"help": "automatically read sequences"}],
    ["S:", "speciesmap=", "speciesmap", "<species map>"],
    ["P:", "paths=", "paths", "<data paths>", 
        {"default": ".",
         "single": True}]
]

conf = util.parseOptions(sys.argv, options, quit=True)


# globals
# m    - matching
# seqs - sequences
# vis  - visualization
# conf - visualization configuration
#


env.addPaths(conf["paths"])
env.addEnvPaths("DATAPATH")



# execute
util.tic("read")

# get species map
if "speciesmap" in conf:
    gene2species = genomeutil.readGene2species(* map(env.findFile, conf["speciesmap"]))
else:
    gene2species = genomeutil.gene2species

# read genomes
m = genomeutil.Matching()
genomes = conf["genomes"].split(",")
for gffFile in conf["gff"]:
    util.tic("read '%s'" % gffFile)
    regions = gff.readGff(gffFile, format=gff.GFF3)
    m.addRegions(regions, gene2species)
m.autoconf(genomes)



util.tic("read syntenic orthologs")
for orth in conf["comps"]:
    genomeio.readSimpleMatches(m, env.findFile(orth))
util.toc()

util.tic("read synteny blocks")
for block in conf["synteny"]:
    genomeio.readBlockDimensions(m, env.findFile(block), True)
util.toc()

util.tic("read sequences")    
seqs = fasta.FastaDict()
if "sequence" in conf:
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
conf = dotplot.initConf()
vis = dotplot.Dotplot2(conf, m, genomes)
vis.setSeqs(seqs)
vis.show()
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

