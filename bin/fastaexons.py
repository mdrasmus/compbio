#!/usr/bin/env python

from rasmus import fasta, genomeio, genomeutil, util
import sys


options = [
#    ["a:", "align=", "align", "AUTO<fasta alignment>"]
    ["g:", "genomes=", "genomes", "AUTO<genome1>,<genome2>,<...>"]
]

#try:
#    param, rest = util.parseArgs(sys.argv, options)
#except:
#    sys.exit(1)


param = {}
param["genomes"] = ["dog,human,mouse,rat"]
rest = ["../mammalian/2005-07-09/out/1.align"]


if True:
    # read exon information
    genomes = {}
    for genome in param["genomes"][-1].split(","):
        util.tic("read %s exons" % genome)
        genomes[genome] = genomeutil.Genome(genome)

        genomeio.readEnsmartExons(genomes[genome], 
                                  genomeio.structfile(genome))
        genomes[genome].autoconf()
        util.toc()

genes = {}
for genome in genomes.values():
    genes.update(genome.genes)


# process each fasta file
for filename in rest:
    util.log("processing '%s'" % filename)
    aln = fasta.readFasta(filename)

    for gene in aln:
        print len(genes[gene].trans)
    



