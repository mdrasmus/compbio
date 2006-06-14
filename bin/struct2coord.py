#!/usr/bin/env python
from rasmus import util, genomeio, genomeutil
import sys



filename = sys.argv[1]

genome = genomeutil.Genome()
genomeio.readEnsmartGenes(genome, filename)

for gene in genome.genes.values():
    print "\t".join([gene.name, gene.chrom.name, 
                    str(gene.start), str(gene.end), str(gene.direction)])

