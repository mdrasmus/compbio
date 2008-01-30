#!/usr/bin/env python

import os, sys

from rasmus import env, util
from rasmus.bio import fasta, genomeutil, alignlib


options = [
    ["u", "ungap", "ungap", "",
        {"single": True}],
    ["c", "noconserved", "noconserved", "",
        {"single": True}],
    ["S:", "smap=", "smap", "<gene2species mapping>"],
    ["P:", "path=", "paths", "<data paths>",
        {"default": []}]
]

conf = util.parseOptions(sys.argv, options, 
                         resthelp="<alignments> ...", quit=True)

# read options
env.addEnvPaths("DATAPATH")
env.addPaths(* conf["paths"])
gene2species = genomeutil.readGene2species(* map(env.findFile, conf["smap"]))


if len(conf["REST"]) > 0:
    alns = map(fasta.readFasta, conf["REST"])
else:
    alns = []
    for line in sys.stdin:
        alns.append(fasta.readFasta(line.rstrip()))

fullaln = fasta.FastaDict()

# setup keys
for gene in alns[0]:
    fullaln[gene2species(gene)] = ""

# concat sequence
for aln in alns:
    if conf["ungap"]:
        aln = alignlib.removeGappedColumns(aln)
    
    if conf["noconserved"]:
        cons = alignlib.calcConservation(aln)
        ind = util.findneq(1.0, cons)
        aln = alignlib.subalign(aln, ind)
        
    
    #for gene, seq in aln.iteritems():
    #    fullaln[gene2species(gene)] += seq
    
    species2gene = {}
    for gene in aln.keys():
        species2gene[gene2species(gene)] = gene
    
    for sp in fullaln:
        if sp not in species2gene:
            seq = "N" * aln.alignlen()
        else:
            seq = aln[species2gene[sp]]
        
        fullaln[sp] += seq
        

# write full alignment
fullaln.write()
