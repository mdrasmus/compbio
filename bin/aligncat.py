#!/usr/bin/env python

import os, sys

from rasmus import env, util, fasta, genomeutil, alignlib


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


alns = map(fasta.readFasta, conf["REST"])

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
    
    for gene in fullaln:
        if gene not in aln:
            seq = "N" * aln.alignlen()
        else:
            seq = aln[gene2species(gene)]
        
        fullaln[gene] += seq
        

# write full alignment
fullaln.write()
