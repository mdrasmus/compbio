#!/usr/bin/env python

# python libs
import sys

# rasmsu libs
from rasmus import alignlib
from rasmus import bionj
from rasmus import fasta
from rasmus import muscle
from rasmus import util
from rasmus import env
from rasmus import phylip
from rasmus import phyml
from rasmus import treelib


options = [
  ["p:", "prog=", "prog", "<phylip program>"],
  ["t:", "usertree=", "usertree", "<user supplied tree",
    {"single": True}],
  ["b:", "bootiters=", "bootiters", "<# iterations>"],
  ["v", "verbose", "verbose", "",
   {"single": True}],
  ["a:", "args=", "args", "<prog args>", 
    {"default": "y",
     "single": True}],
  ["A:", "alignext=", "alignext", "<align extension>"],
  ["T:", "treeext=", "treeext", "<tree extension>"],
  ["D:", "distext=", "distext", "<distance matrix extension>"],
  ["k", "keep", "keep", ""]
]
    
param = util.parseOptions(sys.argv, options, resthelp="<alignments> ...", quit=True)



def run(param, infile):
    if "keep" in param:
        outExtra = infile.replace(param["alignext"][-1], ".output")
    else:
        outExtra = ""
    

    if "usertree" in param:
        usertree = treelib.readTree(param["usertree"])
    else:
        usertree = None
    
       
    
    if param["prog"][-1] == "proml":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile)
    
        if "bootiters" not in param:
            tree = phylip.proml(aln,
                                usertree=usertree, 
                                args=param["args"],
                                verbose=param["verbose"],
                                saveOutput=outExtra)
        else:
            # TODO: out of date
            trees = phylip.bootProml(aln, int(param["bootiters"][-1]), jumble=1,
                                     verbose=param["verbose"])
            tree = phylip.consense(trees, verbose=param["verbose"])


        tree.writeNewick(outfile)
    
    elif param["prog"][-1] == "dnaml":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile)
    
        if "bootiters" not in param:
            tree = phylip.dnaml(aln,
                                usertree=usertree, 
                                args=param["args"], 
                                verbose=param["verbose"],
                                saveOutput=outExtra)
        else:
            # TODO: out of date
            trees = phylip.bootDnaml(aln, int(param["bootiters"][-1]), jumble=1,
                                     verbose=param["verbose"])
            tree = phylip.consense(trees, verbose=param["verbose"])


        tree.writeNewick(outfile)
    
    elif param["prog"][-1] == "phyml_dna":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile) 
        
        if "args" in param:
            tree = phyml.phyml(aln, 
                               verbose=param["verbose"], 
                               seqtype="dna",
                               args=param["args"][-1],
                               saveOutput=outExtra)
        else:
            tree = phyml.phyml(aln, 
                               verbose=param["verbose"], 
                               seqtype="dna",
                               saveOutput=outExtra)
        tree.writeNewick(outfile)
    
    elif param["prog"][-1] == "phyml_pep":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile) 
        tree = phyml.phyml(aln, 
                           verbose=param["verbose"], 
                           seqtype="pep",
                           saveOutput=outExtra)
        tree.writeNewick(outfile)
        
    
    elif param["prog"][-1] == "bionj_pep":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile)
        tree = bionj.bionj(aln, seqtype="pep", verbose=param["verbose"])
        tree.writeNewick(outfile)

    elif param["prog"][-1] == "bionj_dna":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile)
        tree = bionj.bionj(aln, seqtype="dna", verbose=param["verbose"])
        tree.writeNewick(outfile)
    
    elif param["prog"][-1] == "bionj_fourfold":
        outfile = infile.replace(param["alignext"][-1], param["treeext"][-1])
        aln = fasta.readFasta(infile)
        mat = alignlib.calcFourFoldDistMatrix(aln)
        
        tree = bionj.bionj(aln, distmat=mat, 
                           seqtype="dna", verbose=param["verbose"])
        tree.writeNewick(outfile)
    
    
    elif param["prog"][-1] == "dnadist":
        outfile = infile.replace(param["alignext"][-1], param["distext"][-1])
        aln = fasta.readFasta(infile)
        phylip.dnadist(aln, outfile, verbose=param["verbose"])
    
    elif param["prog"][-1] == "protdist":
        outfile = infile.replace(param["alignext"][-1], param["distext"][-1])
        aln = fasta.readFasta(infile)
        phylip.protdist(aln, outfile, verbose=param["verbose"])
    
    else:
        raise "unknown program '%s'" % param["prog"][-1]
    


def main(param):
    # parse param
    files = param[""]
    
    
    for f in files:
        util.tic("%s on %s" % (param["prog"][-1], f))
        run(param, f)
        util.toc()
        
        
        
    
    
main(param)
