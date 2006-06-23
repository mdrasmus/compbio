#!/usr/bin/env python

# python libs
import copy
import os
import sys
import StringIO


# rasmsu libs
from rasmus import alignlib
from rasmus import bionj
from rasmus import clustalw
from rasmus import depend
from rasmus import fasta
from rasmus import muscle
from rasmus import util
from rasmus import env
from rasmus import paml
from rasmus import phylip
from rasmus import phyml
from rasmus import puzzletree
from rasmus import sindirlib
from rasmus import treelib




options = [
  ["p:", "prog=", "prog", "<program1>,<program2>,...",
    {"req": True,
     "single": True}],
  ["t:", "usertree=", "usertree", "<user supplied tree",
    {"single": True,
     "default": None,
     "parser": treelib.readTree}],
  ["b:", "bootiter=", "bootiter", "<# iterations>",
    {"parser": int,
     "default": 1,
     "single": True}],
  ["a:", "args=", "args", "<prog args>", 
    {"default": []}],
  ["v", "verbose", "verbose", "",
   {"single": True}],
  
  "Distributed arguments",
  ["g:", "groupsize=", "groupsize", "<group size>",
    {"single": True,
     "default": 20,
     "parser": int}],
  ["n:", "nproc=", "nproc", "<max number of processors>",
    {"single": True,
     "default": 10,
     "parser": int}],
  ["f", "force", "force", "",
    {"help": "force execution",
     "single": True}],

  "Output extensions",     
  ["F:", "fastaext=", "fastaext", "<fasta extension",
    {"default": [".fa", ".fasta"]}],
  ["A:", "alignext=", "alignext", "<align extension>",
    {"default": [".aln", ".afa", ".align"]}],
  ["T:", "treeext=", "treeext", "<tree extension>",
    {"default": [".tree"]}],
  ["U:", "usertreeext=", "usertreeext", "<user tree extension>",
    {"single": True,
     "default": None}],
  ["D:", "distext=", "distext", "<distance matrix extension>",
    {"default": [".dist"]}],
  ["O:", "extraoutput=", "extraoutput", "",
    {"single": True,
     "default": ""}]
    
]

conf = util.parseOptions(sys.argv, options, 
                         resthelp="<alignments> ...", quit=True)


fasta2alignProgs = ["clustalw", "muscle"]

align2treeProgs = ["proml", "dnaml", "protpars", "dnapars",
             "phyml_dna", "phyml_pep"]
             
dist2treeProgs = [ "bionj", "lse" ]

align2distProgs = ["dnadist", "protdist", "fourfold", "ds", "dn", "lapd",
                   "puzzledist"]



# filenames
def getAlignFile(conf, basename):
    return basename + conf["alignext"][-1]

def getDistFile(conf, basename):
    return basename + conf["distext"][-1]

def getTreeFile(conf, basename):
    return basename + conf["treeext"][-1]

def getUserTreeFile(conf, basename):
    return basename + conf["usertreeext"]

def getLabelFile(conf, basename):
    return getAlignFile(conf, basename)


def getUserTree(conf, basename):
    if conf["usertreeext"] != None:
        return treelib.readTree(getUserTreeFile(conf, basename))
    else:
        return conf["usertree"]


def getFileType(conf, infile):
    """Determine the file type of 'infile'"""
    
    for ext in conf["fastaext"]:
        if infile.endswith(ext):
            return "fasta", infile.replace(ext, "")
    
    for ext in conf["alignext"]:
        if infile.endswith(ext):
            return "align", infile.replace(ext, "")

    for ext in conf["treeext"]:
        if infile.endswith(ext):
            return "tree", infile.replace(ext, "")
    
    for ext in conf["distext"]:
        if infile.endswith(ext):
            return "dist", infile.replace(ext, "")
    
    raise "unknown file type '%s'" % infile



def run(conf, infile):
    # determine infile type
    infileType, basename = getFileType(conf, infile)

    # parse common options
    if conf["extraoutput"] != "":
        conf["extraoutput"] = basename + conf["extraoutput"]
    
    if infileType == "fasta":
        fastafile = infile
        alignfile = getAlignFile(conf, basename)
        distfile  = getDistFile(conf, basename)
        labelfile = getLabelFile(conf, basename)
    elif infileType == "align":
        fastafile = None
        alignfile = infile
        distfile  = getDistFile(conf, basename)
        labelfile = getLabelFile(conf, basename)
    elif infileType == "dist":
        fastafile = None
        alginfile = None
        distfile  = infile
        labelfile = getLabelFile(conf, basename)
        
        if not os.path.exists(labelfile):
            labelfile = None
    
    progs = conf["prog"].split(",")
    
    for prog, args in zip(progs, conf["args2"]):
        conf["args"] = args
    
        if prog in fasta2alignProgs:
            assert fastafile != None, "fasta required"
            fasta2align(conf, prog, fastafile, basename)
        
        elif prog in align2treeProgs:
            assert alignfile != None, "alignment required"
            align2tree(conf, prog, alignfile, basename)
        
        elif prog in align2distProgs:
            assert alignfile != None, "alignment required"
            align2dist(conf, prog, alignfile, basename)
        
        elif prog in dist2treeProgs:
            dist2tree(conf, prog, distfile, labelfile, basename)
            
        else:
            raise "unknown program '%s'" % prog


def fasta2align(conf, prog, fastafile, basename):
    alignfile = getAlignFile(conf, basename)
    seqs = fasta.readFasta(fastafile)
    
    if prog == "muscle":
        aln = muscle.muscle(seqs, verbose=conf["verbose"])
        
    elif prog == "clustalw":
        aln = clustalw.clustalw(seqs, verbose=conf["verbose"])
    
    aln.write(alignfile)
    
    
def align2tree(conf, prog, alignfile, basename):
    treefile = getTreeFile(conf, basename)
    distfile = getDistFile(conf, basename)
    aln = fasta.readFasta(alignfile)
    usertree = getUserTree(conf, basename)
    
        
    if prog in ["proml", "dnaml", "protpars", "dnapars"]:
        
        tree = phylip.align2tree(prog,
                                 aln,
                                 usertree=usertree, 
                                 args=conf["args"],
                                 bootiter=conf["bootiter"],
                                 verbose=conf["verbose"],
                                 saveOutput=conf["extraoutput"])
        
        if conf["bootiter"] > 1:
            tree = phylip.consense(tree, verbose=conf["verbose"])
    
    elif prog in ["phyml_dna", "phyml_pep"]:
        if prog == "phyml_dna":
            seqtype = "dna"
        else:
            seqtype = "pep"
    
        tree = phyml.phyml(aln, 
                           usertree=usertree, 
                           verbose=conf["verbose"], 
                           seqtype=seqtype,
                           bootiter=conf["bootiter"],
                           args=conf["args"],
                           saveOutput=conf["extraoutput"])
    

    else:
        raise "unknown phylogeny program '%s'" % prog
    
    
    tree.writeNewick(treefile)


def dist2tree(conf, prog, distfile, labelfile, basename):
    labels, mat = phylip.readDistMatrix(distfile)
    treefile = getTreeFile(conf, basename)
    usertree = getUserTree(conf, basename)
    
    
    if labelfile != None:
        labels = fasta.readFasta(labelfile).keys()
    
    if prog == "bionj":
        tree = bionj.bionj(labels=labels, distmat=mat, verbose=conf["verbose"])
        
    elif prog == "lse":
        if usertree == None:
            raise "Must supply usertree with 'lse'"
        
        sindirlib.setTreeDistances({"debug": 2}, 
                                   usertree, mat, labels)
        tree = usertree
    else:
        raise "unknown program '%s'" % prog
    
    tree.writeNewick(treefile)
    

def align2dist(conf, prog, alignfile, basename):
    distfile = getDistFile(conf, basename)
    aln = fasta.readFasta(alignfile)

    
    if prog == "dnadist":
        phylip.dnadist(aln, distfile, 
                       verbose=conf["verbose"], 
                       args=conf["args"])
    
    elif prog == "protdist":
        phylip.protdist(aln, distfile, 
                        verbose=conf["verbose"], 
                        args=conf["args"])
    
    elif prog == "fourfold":
        mat = alignlib.calcFourFoldDistMatrix(aln)
        phylip.writeDistMatrix(mat, out=distfile)
    
    elif prog == "ds":
        dn, ds = paml.dndsMatrix(aln, verbose=conf["verbose"])
        phylip.writeDistMatrix(ds[1], out=distfile)

    elif prog == "dn":
        dn, ds = paml.dndsMatrix(aln, verbose=conf["verbose"])
        phylip.writeDistMatrix(dn[1], out=distfile)
    
    elif prog == "lapd":
        if conf["args"] == None:
            args = ""
        else:
            args = conf["args"]
        
        os.system("lapd -sd %s '%s' > '%s'" % (args, alignfile, distfile))
    
    elif prog == "puzzledist":
        puzzletree.getDistMatrix(aln, args=conf["args"], output=distfile,
                                 verbose=conf["verbose"])
    
    else:
        raise "unknown program '%s'" % prog

    


def main(conf):
    # parse conf
    files = conf[""]
    
    # save arguments for each program
    allArgs = conf["args"]
    progs = conf["prog"].split(",")
    
    while len(allArgs) < len(progs):
        allArgs.append(" ")
    
    for i in range(len(allArgs)):
        if len(allArgs[i].replace(" ", "")) == 0:
            allArgs[i] = None
    
    conf["args2"] = allArgs
    
    
    if len(files) > conf["groupsize"] and depend.hasLsf():
        # distribute the work

        pipeline = depend.Pipeline("phyloall")
        pipeline.setMaxNumProc(conf["nproc"])
        
        prefix = " ".join(sys.argv[:-len(files)]) + " "
        jobs = []
        
        for i in range(0, len(files), conf["groupsize"]):
            jobs.append(pipeline.add("job%d" % i, 
                         prefix + " ".join(files[i:i+conf["groupsize"]])))
        pipeline.add("all", "echo all jobs complete", jobs)
        
        if conf["force"]:
            pipeline.reset()
        
        pipeline.run("all")
        pipeline.process()
    else:
        for f in files:
            util.tic("%s on %s" % (conf["prog"], f))
            run(conf, f)
            util.toc()
        
        
        
    
    
main(conf)
