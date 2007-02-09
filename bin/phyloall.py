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
from rasmus import env
from rasmus import fasta
from rasmus import mrbayes
from rasmus import muscle
from rasmus import paml
from rasmus import phylip
from rasmus import phylo
from rasmus import phyml
from rasmus import puzzletree
from rasmus import tablelib
from rasmus import treelib
from rasmus import util



options = [
  ["p:", "prog=", "prog", "<program1>,<program2>,...",
    {"single": True,
     "help": "use --proghelp to see all supported programs"}],
  ["t:", "usertree=", "usertree", "<user supplied tree>",
    {"single": True,
     "default": None,
     "parser": treelib.readTree}],
  ["", "opttree=", "opttree", "",
    {"single": True,
     "default": True,
     "parser": util.str2bool,
     "help": "optimize tree topology"}],
  ["b:", "bootiter=", "bootiter", "<# iterations>",
    {"parser": int,
     "default": 1,
     "single": True}],
  ["a:", "args=", "args", "prog:<prog args>", 
    {"default": []}],
  ["v", "verbose", "verbose", "",
   {"single": True}],
  ["m:", "minsize=", "minsize", "<minimum gene family size>",
    {"single": True,
     "default": 3,
     "parser": int,
     "help": "minimum gene family size to reconstruct"}],
  ["M:", "maxsize=", "maxsize", "<maximum gene family size>",
    {"single": True,
     "default": 1000000000, # no limit (effectively)
     "parser": int,
     "help": "maximum gene family size to reconstruct"}],
  ["", "proghelp", "proghelp", "",
    {"single": True}],
  ["r", "resume", "resume", "",
    {"help": "skip over files that already exist",
     "single": True}],    
  ["i", "stdin", "stdin", "",
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
  ["", "statusdir=", "statusdir", "<status directory>",
    {"single": True, 
     "default": "phyloall"}],

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
  ["O:", "extraoutputext=", "extraoutputext", "",
    {"single": True,
     "default": ""}],

  "Misc options",
  ["", "status=", "status", "<ext1>,<ext2>,...",
    {"single": True,
     "help": "report the existance of each file with ext1, ext2, etc."}],
  ["", "stats=", "stats", "<stats file>",
    {"single": True,
     "help": "report the basic stats about alignments and trees"}],

]

conf = util.parseOptions(sys.argv, options, 
                         resthelp="<alignments> ...", quit=True)


# supported programs
fasta2alignProgs = ["clustalw", "muscle"]

align2treeProgs = ["proml", "dnaml", "protpars", "dnapars",
             "phyml_dna", "phyml_pep", "mrbayes_dna", "mrbayes_pep"]
             
dist2treeProgs = ["bionj", "lse", "wls", "nj"]

align2distProgs = ["dnadist", "protdist", "fourfold", "ds", "dn", "lapd",
                   "puzzledist"]

tree2distProgs = ["tree2dist"]



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
            return "fasta", util.replaceExt(infile, ext, "")
    
    for ext in conf["alignext"]:
        if infile.endswith(ext):
            return "align", util.replaceExt(infile, ext, "")

    for ext in conf["treeext"]:
        if infile.endswith(ext):
            return "tree", util.replaceExt(infile, ext, "")
    
    for ext in conf["distext"]:
        if infile.endswith(ext):
            return "dist", util.replaceExt(infile, ext, "")
    
    raise "unknown file type '%s'" % infile


def checkFamilySize(conf, size, filename):
    if size < conf["minsize"] or size > conf["maxsize"]:
        print "skipping '%s'; family size %d outside [%d, %d]" % \
            (filename, size, conf["minsize"], conf["maxsize"])
        return False
    return True


def checkFileExists(conf, infilename, outfilename):
    if infilename == None:
        raise Exception("required file not given")
    
    if not os.path.exists(infilename):
        print "skipping '%s'; input file does not exist" % infilename
        return False
    elif os.path.exists(outfilename) and conf["resume"]:
        print "skipping '%s'; output file already exists" % outfilename
        return False
    else:
        return True


def getDataFiles(conf, infile):
    infileType, basename = getFileType(conf, infile)

    if infileType == "fasta":
        fastafile = infile
        alignfile = getAlignFile(conf, basename)
        distfile  = getDistFile(conf, basename)
        labelfile = getLabelFile(conf, basename)
        treefile  = getTreeFile(conf, basename)
    elif infileType == "align":
        fastafile = None
        alignfile = infile
        distfile  = getDistFile(conf, basename)
        labelfile = getLabelFile(conf, basename)
        treefile  = getTreeFile(conf, basename)
    elif infileType == "dist":
        fastafile = None
        alignfile = None
        distfile  = infile
        labelfile = getLabelFile(conf, basename)
        treefile  = getTreeFile(conf, basename)
        
        if not os.path.exists(labelfile):
            labelfile = None
    
    return fastafile, alignfile, distfile, labelfile, treefile
    

def run(conf, infile, test=False):
    # determine infile type
    infileType, basename = getFileType(conf, infile)

    # parse common options
    if conf["extraoutputext"] != "":
        conf["extraoutput"] = basename + conf["extraoutputext"]
    else:
        conf["extraoutput"] = ""
    
    
    fastafile, alignfile, distfile, labelfile, treefile = \
        getDataFiles(conf, infile)
    if not test:
        util.logger("fasta:", fastafile)
        util.logger("align:", alignfile)
        util.logger("dist: ", distfile)
        util.logger("label:", labelfile)
        util.logger("tree: ",  treefile)
    
    progs = conf["prog"].split(",")
    
    # set skip flag
    executed = False
    
    # run each program
    for prog in progs:
        conf["args"] = conf["argsLookup"][prog]
        
        if prog in fasta2alignProgs:
            if not checkFileExists(conf, fastafile, alignfile): 
                continue
            if not test:
                fasta2align(conf, prog, fastafile, basename)
        
        elif prog in align2treeProgs:
            if not checkFileExists(conf, alignfile, treefile): 
                continue
            if not test:
                align2tree(conf, prog, alignfile, basename)
        
        elif prog in align2distProgs:
            if not checkFileExists(conf, alignfile, distfile): 
                continue
            if not test:
                align2dist(conf, prog, alignfile, basename)
        
        elif prog in dist2treeProgs:
            if not checkFileExists(conf, distfile, treefile): 
                continue
            if not test:
                dist2tree(conf, prog, distfile, labelfile, basename)
        
        elif prog in tree2distProgs:
            if not checkFileExists(conf, treefile, distfile): 
                continue
            if not test:
                tree2dist(conf, prog, treefile, labelfile, distfile, basename)
            
        else:
            raise "unknown program '%s'" % prog
        
        # if we are here, then a program has run and we should not
        # skip this file
        if test:
            return True
        executed = True
    
    return executed
            


def fasta2align(conf, prog, fastafile, basename):
    alignfile = getAlignFile(conf, basename)
    seqs = fasta.readFasta(fastafile)
    
    # check family size
    if not checkFamilySize(conf, len(seqs), fastafile):
        return

    
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
    
    # check family size
    if not checkFamilySize(conf, len(aln), alignfile):
        return

        
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
                           opttree=conf["opttree"],
                           verbose=conf["verbose"], 
                           seqtype=seqtype,
                           bootiter=conf["bootiter"],
                           args=conf["args"],
                           saveOutput=conf["extraoutput"])
    
    elif prog in ["mrbayes_dna", "mrbayes_pep"]:
        if prog == "mrbayes_dna":
            seqtype = "dna"
        else:
            seqtype = "pep"
        
        tree = mrbayes.mrbayes(aln, 
                               seqtype=seqtype,
                               verbose=conf["verbose"],
                               saveOutput=conf["extraoutput"])
        
    else:
        raise "unknown phylogeny program '%s'" % prog
    
    tree.writeNewick(treefile)



def dist2tree(conf, prog, distfile, labelfile, basename):
    labels, mat = phylip.readDistMatrix(distfile)
    treefile = getTreeFile(conf, basename)
    usertree = getUserTree(conf, basename)
    
    # check family size
    if not checkFamilySize(conf, len(labels), distfile):
        return

    
    if labelfile != None:
        labels = fasta.readFasta(labelfile).keys()
    
    if prog == "bionj":
        tree = bionj.bionj(labels=labels, distmat=mat, 
                           verbose=conf["verbose"])
        
    elif prog == "lse":
        if usertree == None:
            raise "Must supply usertree with 'lse'"
        
        lse = phylo.leastSquareError(usertree, mat, labels, weighting=False)
        tree = usertree
    
    elif prog == "wls":
        if usertree == None:
            raise "Must supply usertree with 'wls'"
        
        wls = phylo.leastSquareError(usertree, mat, labels, weighting=True)
        tree = usertree

    elif prog == "nj":
        tree = phylo.neighborjoin(mat, labels, usertree=usertree)
    
    else:
        raise "unknown program '%s'" % prog
    
    tree.writeNewick(treefile)
    


def align2dist(conf, prog, alignfile, basename):
    distfile = getDistFile(conf, basename)
    aln = fasta.readFasta(alignfile)
    
    # check family size
    if not checkFamilySize(conf, len(aln), alignfile):
        return

    
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


def tree2dist(conf, prog, treefile, labelfile, distfile, basename):
    tree = treelib.readTree(treefile)
    labels = fasta.readFasta(labelfile).keys()
    
    if prog == "tree2dist":
        mat = phylo.tree2distmat(tree, labels)
        
    phylip.writeDistMatrix(mat, out=distfile)
    

def displayHelp():
    print >>sys.stderr, "SUPPORTED PROGRAMS\n"

    print >>sys.stderr, "Fasta to alignment"
    print >>sys.stderr, " ", " ".join(fasta2alignProgs)
    print >>sys.stderr

    print >>sys.stderr, "Alignment to tree"
    print >>sys.stderr, " ", " ".join(align2treeProgs)
    print >>sys.stderr

    print >>sys.stderr, "Alignment to distance matrix"
    print >>sys.stderr, " ", " ".join(align2distProgs)
    print >>sys.stderr

    print >>sys.stderr, "Distance matrix to tree"
    print >>sys.stderr, " ", " ".join(dist2treeProgs)
    print >>sys.stderr    

    print >>sys.stderr, "Tree to distance matrix"
    print >>sys.stderr, " ", " ".join(tree2distProgs)
    print >>sys.stderr    
    


def reportStatus(conf, infiles):
    exts = conf["status"].split(",")
    
    tab = tablelib.Table(headers=["name"])
    for ext in exts:
        tab.headers.append(ext)
    
    for infile in infiles:
        infileType, basename = getFileType(conf, infile)
        
        row = {"name": basename}
        
        for ext in exts:
            row[ext] = os.path.exists(basename + ext)
        tab.append(row)
    
    tab.write(sys.stdout)
    return


def reportStats(conf, statsfile, infiles):
    """Report basic stats about alignments and trees"""
    
    ## TODO: make tables stream
    
    stats = tablelib.Table(headers=["name", 
                                    "nseqs",
                                    "alignlen", 
                                    "alignlen_ungapped",
                                    "median_genelen",
                                    "percid",
                                    "treelen"])
    
    util.tic("reporting stats")
    
    for infile in infiles:
        fastafile, alignfile, distfile, labelfile, treefile = \
            getDataFiles(conf, infile)
        infileType, basename = getFileType(conf, infile)
        
        util.logger(basename)
        
        
        # default values
        nseqs = 0
        alignlen = 0
        alignlen_ungapped = 0
        percid = 0.0
        treelen = 0.0
        
        # get align stats
        if os.path.exists(alignfile):
            seqs = fasta.readFasta(fastafile)
            aln = fasta.readFasta(alignfile)
            aln2 = alignlib.removeGappedColumns(aln)
            nseqs = len(aln)
            alignlen = aln.alignlen()
            alignlen_ungapped = aln2.alignlen()
            median_genelen = stats.median(map(len, seqs.itervalues()))
            cons = alignlib.calcConservation(aln)
            percid = util.counteq(1.0, cons) / float(alignlen)
        
        # get tree stats
        if os.path.exists(treefile):
            tree = treelib.readTree(treefile)
            treelen = sum(x.dist for x in tree)
            if nseqs > 0:
                assert nseqs == len(tree.leaves())
            else:
                nseqs = len(tree.leaves())
        
        # record stats
        stats.add(name = os.path.basename(basename),
                  nseqs = nseqs,
                  alignlen = alignlen,
                  alignlen_ungapped = alignlen_ungapped,
                  median_genelen =  median_genelen,
                  percid = percid,
                  treelen = treelen)
    
    stats.write(statsfile)
    
    util.toc()
        
                  
        

def parseArguments(conf):
    """Determines phylo-program specific arguments"""

    # save arguments for each program
    progs = conf["prog"].split(",")
    args = conf["args"]
    conf["argsLookup"] = {}

    for arg in args:
        if ":" not in arg:
            Exception("must specify program name. prog:args")
        i = arg.index(":")
        prog = arg[:i]
        arg = arg[i+1:]
        
        conf["argsLookup"][prog] = arg
    
    # set default arguments
    for prog in progs:
        if prog not in conf["argsLookup"]:
            conf["argsLookup"][prog] = None


def main(conf):
    # print program help    
    if conf["proghelp"]:
        displayHelp()
        return 1
    

    # determine input files
    files2 = copy.copy(conf["REST"])
    if conf["stdin"]:
        for line in sys.stdin:
            files2.append(line.rstrip())
    

    # report file status
    if "status" in conf:
        reportStatus(conf, files2)
        return 0
    
    # report stats
    if "stats" in conf:
        reportStats(conf, conf["stats"], files2)
        return 0
    
    
    if "prog" not in conf:
        sys.stderr.write("phyloall.py: no programs (--prog) given.  Quiting...")
        return 1
    
    parseArguments(conf)
    
    # filter input files that are already complete
    files = filter(lambda x: run(conf, x, test=True), files2)
    
    util.logger("will process %d input files" % len(files))
    util.logger("will skip %d input files" % (len(files2) - len(files)))
    
    
    util.tic("phyloall")
    if len(files) > conf["groupsize"] and depend.hasLsf():
        # distribute the work
        util.logger("distributing work")

        pipeline = depend.Pipeline(conf["statusdir"])
        pipeline.setMaxNumProc(conf["nproc"])
        
        args = copy.copy(sys.argv)
        if len(conf["REST"]) > 0:
            args = args[:-len(conf["REST"])]
        
        # remove arguments that give dispatched jobs trouble
        if "-i" in args:
            args.remove("-i")
        if "--stdin" in args:
            args.remove("--stdin")
        prefix = " ".join(args) + " "
        jobs = []
        
        
        for i in range(0, len(files), conf["groupsize"]):
            jobs.append(pipeline.add("job%d" % i, 
                         prefix + " ".join(files[i:i+conf["groupsize"]])))
        pipeline.add("all", "echo all jobs complete", jobs)
        
        
        pipeline.reset()
        
        pipeline.run("all")
        pipeline.process()
    else:
        # run locally
        
        for f in files:
            util.tic("%s on %s" % (conf["prog"], f))
            run(conf, f)
            util.toc()
    
    util.toc()
    
    return 0
        
        
        
    
if __name__ == "__main__":
    sys.exit(main(conf))
