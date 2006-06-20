#!/usr/bin/env python

from rasmus import algorithms, stats, util, phylip, fasta
from rasmus import ensembl, phyloutil, genomeutil, env
import sys, os
import math, StringIO, copy, random

from compbio.tools import pp

from rasmus import depend





options = [
 ["o:", "outdir=", "outdir", "<output directory>",
    {"help": "directory to place output", 
     "single": True}],
 ["s:", "stree=", "stree", "<species tree>"],
 ["S:", "smap=", "smap", "<gene2species mapping>"],
 ["f:", "ext=", "ext", "<input file extension>"],
 ["t:", "treeext=", "treeext", "<correct tree extension>"],
 ["n:", "num=", "num", "<max number of trees to eval>",
    {"single": True,
     "default": 1000000}],
 ["N:", "nproc=", "nproc", "<max number of processors to use>",
    {"single": True,
     "default": 1000}],
 ["i:", "start=", "start", "<starting tree>", 
    {"single": True,
     "default": 0}],
 ["e:", "exec=", "exec", "<command to exec>"],     
 ["r", "results", "results", "", 
    {"help": "just compute results"}],
 ["g:", "groups=", "groups", "<number of exec per group>",
    {"default": [4]}],
 ["P:", "statusdir=", "statusdir", "<status directory>",
    {"default": "sindir-test-status",
     "single": True}],
 ["L", "local", "local", "",
    {"help": "Do not distribute jobs"}],
 ["F", "force", "force", "",
    {"help": "Force rerun of all jobs"}]
]


conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<input files>")


# Pipeline setup
# decide whether to use default dispatch (LSF,BASH) or force BASH
if "local" in conf:
    pipeline = depend.Pipeline(conf["statusdir"], False, depend.BASH_DISPATCH)
else:
    pipeline = depend.Pipeline(conf["statusdir"])
pipeline.setLogOutput()


#
# filename conventions
#

def getBasenames(conf, infile):
    basename = infile.replace(conf["ext"][-1], "")
    return os.path.dirname(basename), os.path.basename(basename)

def getCorrectTree(conf, infile):
    basedir, basefile = getBasenames(conf, infile)
    return os.path.join(basedir, basefile + conf["treeext"][-1])

def getOutputTree(conf, infile):
    basedir, basefile = getBasenames(conf, infile)
    return os.path.join(conf["outdir"], basefile + ".tree")



def main(conf):
    env.addEnvPaths("DATAPATH")
    
    print "Pipeline is using dispatch: '%s'", pipeline.dispatch
    
    if "results" in conf:
        makeReport(conf)
    else:
        testAll(conf)
        makeReport(conf)


def testAll(conf):
    util.tic("testing")
    
    files = conf[""]
    
    jobs = []
    start = int(conf["start"])
    num   = int(conf["num"])
    for infile in files[start:start+num]:
        jobs.append(runJob(conf, infile))
    jobs = filter(lambda x: x != None, jobs)
    
    groups = pipeline.addGroups("testgroup", jobs, int(conf["groups"][-1]))
    alljobs = pipeline.add("testall", "echo all", groups)
    
    pipeline.reset()
    pipeline.run("testall")
    pipeline.process()
    
    util.toc()


def runJob(conf, infile):
    basedir, basefile = getBasenames(conf, infile)
    
    # skip tests when output tree already exists
    if "force" not in conf and \
       os.path.exists("%s/%s.tree" % (conf["outdir"], basefile)):
        util.log("skip '%s/%s.tree' output exists " % (conf["outdir"], basefile))
        return None
    
    cmd = conf["exec"][-1]
    cmd = cmd.replace("$FILE", basefile)
    cmd = cmd.replace("$DIR", basedir)
    
    jobname = pipeline.add("job_" + basefile, cmd, [])
    
    return jobname


def checkOutput(conf, infile, stree, gene2species):
    basedir, basefile = getBasenames(conf, infile)
    outfile = getOutputTree(conf, infile)
    correctTreefile = getCorrectTree(conf, infile)

    if not os.path.exists(outfile):
        return None, None 

    tree1 = algorithms.readTree(outfile)
    tree2 = algorithms.readTree(correctTreefile)

    tree1 = phyloutil.reconRoot(tree1, stree, gene2species)
    tree2 = phyloutil.reconRoot(tree2, stree, gene2species)

    hash1 = phyloutil.hashTree(tree1)
    hash2 = phyloutil.hashTree(tree2)
    
    return tree1, tree2


def makeReport(conf):
    util.tic("make report")

    gene2species = genomeutil.readGene2species(* map(env.findFile, 
                                                     conf["smap"]))
    stree = algorithms.readTree(env.findFile(conf["stree"][-1]))
    

    outfiles = util.listFiles(conf["outdir"], ".tree")
    infiles = conf[""]
    
    results = []
    counts = util.Dict(1, 0)
    orths = [0, 0, 0, 0]
    
    
    for infile in infiles:
        basedir, basefile = getBasenames(conf, infile)
        tree1, tree2 = checkOutput(conf, infile, stree, gene2species)
        
        if tree1 == None:
            continue
        
        hash1 = phyloutil.hashTree(tree1)
        hash2 = phyloutil.hashTree(tree2)
        
        # make a species hash
        shash1 = phyloutil.hashTree(tree1, gene2species)
        shash2 = phyloutil.hashTree(tree2, gene2species)
        counts[(shash1,shash2)] += 1
        
        if hash1 == hash2:
            results.append([basefile, True])
        else:
            results.append([basefile, False])
        
        
        orths = util.vadd(orths, testOrthologs(tree1, tree2, stree, gene2species))
    
    
    # print final results    
    reportfile = os.path.join(conf["outdir"], "results")
    out = file(reportfile, "w")

    total = len(results)
    ncorrect = util.counteq(True, util.cget(results, 1))
    nwrong = util.counteq(False, util.cget(results, 1))


    print >>out, "total:      %d" % total
    print >>out, "#correct:   %d (%f%%)" % (ncorrect, 100*ncorrect / float(total))
    print >>out, "#incorrect: %d (%f%%)" % (nwrong, 100*nwrong / float(total))
    print >>out
    print >>out


    util.printcols(results, out=out)
    
    # print out shash counts
    mat = []
    items = counts.items()
    items.sort(key=lambda x: -x[1])
    tot = float(sum(counts.values()))
    for (tree1, tree2), num in items:
        if tree1 == tree2: c = "*"
        else: c = " "
        mat.append([c, tree1, num, num/tot])
        mat.append([c, tree2, "", ""])
        mat.append(["", "", "", ""])
    
    util.printcols(mat, out=out)
    
    
    # find ortholog sn, sp
    [tp, fn, fp, tn] = orths
    print >>out
    print >>out, "ortholog detection:"
    print >>out, "sensitivity:", tp / float(tp + fn)
    print >>out, "specificity:", tn / float(fp + tn)

    util.toc()



def testOrthologs(tree1, tree2, stree, gene2species):
    recon1 = phyloutil.reconcile(tree1, stree, gene2species)
    recon2 = phyloutil.reconcile(tree2, stree, gene2species)
    
    orths1 = phyloutil.findAllOrthologs(tree1, stree, recon1)
    orths2 = phyloutil.findAllOrthologs(tree2, stree, recon2)
    
    set1 = util.makeset(map(tuple, orths1))
    set2 = util.makeset(map(tuple, orths2))
    
    ngenes = len(tree2.leaves())
    nonorths = ngenes * (ngenes + 1) / 2 - len(set2)
    
    # sensitivity and specificity
    overlap = util.intersect(set1, set2)
    tp = len(overlap)
    fp = len(set1) - len(overlap)
    fn = len(set2) - len(overlap)
    tn = nonorths - fp
    
    return [tp, fn, fp, tn]



"sindir.py -p $PARAMS -s $STREE -S $SMAP -d $DIST -l $ALIGN -T $CORRECT_TREE -o $OUT_TREE $EXTRA &>$OUT_DEBUG"

"""
def runPhylip(param, testdir, alnfile, gene2species):
    distext = param["distext"][-1]
    labelext = param["labelext"][-1]
    treeext = param["treeext"][-1]
        
    
    basefile = os.path.basename(distfile.replace(distext, ""))

    labelfile = distfile.replace(distext, labelext)
    outtree = os.path.join(param["outdir"], basefile + ".tree")
    correctTreefile = distfile.replace(distext, treeext)

    debugfile = os.path.join(param["outdir"], basefile + ".debug")

    cmd = "sindir.py -p %s -s %s -S %s -d %s -l %s -T %s -o %s %s &>%s" % (
              param["param"], param["stree"], param["smap"][-1], 
              distfile, labelfile, correctTreefile, outtree, 
              param["extra"], debugfile)
    
    jobname = pipeline.add(os.path.basename(distfile), cmd, [])
    
    return jobname
"""




main(conf)
