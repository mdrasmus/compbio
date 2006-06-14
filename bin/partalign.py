#!/usr/bin/env python

import os
import sys

from rasmus import bionj
from rasmus import clustalw
from rasmus import ensembl
from rasmus import env
from rasmus import fasta
from rasmus import genomeutil
from rasmus import grid
from rasmus import muscle
from rasmus import phylip
from rasmus import stats
from rasmus import util


SCRIPT = sys.argv[0]

options = [
    ["p:", "part=", "part", "<part file>"],
    ["f:", "fasta=", "fasta", "<fasta file>"],    
    ["o:", "outdir=", "outdir", "<output directory>"],
    ["O:", "one2one=", "one2one", "<genome1>,<genome2>,..."],
    ["i:", "index=", "index", "<starting index>"],
    ["e:", "end=", "end", "<ending index>"],
    ["t:", "tree=", "tree", "<filename suffix>"],
    ["d", "dist", "dist", "", 
        {"help": "Produce a distance file"}],
    ["m:", "minsize=", "minsize", "<min size part>"],
    ["S:", "smap=", "smap", "<gene2species map file>"],
    ["P:", "paths=", "paths", "<data path1>:<data path2>:...", 
        {"default": "."}],
    ["T:", "treeprog=", "treeprog", "<tree building program>", 
        {"default": "proml"}]
]


param = util.parseOptions(sys.argv, options, quit=True)


def main(param):
    # setup paths
    env.addEnvPaths("DATAPATH")
    env.addPaths(param["paths"])

    # read species map
    if "smap" in param:
        gene2species = genomeutil.readGene2species(* map(env.findFile, param["smap"]))
    else:
        gene2species = genomeutil.gene2species

    # parse parameters
    if "index" in param:
        param["start"] = int(param["index"][-1])
    else:
        param["start"] = 0
    
    if "end" in param:
        param["end"] = int(param["end"][-1])
    else:
        param["end"] = -1

    if "minsize" in param:
        param["minsize"] = int(param["minsize"][-1])
    else:
        param["minsize"] = 3
    
    
    if "part" in param:
        seqs = fasta.FastaDict()
        nseqs = 0
        util.tic("read sequences")
        for f in param["fasta"]:
            util.tic("reading '%s'" % f)
            seqs.read(env.findFile(f))
            util.log("read %d sequences" % (len(seqs) - nseqs))
            nseqs = len(seqs)
            util.toc()
        util.toc()
        
        parts = util.readDelim(env.findFile(param["part"][-1]))
        partAlign(param, parts, seqs, gene2species)
    else:
        for filename in param["fasta"]:
            runAlign(param, filename)


def buildSubprocParams(param):
    cmd = ""
    
    if "tree" in param:
        cmd += "-t %s " % param["tree"][-1]
        
    if "treeprog" in param:
        cmd += "-T %s " % param["treeprog"][-1]
    
    if "dist" in param:
        cmd += "-d "
    
    return cmd
    

def partAlign(param, parts, seqs, gene2species):
    batch = 10
    batches = []
    jobs = []
    
    if "one2one" in param:
        genomes = tuple(util.sort(param["one2one"][-1].split(",")))
    
    if param["end"] == -1:
        param["end"] = len(parts)

    # code for finishing an alignment job
    def onJobDone(group):
        for job in group.jobs:
            util.log("completed partition %d" % job.partid)


    # submit jobs
    util.tic("create jobs for %d partitions" % len(parts))
    for i in xrange(param["start"], param["end"]):
        part = parts[i]
        
        # filter for one2one is requested
        if "one2one" in param:
            genomes2 = tuple(util.sort(map(gene2species, part)))
            if genomes != genomes2:
                util.log("part %d of %d: SKIP - NOT ONE2ONE" % 
                         (i, len(parts)))
                continue

        
        # ensure we have minimum number of sequences
        seqs2 = util.subdict(seqs, part)        
        if len(seqs2) < param["minsize"]:
            util.log("part %d of %d: SKIP - NOT ENOUGH SEQUENCES %d" % 
                     (i, len(parts), len(seqs2)))
            continue
        
        # create a job
        seqfile = "%s/%d.fasta" % (param["outdir"][-1], i)
        fasta.writeFasta(seqfile, seqs2)
        
        util.log("part %d of %d: create job for %d sequences" % 
                 (i, len(parts), len(part)))
        jobs.append(grid.Job("%s -f %s -o %s %s" %
                             (SCRIPT, seqfile, param["outdir"][-1], 
                             buildSubprocParams(param))))
        jobs[-1].partid = i
        
        
        # batch jobs together
        if len(jobs) >= batch:
            util.log("submit group of %d jobs" % len(jobs))
            batches.append(grid.JobGroup(jobs))
            grid.submitJob(batches[-1])
            jobs = []
        
        
        # finish up done jobs
        done = grid.getDoneJobs(batches)
        for job in done:
            onJobDone(job)
            batches.remove(job)
    util.toc()
    
    # submit unbatched jobs
    util.log("submit group of %d jobs" % len(jobs))
    batches.append(grid.JobGroup(jobs))
    grid.submitJob(batches[-1])
    jobs = []
    
    # wait for all jobs to stop
    util.tic("wait for all jobs to finish")
    grid.waitForJobs(batches, onDone=onJobDone)
    util.toc()
    
    

def runAlign(param, seqFile):
    alnFile = seqFile.replace(".fasta", ".align")
    
    # make alignment
    if not os.path.exists(alnFile):
        try:
            ret = os.system("muscle -in %s -out %s" % (seqFile, alnFile))
        except:
            # sometimes muscle runs out of memory
            # if so, use clustalw instead
            seqs = fasta.readFasta(seqFile)
            aln = clustalw.clustalw(seqs)
            fasta.writeFasta(alnFile, aln)
    
    # make distances
    if "dist" in param:
        distFile = seqFile.replace(".fasta", ".dist")
        if not os.path.exists(distFile):
            aln = fasta.readFasta(alnFile)
            labels = phylip.protdist(aln, distFile)
    
    # make tree
    if "tree" in param:
        treeFile = seqFile.replace(".fasta", "." + param["tree"][-1])
        if not os.path.exists(treeFile):
            aln = fasta.readFasta(alnFile)
            
            if param["treeprog"][-1] == "proml":
                tree = phylip.proml(aln)
                tree.writeNewick(treeFile)
            elif param["treeprog"][-1] == "bionj":
                tree = bionj.bionj(aln)
                tree.writeNewick(treeFile)
            else:
                util.error("Unknown tree program '%s'" % param["treeprog"][-1])


main(param)
