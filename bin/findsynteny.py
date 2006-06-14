#!/usr/bin/python

import getopt
import sys
from rasmus import util, synteny, genomeio, genomeutil

options = [
    ["g:", "genome=", "genome", "AUTO<genome>"],
    ["m:", "match=", "match", "AUTO<match file>"],
    ["o:", "outprefix=", "outprefix", "AUTO<output prefix>"],
    ["l:", "log=", "log", "AUTO<log file>"],
    ["R", "refine", "refine", "AUTO"],
    ["s:", "synteny=", "synteny", "AUTO<synteny>"],
    ["c:", "syncomps=", "syncomps", "AUTO<syncomps>"],
    ["f:", "fields=", "fields", "AUTO<gene1col>,<gene2col>,<scorecol>"],
    ["n:", "near=", "near", ""],
    ["r:", "range=", "range", ""],
    ["b:", "minblock=", "minblock", "AUTO<min block size>"],
    ["S:", "sigscore=", "sigscore", "AUTO<sig score>"],
    ["d:", "dist=", "dist", "AUTO<distance in number of genes to merge>"]
]


try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


data = []


cols = map(int, param["fields"][-1].split(","))

if "match" in param: 
    for i in range(len(param["match"])):
        data.append([param["genome"][2*i], param["genome"][2*i+1], 
                     genomeio.structfile(param["genome"][2*i]),
                     genomeio.structfile(param["genome"][2*i+1]),
                     param["match"][i]])


if not "outdir" in param:
    param["outdir"] = ["."]

if "log" in param:
    logfile = file(param["log"][-1], "w")
else:
    logfile = file("synteny.log", "w")


conf = synteny.initConf()
conf["cols"] = cols
conf["minscore"] = 60

if "minblock" in param:
    conf['minBlockSize'] = int(param["minblock"][-1])

if "sigscore" in param:
    conf["sigScore"] = float(param["sigscore"][-1])

if "near" in param:
    conf["nearRange"] = int(param["near"][-1])

if "range" in param:
    conf["rangeGenes"] = int(param["range"][-1])


if "refine" not in param:
    if "dist" in param:
        conf["rangeGenes"] = param["dist"]
    
    comps, matching = synteny.findMultiSynteny(
                            data, conf, param["outprefix"][-1], logfile)
    
    synteny.refineSynteny(matching, comps, data, conf, 
                            param["outprefix"][-1], logfile)
else:
    matching = genomeutil.Matching()
    genomes = util.uniq(param["genome"])
    genomeio.readGenomes(matching, genomes)
    comps = util.readDelim(param["syncomps"][-1])
    genomeio.readSimpleMatches(matching, param["syncomps"][-1])
    genomeio.readBlockDimensions(matching, param["synteny"][-1], True)

    synteny.refineSynteny(matching, comps, data, conf, param["outprefix"][-1], logfile)
