#!/usr/bin/env python

#!/usr/bin/python

import sys
from rasmus import util

options = [
  ["g:", "genomes=", "genomes", "AUTO<genomes>"],
  ["n:", "name=", "name", "AUTO<prefix name>"],
  ["p:", "pairs=", "pairs", "AUTO<genome pairs file>"],
  
  ["a", "all", "all", "AUTO"],
  ["b:", "blastp=", "blastp", "AUTO<blast options>"],
  ["m", "matches", "matches", "AUTO"],
  ["s:", "synteny=", "synteny", "AUTO<synteny options>"],

  ["P", "part", "part", "AUTO"],
  ["r:", "report=", "report", "AUTO<report directory>"],
  ["M", "reportMatches", "reportMatches", "AUTO"],
  ["A:", "buildAlign=", "buildAlign", "AUTO<options>"],
  ["T:", "buildTrees=", "buildTrees", "AUTO<options>"]
]


try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


logfile = "synphyl.log"

genomes = param["genomes"][-1].split(",")
prefix = param["name"][-1]

# determine genome pairs
if "pairs" in param:
    pairs = util.readDelim(param["pairs"][-1])
else:
    pairs = []
    for i in range(len(genomes)):
        for j in range(i+1, len(genomes)):
            pairs.append([genomes[i], genomes[j]])





def mkblastpfile(genome1, genome2):
    return genome1 + "_" + genome2 + ".blastp"

def mkmatchfile(genome1, genome2):
    return genome1 + "_" + genome2 + ".match"
    
def getMatchFiles(genomes):
    matches = []
    
    if "matches" in param:
        # if we make our own matches then create all possible
        for i in range(len(genomes)):
            for j in range(i, len(genomes)):
                matches.append(mkmatchfile(genomes[i], genomes[j]))
    else:
        # we only use matches that are in our pairs
        for genome1, genome2 in pairs:
            matches.append(mkmatchfile(genome1, genome2))
    return matches


def logPrint(msg):
    print "\techo `date` %s >> %s" % (msg, logfile)

def logStart(msg):
    print "\techo `date` start %s >> %s" % (msg, logfile)

def logEnd(msg):
    print "\techo `date` end   %s >> %s" % (msg, logfile)


print "#"
print "# synphyl Makefile"
print "#"
print

if "blastp" in param or "all" in param:
    print
    print "# make blastps"
    
    if param["blastp"][-1] == "":
        opt = "-e 1e-3"
    else:
        opt = param["blastp"][-1]
    
    blastps = []
    for genome1, genome2 in pairs:
        name = mkblastpfile(genome1, genome2)

        print "%s.blastp: " % (name)
        logStart(name)
        print "\tblastall -p blastp %s -i %s.fasta -d %s.fasta -o %s" % \
            (opt, genome1, genome2, name)
        logEnd(name)

    print "all_blastp: ", " ".join(blastps)
    

if "matches" in param or "all" in param:
    print
    print "# make matches"
    
    matches = []
    for genome1, genome2 in pairs:
        name = mkmatchfile(genome1, genome2)
        name2 = mkblastpfile(genome1, genome2)
        matches.append(name)

        print "%s.match: %s.blastp" % (name, name)
        logStart(name)        
        print "\tblastp2match.py -n 60 %s" % name2
        logEnd(name)
    print "all_match: ", " ".join(matches)


            
if "synteny" in param or "all" in param:
    print
    print "# make synteny"
    
    if param["synteny"][-1] == '':
        param["synteny"][-1] = "-S 1000 -b 1"
    
    cmd = ""
    deps = ""
    for genome1, genome2 in pairs:
        cmd += " -g %s -g %s -m %s \\\n" % \
            (genome1, genome2, mkmatchfile(genome1, genome2))
        deps += " %s" % mkmatchfile(genome1, genome2)
            
    print prefix + ".refine.syncomps: " + deps
    logStart(prefix+".refine.syncomps")
    print "\tfindsynteny.py -f 0,1,2 %s -o %s \\\n%s" % \
        (param["synteny"][-1], prefix, cmd)
    logEnd(prefix+".refine.syncomps")



if "part" in param or "all" in param:
    print 
    print "# make master part file"
    
    print "%s.bbh.part:" % prefix    
    logStart("%s.bbh.part" % prefix)
    bbhs = ""
    for genome1, genome2 in pairs:
        print "\tbestbidir.py -g %s.genes -m %s_%s.match > %s_%s.bbh.match" % \
            (prefix, genome1, genome2, genome1, genome2)
        bbhs += "-m %s_%s.bbh.match " % (genome1, genome2)
    print "\tunionpart.py "+bbhs+" > %s.bbh.part" % prefix
    logEnd("%s.bbh.part" % prefix)
    
    print
    print "%s.syn.bbh.part: %s.refine.syncomps %s.bbh.part " % \
        (prefix, prefix, prefix)
    logStart("%s.syn.bbh.part" % prefix)
    print "\tunionpart.py -p %s.refine.syncomps -p %s.bbh.part > %s.syn.bbh.part" % \
        (prefix, prefix, prefix)
    logEnd("%s.syn.bbh.part" % prefix)
    
    matches = getMatchFiles(genomes)
    cmd = ""
    for match in matches:
        cmd += "-m "+match+" "
    
    print
    print "%s.part: %s.syn.bbh.part" % (prefix, prefix)
    logStart("%s.part" % prefix)
    print "\tgenepull.py -l .5,150 -g %s.genes -O %s.filter.genes -f %s.fasta -p %s.syn.bbh.part %s > %s.part" % \
        (prefix, prefix, prefix, prefix, cmd, prefix)
    logEnd("%s.part" % prefix)


if "report" in param:
    print
    print "# make html report"
    
    print "report: %s.part" % prefix
    logStart("report")
    print "\tpartstats.py -p %s.part -g %s -z -d %s.desc -x %s" % \
        (prefix, ",".join(genomes), prefix, param["report"][-1])
    print "\txsltproc -o %s/main.html %s/partstats.xsl %s/main.xml" % \
        (param["report"][-1], param["report"][-1], param["report"][-1])
    logEnd("report")


if "reportMatches" in param:
    print
    print "# make html matches report"
    
    print "%s-part-matches: %.part" % (prefix, prefix)
    logStart("%s-part-matches" % prefix)
    print "\trm -rf %s-part-matches" % prefix
    print "\tmkdir %s-part-matches" % prefix
    print "\tpartmatches.py -p %s.part -c 0,1,2 -o %s-part-matches \\\n  %s" % \
        (prefix, prefix, " \\\n  ".join(getMatchFiles(genomes)))
    logEnd("%s-part-matches" % prefix)
    
    print
    print "report_matches: %s-part-matches" % prefix
    logStart("report_matches")
    print "\tfind %s-part-matches -type f | xargs matchreport.py -g %s -c 0,1,2 -o %s " % \
        (prefix, ",".join(genomes), param["report"][-1])
    logEnd("report_matches")
    

if "buildAlign" in param:
    if param["buildAlign"][-1] == "":
        opt = ""
    else:
        opt = param["buildAlign"][-1]

    print "%s.aligns: %s.part" % (prefix, prefix)
    logStart("%s.aligns" % prefix)
    print "\tpartalign.py %s -p %s.part -f %s.fasta -o www/%s" % \
        (opt, prefix, prefix, prefix)
    print "\ttouch %s.aligns" % prefix
    print 
    logEnd("%s.aligns" % prefix)


if "buildTrees" in param:
    if param["buildTrees"][-1] == "":
        opt = ""
    else:
        opt = param["buildTrees"][-1]

    print "%s.trees: %s.aligns" % (prefix, prefix)
    logStart("%s.trees" % prefix)
    print "\tbuildtrees.py %s" % opt
    print "\ttouch %s.trees" % prefix
    logEnd("%s.trees" % prefix)
