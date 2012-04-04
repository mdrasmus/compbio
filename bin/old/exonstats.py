#!/usr/bin/env python


from rasmus import fasta, util, genomeio, genomeutil, exonutil
import sys


options = [
 ["a:", "align=", "align", "AUTO<fasta alignment>"],
 ["b:", "batch=", "batch", "AUTO<genoms>"],
 ["m", "missed", "missed", "AUTO"],
 ["s", "stats", "stats", "AUTO"]
]


try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)



   

def printHeader(fields):
    print "\t".join(fields)

def printRow(stats, fields):
    values = [str(stats[x]) for x in fields]
    print "\t".join(values)



def initStats(genomes):
    stats = {
    "perfect": 0,
    "movement": 0,
    "conserved": 0,
    #"alt": 0,
    "exon gain": 0,
    "exon loss": 0,
    "intron gain": 0,
    "intron loss": 0,
    "total": 0
    }
    
    fields = [
        "file",
        "total",
        "conserved",
        "perfect",
        "movement",
        "exon gain",
        "exon loss",
        "intron gain",
        "intron loss"]
    
    for genome in genomes:
        stats["exon gain (%s)" % genome] = 0
        fields.append("exon gain (%s)" % genome)
    
        stats["exon loss (%s)" % genome] = 0
        fields.append("exon loss (%s)" % genome)
        
        stats["intron gain (%s)" % genome] = 0
        fields.append("intron gain (%s)" % genome)
    
        stats["intron loss (%s)" % genome] = 0
        fields.append("intron loss (%s)" % genome)
        
        stats["exon miss (%s)" % genome] = 0
        fields.append("exon miss (%s)" % genome)
    
    #fields.append("alt")
    #for genome in genomes:
    #    stats["alt (%s)" % genome] = 0
    #    fields.append("alt (%s)" % genome)

    return stats, fields

    

def analyze(filename, aln, comps, exonData, stats = None):
    genomes = [x.genome for x in aln]
    genomes.sort()

    if stats == None:
        stats, fields = initStats(genomes)
    
    stats["total"] += len(comps)
    
    # now report on each aligned column of exons
    for comp in comps:
        # sort by row
        #comp.sort(lambda a,b: cmp(exonData[a]["row"], exonData[b]["row"]))
        
        # prune alternative spliced exons
        comp2 = []
        for exon1 in comp:
            ov = False
            for exon2 in comp2:
                if exonData[exon1]["row"] == exonData[exon2]["row"] and \
                   util.overlap(exon1.start, exon1.end, exon2.start, exon2.end):
                    ov = True
                    print exon1, exon2
            if not ov:
                comp2.append(exon1)
        comp = comp2
        
        rows = exonutil.getExonRows(comp, exonData)
        
        if exonutil.isConserved(comp, exonData, len(aln)):
            stats["conserved"] += 1        
            if exonutil.isPerfectConserved(comp, exonData):
                stats["perfect"] += 1
            else:
                if max([len(x) for x in rows]) == 1:
                    stats["movement"] += 1
                else:
                    print rows, comp
                    # multiple exons in same row
                    
                    # if multiple rows have multiple exons then intron loss
                    if len(filter(lambda x: len(x) > 1, rows)) > 1:
                        stats["intron loss"] += 1
                        
                        for i in xrange(len(rows)):
                            if len(rows[i]) == 1:
                                stats["intron loss (%s)" % aln[i].genome] += 1
                    else:
                        stats["intron gain"] += 1
                        
                        for i in xrange(len(rows)):
                            if len(rows[i]) > 1:
                                stats["intron gain (%s)" % aln[i].genome] += 1

                    
        else:
            if len(filter(lambda x: len(x) > 0, rows)) == 1:
                row = util.find(lambda x: len(x) > 0, rows)[0]
                stats["exon gain"] += 1
                stats["exon gain (%s)" % aln[row].genome] += 1
            else:
                stats["exon loss"] += 1
                
                for i in xrange(len(aln)):
                    if i not in rows:                        
                        # is there sequence there?
                        # find an exon that is in the comp
                        data = exonData[comp[0]]
                        seq = aln[i].seq[data["start"]:data["end"]]
                        ngaps = util.counteq("-", seq)
                        if ngaps < len(seq) / 2 and \
                           exonutil.isWellAligned(comp, exonData):
                            stats["exon miss (%s)" % aln[i].genome] += 1
                            
                            if "missed" in param:
                                localLookup = genomeutil.mkLocalLookup(aln[i].seq)
                                
                                start = genomeutil.align2global(data["start"], 
                                        aln[i].start, aln[i].end, aln[i].strand, 
                                        localLookup)
                                
                                end = genomeutil.align2global(data["end"], 
                                        aln[i].start, aln[i].end, aln[i].strand, 
                                        localLookup)
                                
                                if start > end:
                                    tmp = start; start = end; end = tmp
                                
                                print filename, data["start"], data["end"], \
                                      aln[i].genome, aln[i].chrom, start, end
                                    
                        else:
                            stats["exon loss (%s)" % aln[i].genome] += 1
        
    
    #if max(rows.values()) > 1:
    #    altrows = {}
    #    for i in rows:
    #        if rows[i] > 1:
    #            altrows[i] = rows[i]
    #    stats["alt"] += len(altrows)
    #
    #    for i in altrows:
    #        stats["alt (%s)" % aln[i].genome] += 1

    return stats, fields

    


if "align" in param:
    # read alignment
    util.tic("read input")
    util.tic("read alignment")
    alignFile = param["align"][-1]
    aln = fasta.readAlignment(alignFile)
    util.toc()

    # read exon information
    genomes = {}
    for row in aln:
        util.tic("read %s exons" % row.genome)
        genomes[row.genome] = genomeutil.Genome(row.genome)

        genomeio.readEnsmartExons(genomes[row.genome], 
                                  genomeio.structfile(row.genome, row.chrom))
        genomes[row.genome].autoconf()
        util.toc()
    util.toc()

    comps, exonData = exonutil.findAlignedExons(aln, genomes)
    
    stats, fields = analyze(alignFile, aln, comps, exonData)
    stats["file"] = alignFile
    
    if "stats" in param:
        for key in fields:
            print "%s\t%s" % (key, str(stats[key]))

    
elif "batch" in param:
    # read exon information
    genomes = {}
    for genome in param["batch"][-1].split(","):
        util.tic("read %s exons" % genome)
        genomes[genome] = genomeutil.Genome(genome)

        genomeio.readEnsmartExons(genomes[genome], 
                                  genomeio.structfile(genome))
        genomes[genome].autoconf()
        util.toc()
    
    keys = genomes.keys()
    keys.sort()
    stats, fields = initStats(keys)
    
    if "stats" in param:
        printHeader(fields)
    
    for alignFile in rest:
        # read alignment
        util.log(alignFile)
        aln = fasta.readAlignment(alignFile)
        
        comps, exonData = exonutil.findAlignedExons(aln, genomes)
        stats, fields2 = analyze(alignFile, aln, comps, exonData)
        stats["file"] = alignFile
        
        if "stats" in param:
            printRow(stats, fields)
            sys.stdout.flush()
    
    
