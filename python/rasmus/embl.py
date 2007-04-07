#
# NOTE: this module is pretty much deprecated
#

import sys
import os

from rasmus import util

timer = util.Timer()

# common prefixes in EMBL flat files
chromPrefix    = "AC   "
featurePrefix  = "FT   "
genePrefix     = "FT   gene"
geneNamePrefix = "FT                   /gene="
varPrefix      = "FT                   /"
mrnaPrefix     = "FT   mRNA"
continuePrefix = "FT      "
seqPrefix      = " "
goPrefix       = "FT                   /db_xref=\"GO:"

class Gene:
    name = ""
    start = 0
    end = 0
    coords = ""
    
    def __init__(self, name, coords):
        self.name   = name
        self.coords = coords

class Chrom:
    name = ""
    start = 0
    end = 0
    
    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end


class Seeker:
    files = []
    index = None
    
    def __init__(self, filenames, index):
        self.files = {}
        self.index = index
        
        for f in filenames:
            self.files[f] = file(f)
    
    def seek(self, gene):
        infile = self.files[self.index[gene][0]]
        infile.seek(self.index[gene][1])
        return infile


def prunePrefix(line):
    return line[21:]

def seekNextGene(infile):
    for line in util.SafeReadIter(infile):
        if line.startswith(genePrefix):
            # backup to beginning of line
            infile.seek(-len(line), 1)
            return True
    return False

def seekNextChrom(infile):
    for line in util.SafeReadIter(infile):
        if line.startswith(chromPrefix):
            # backup to beginning of line
            infile.seek(-len(line), 1)
            return True
    return False


def readGene(infile, offset = 0):
    coordsStr = ""
    
    # read through gene section
    for line in util.SafeReadIter(infile):
        if line.startswith(geneNamePrefix):
            name = line[len(geneNamePrefix):].rstrip()
            break
        else:
            coordsStr += prunePrefix(line.rstrip())
    
    # parse coords
    coords = parseCoords(coordsStr, offset)
    
    return Gene(name, coords)


def readChrom(infile):
    line = infile.readline()
    return parseChrom(line)
    

def parseChrom(line):
    (chrom, start, end) = ("NONE", 0, 0)
    try:
        (chrom, start, end) = line.split(":")[2:5]
    except ValueError:
        if len(line.split(":")) == 1:
            chrom = line.replace(chromPrefix, "").rstrip()
            start = 1
            end = 1
        else:
            raise ValueError
    return Chrom(chrom, start, end)  


def readGos(infile):
    gos = []
    
    # skip first line
    infile.readline()
    
    # iter through rest of gene
    for line in util.SafeReadIter(infile):
        if line.startswith(genePrefix) or \
           not line.startswith(featurePrefix):
            break
        if line.startswith(goPrefix):
            gos.append(line[len(goPrefix):].replace("\n", "").replace("\"", ""))
    return gos

def readExons(infile, offset = 0):
    exonList = []
    
    # skip gene start
    infile.readline()
    
    # keep track whether we're in mRNA section
    mrna = False
    coords = ""
    
    for line in util.SafeReadIter(infile):
        if line.startswith(genePrefix) or \
            not line.startswith(featurePrefix):
                break
        if line.startswith(mrnaPrefix):
            mrna = True
            coords += prunePrefix(line).rstrip()
        elif mrna:
            if not line.startswith(continuePrefix):
                exonList.append(parseCoords(coords, offset))
                mrna = False
                coords = ""
            elif not line.startswith(varPrefix):
                coords += prunePrefix(line.rstrip())
    return exonList

        
def parseCoords(expr, start):
    coords = []
    
    if expr[0:5] == "join(":
        end = expr.rfind(")")
        for exp in expr[5:end].split(","):
            coords += parseCoord(exp, start)
    else:
        coords += parseCoord(expr, start)
    return coords


def parseCoord(expr, start, flip = 1):
    if expr == "":
        return []

    if expr[0:11] == "complement(":
        end = expr.rfind(")")
        return parseCoord(expr[11:end], start, -1 * flip)
    elif expr[0] in "0123456789":
        tokens = expr.split("..")
        return [int(tokens[0]) + start, int(tokens[1]) + start, flip]
    else:
        return []


def readIndex(filename):
    index = {}
    for line in file(filename):
        (gene, filename, pos) = line.split()
        index[gene] = (filename, int(pos))
    return index


def readChromLookup(files, chromIndexFile):
    # read in chrom indexes
    lookup = {}
    for f in files:
        lookup[f] = []
    
    # build chrom lookup table
    for line in file(chromIndexFile):
        (chrom, start, end, filename, pos) = line.split()
        lookup[filename].append((int(pos), chrom, int(start), int(end)))
    return lookup


def readChromIndex(files, chromIndexFile):
    # read in chrom indexes
    index = {}
    
    class ChromIndexEntry:
        def __init__(self, start, end, filename, pos):
            self.start = start
            self.end = end
            self.filename = filename
            self.pos = pos
            
    
    # build chrom lookup table
    for line in file(chromIndexFile):
        (chrom, start, end, filename, pos) = line.split()
        util.case(index, chrom, []).append(
            ChromIndexEntry(int(start), int(end), filename, int(pos)))
    
    # sort index
    for chrom in index:
        index[chrom].sort(lambda a,b: a.start - b.start)
    
    return index
    

def makeGeneIndexFile(emblFilenames, indexFilename):
    outfile = file(indexFilename, "w")

    for emblFilename in emblFilenames:
        infile = file(emblFilename)


        timer.start("make index for " + emblFilename)
        while seekNextGene(infile):
            index = infile.tell()
            gene = readGene(infile)
            print >>outfile, gene.name, emblFilename, index
        timer.stop()

def makeChromIndexFile(emblFilenames, indexFilename):
    outfile = file(indexFilename, "w")

    for emblFilename in emblFilenames:
        infile = file(emblFilename)

        timer.start("make index for " + emblFilename)
        while seekNextChrom(infile):
            index = infile.tell()
            chrom = readChrom(infile)
            print >>outfile, chrom.name, chrom.start, chrom.end, \
                             emblFilename, index
        timer.stop()



def makeGoFile(emblFilenames, index, goFilename):
    seeker = Seeker(emblFilenames, index)        
    out = file(goFilename, "w")
    
    for gene in index:
        infile = seeker.seek(gene)
        gos = readGos(infile)
        print >>out, gene, 
        for go in gos:
            print >>out, go,
        print >>out


def lookupChrom(geneName, lookup, geneIndex):
    # lookup chrom for gene
    (filename, pos) = geneIndex[geneName]
    lst = lookup[filename]
    i = 0
    for i in range(len(lst)):
        if lst[i][0] > pos:
            x = lst[i-1][1:]
            return (x[0], x[1]-1, x[2]-1)
    return (None, None, None)

    
def makeGenesFile(emblFilenames, chromIndexFile, geneIndex, genomeFilename):
    seeker = Seeker(emblFilenames, geneIndex)
    outfile = file(genomeFilename, "w")
    
    # read in chrom indexes
    lookup = readChromLookup(emblFilenames, filename)
    
    for geneName in geneIndex:
        infile = seeker.seek(geneName)
        
        (chrom, start, end) = lookupChrom(geneName, lookup, geneIndex)
        
        if chrom != None:
            gene = readGene(infile, start)
            if len(gene.coords) > 0:
                print >>outfile, gene.name, chrom, \
                                 gene.coords[0], gene.coords[1], gene.coords[2]
    
def makeExonsFile(emblFilenames, chromIndexFile, geneIndex, exonFilename):
    seeker = Seeker(emblFilenames, geneIndex)
    outfile = file(exonFilename, "w")
    
    # read in chrom indexes
    lookup = readChromLookup(emblFilenames, filename)
    
    for geneName in geneIndex:
        infile = seeker.seek(geneName)
        
        (chrom, start, end) = lookupChrom(geneName, lookup, geneIndex)
        
        if chrom != None:
            exonList = readExons(infile, start)
            if len(exonList) > 0:
                for exons in exonList:
                    print >>outfile, geneName, 
                    for i in exons:
                        print >>outfile, i,
                    print >>outfile

def makeSeqFiles(emblFilenames, chromsIndex, seqFilePrefix, conf):   
    index = readChromIndex(emblFilenames, chromsIndex)
    
    # make infiles
    infiles = {}
    for f in emblFilenames:
        infiles[f] = file(f)
    
    for chrom in index:
        if conf["exclude_NT"] and chrom.find("NT") != -1:
            continue
    
        out = file(seqFilePrefix + chrom + ".seq", "w")
        
        i = 1
        for entry in index[chrom]:
            if entry.start != i: 
                print "ERROR!"
                return
            
            # start reading
            infile = infiles[entry.filename]
            infile.seek(entry.pos)
            
            for line in util.SafeReadIter(infile):
                if line.startswith(seqPrefix):
                    out.write(line[5:71].replace(" ", ""))
                    break
            
            for line in util.SafeReadIter(infile):
                if not line.startswith(seqPrefix):
                    break
                else:
                    out.write(line[5:71].replace(" ", ""))
            
            i = entry.end + 1
        out.close()
    return

    # populate sequence files
    chrom = ""
    out = None
    
    for f in emblFilenames:
        infile = file(f)
        for line in infile:
            if line.startswith(chromPrefix):
                chrom = parseChrom(line).name
                if out != None:
                    out.close()
                    out = None
                if (not conf["exclude_NT"]) or (not chrom.startswith("NT")):
                    out = file(seqFilePrefix + chrom + ".seq", "a")
            if line.startswith(seqPrefix) and out != None:
                out.write(line[5:71].replace(" ", ""))


def processConfInit(conf):
    conf["make_gene_index"]  = True
    conf["make_chrom_index"] = True
    conf["make_genome"]      = True
    conf["make_exons"]       = True
    conf["make_go"]          = True
    conf["make_seq"]         = True
    conf["exclude_NT"]       = False
    

def processAll(species, files, prepDir, outputDir, conf):
    timer = util.Timer()

    # setup output files
    genesIndex  = prepDir + species + "_genes.index"
    chromsIndex = prepDir + species + "_chroms.index"
    goFile      = outputDir + species + ".go"
    genomeFile  = outputDir + species + ".genome"
    exonFile    = outputDir + species + ".exons"
    seqPrefix   = outputDir + species + "_"

    print genesIndex

    if conf["make_gene_index"]:
        makeGeneIndexFile(files, genesIndex)


    if conf["make_chrom_index"]:
        timer.start("make chrom index "+chromsIndex)
        makeChromIndexFile(files, chromsIndex)
        timer.stop()

    if conf["make_genome"]:
        timer.start("make genome file "+genomeFile)
        makeGenesFile(files, 
                      chromsIndex, 
                      readIndex(genesIndex), 
                      genomeFile)
        timer.stop()

    if conf["make_exons"]:
        timer.start("make exon file "+exonFile)
        makeExonsFile(files,
                      chromsIndex, 
                      readIndex(genesIndex), 
                      exonFile)
        timer.stop()

    if conf["make_go"]:
        timer.start("make go "+goFile)
        makeGoFile(files, readIndex(genesIndex), goFile)
        timer.stop()


    if conf["make_seq"]:
        timer.start("make sequence files")
        makeSeqFiles(files,
                     chromsIndex,
                     seqPrefix, conf)
        timer.stop()
