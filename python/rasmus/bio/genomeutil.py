#
# genomeutil.py
# Genome and gene objects
#

# python libs
import copy
import math
import os
import sys

# rasmus libs
from rasmus import algorithms
from rasmus import env
from rasmus import stats
from rasmus import treelib
from rasmus import util

from rasmus.bio import fasta
from rasmus.bio import gff


#
# Wed May  9 21:25:13 EDT 2007
# someday, I hope to delete all this code; its very DEPRECATED
#



#--------------------------------------------------------------------------------
# Common options parsing
#--------------------------------------------------------------------------------

options = [
  ["s:", "stree=", "stree", "<species tree>",
    {"single": True,
     "help": "species tree (newick format)"}],
  ["S:", "smap=", "smap", "<gene2species map>",
    {"help": "mapping of gene names to species names"}],
  ["P:", "paths=", "paths", "<files path>",
    {"help": "colon separated paths used to search for data files",
     "default": "."}],
]

def readOptions(conf):
    """Setup data paths and parse common options"""

    # setup data paths
    env.addPaths(* conf["paths"])
    env.addEnvPaths("DATAPATH")
    
    # read species map
    if "smap" in conf:
        conf["gene2species"] = readGene2species(*map(env.findFile, conf["smap"]))
    else:
        conf["gene2species"] = gene2species
    
    if "stree" in conf:
        conf["stree"] = treelib.readTree(env.findFile(conf["stree"]))




#============================================================================
# gene2species mappings and countings
#
# TODO: I am starting to move these to phylo since I use them more 
# over there.  However, I will not necessarily always use gene2species for
# phylo... still thinking it over...
#

def gene2species(genename):
    """default gene2species mapping"""
    return genename


def makeGene2species(maps):
    # find exact matches and expressions
    exacts = {}
    exps = []
    for mapping in maps:
        if "*" not in mapping[0]:
            exacts[mapping[0]] = mapping[1]
        else:
            exps.append(mapping)
    
    # create mapping function
    def gene2species(gene):
        # eval expressions first in order of appearance
        for exp, species in exps:
            if exp[-1] == "*":
                if gene.startswith(exp[:-1]):
                    return species
            elif exp[0] == "*":
                if gene.endswith(exp[1:]):
                    return species
        
        if gene in exacts:
            return exacts[gene]
        
        raise Exception("Cannot map gene '%s' to any species" % gene)
    return gene2species


def readGene2species(* filenames):
    for filename in filenames:
        if filename == "DEFAULT":
            smap = gene2species
        else:
            maps = []
            for filename in filenames:
                maps.extend(util.readDelim(util.skipComments(
                            util.openStream(filename))))
            smap = makeGene2species(maps)
    
    return smap


# this does not need to be in a module
"""
def genomeComposition(genomes, comp, gene2species=gene2species):
    counts = {}
    for genome in genomes:
        counts[genome] = 0
    for gene in comp:
        genome = gene2species(gene)
        if genome in genomes:
            counts[genome] += 1
    return counts


def componentCompositions(order, comps, gene2species=gene2species):
    compositions = util.Dict(1, 0)
    for comp in comps:
        counts = genomeComposition(order, comp, gene2species)
        key = []
        for genome in order:
            key.append(counts[genome])
        compositions[tuple(key)] += 1
    return compositions.data
"""



# constants
aminoAcids = "FLSYCWLPHQRIMTNKVADEG"
aaPseudoCounts = {"-": 1/21.0}

for i in aminoAcids:
    aaPseudoCounts[i] = 1/21.0

        
#============================================================================
# Gene and genome objects for synteny and blast matches
#


class Gene:
    def __init__(self):
        self.name  = ""
        self.chrom = None
        self.start = 0
        self.end   = 0
        self.direction = 1
        self.description = ""
        self.matches = []
        self.selfMatches = []
        self.pruned = []
        self.trans = {}
        self.order = 0
    
    
    def __repr__(self):
        if self.chrom != None:
            return "<GENE " + self.chrom.genome.name + ":" + \
                         self.chrom.name + ":" + self.name + ">"
        else:
            return "<GENE " + self.name + ">"
    
    def length(self):
        return self.end - self.start + 1
    
    # used by Matching
    def allMatches(self):
        return self.matches + self.selfMatches
    
    def neighbors(self):
        return map(lambda x: x.otherGene(self), self.allMatches())
    
    def synteny(self):
        return filter(lambda x: x.block != None, self.matches)
    
    """
    def prune(self, match):
        if match.isSelfMatch():
            self.selfMatches.remove(match)
        else:
            self.matches.remove(match)
        self.pruned.append(match)

    def unprune(self, match):
        if match.isSelfMatch():
            self.selfMatches.append(match)
        else:
            self.matches.append(match)
        self.pruned.remove(match)
    """
    
    def attach(self, match):
        if match.isSelfMatch():
            self.selfMatches.append(match)
        else:
            self.matches.append(match)

    def sequence(self):
        seqdb = self.chrom.genome.seqdb
        if seqdb:
            return seqdb.geneTranscript(self)
        else:
            return ""

    def protein(self):
        seqdb = self.chrom.genome.seqdb
        if seqdb:
            return seqdb.geneProtein(self)
        else:
            return ""
        



class Transcript:
    def __init__(self, gene = None, name = ""):
        self.exons = []
        self.gene = gene
        self.name = name
        

class Exon:   
    def __init__(self, trans, start, end, direction, codeStart, codeEnd):
        self.name = ""
        self.trans = trans
        self.start = start
        self.end = end
        self.direction = direction
        self.codeStart = codeStart
        self.codeEnd = codeEnd

    
    def __repr__(self):
        if self.trans:
            return "<EXON " + self.trans.gene.chrom.genome.name + ":" + \
                   self.trans.gene.chrom.name + ":" + str(self.start) + "," + \
                   str(self.end) + "," + str(self.direction) + ">"

    def sequence(self):
        return self.trans.gene.chrom.genome.seqdb.exonSeq(self)

      

class Chromosome:
    def __init__(self, name = "", genome = None):
        self.name = name
        self.genes = []
        self.genome = genome
        self.size = 0
    
    def __repr__(self):
        if self.genome:
            return "<CHROM " + self.genome.name + ":" + self.name + ">"



class Genome:
    def __init__(self, name = ""):
        self.name = name
        self.size = 0
        self.genes = {}
        self.chroms = {}

    def read(self, filename):
        infile = file(filename)
    
        for line in infile:
            tokens = line.rstrip().split()
            
            # add chromosome if needed
            chromName = tokens[1]
            if not self.chroms.has_key(chromName):
                chrom = Chromosome(chromName)
                chrom.genome = self
                self.chroms[chromName] = chrom
            
            # create gene
            gene = Gene()
            gene.name  = tokens[0]
            gene.chrom = self.chroms[chromName]
            gene.start = int(tokens[2])
            gene.end   = int(tokens[3])
            if tokens[4] == "+" or tokens[4] == "1":
                gene.direction = 1
            else:
                gene.direction = -1
            self.genes[gene.name] = gene
    
    
    def autoconf(self, chroms=None):
        if chroms == None:
            chroms = self.defaultOrder()
    
        # update chromosomes
        for chrom in self.chroms.itervalues():
            chrom.genes = []
        
        # set sizes
        for gene in self.genes.itervalues():
            gene.chrom.genes.append(gene)
            if gene.end > gene.chrom.size:
                gene.chrom.size = gene.end
        self.size = sum(map(lambda x: x.size, self.chroms.values()))
        
        # sort genes in chroms
        for chrom in self.chroms.itervalues():
            chrom.genes.sort(lambda a, b: cmp(a.start, b.start))
            for i in xrange(len(chrom.genes)):
                chrom.genes[i].order = i
        
        # set chromosome order
        for i in xrange(len(chroms)):
            self.chroms[chroms[i]].order = i
    
    
    def getChromOrder(self):
        chroms = self.chroms.keys()
        chroms.sort(lambda a,b: cmp(self.chroms[a].order,
                                    self.chroms[b].order))
        return chroms
    
    
    def autoconf2(self, chrOrder=None):
        # update chromosomes
        for chrom in self.chroms.itervalues():
            chrom.genes = []
        
        for gene in self.genes.itervalues():
            gene.chrom.genes.append(gene)
            if gene.end > gene.chrom.size:
                gene.chrom.size = gene.end
        
        # sort genes in chroms
        for chrom in self.chroms.itervalues():
            chrom.genes.sort(lambda a, b: cmp(a.start, b.start))
            for i in xrange(len(chrom.genes)):
                chrom.genes[i].order = i
        
        # set order and position
        if chrOrder == None:
            chrOrder = self.defaultOrder()
        self.setOrder(chrOrder)
    
            
    # make a default chromosome order
    # sort by size
    def defaultOrder(self):
        keys = self.chroms.keys()
        keys.sort(lambda a,b: cmp(self.chroms[b].size, self.chroms[a].size))
        return keys

        
        
        

class GeneIter:
    """An iterator that walks down a sorted list of genes"""

    def __init__(self, genes, start, end, index=None):
        self.genes  = genes
        self.start  = start
        self.end    = end
        self.ngenes = len(genes)
        
        if index != None:
            self.index = index
        else:
            self.index = 0
            
            # find starting index by binary search
            low, top = algorithms.binsearch(genes, start-1, 
                                            lambda a,b: cmp(a.start, b))
            
            if top != None:
                self.index = top
            else:
                self.index = self.ngenes
    
    def __iter__(self):
        return self
    
    
    def next(self):
        if (self.index < self.ngenes) and \
           (self.genes[self.index].start < self.end):
            gene = self.genes[self.index]
        else:
            raise StopIteration
        
        self.index += 1
        return gene
        



class Segment:
    def __init__(self, genome, chrom, start, end, direction):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.end = end
        self.direction = direction
    
    
    def __str__(self):
        if isinstance(self.genome, str):
            return "%s\t%s\t%d\t%d\t%d" % \
                   (self.genome, self.chrom, self.start, self.end, self.direction)
        else:
            return "%s\t%s\t%d\t%d\t%d" % \
                   (self.genome.name, self.chrom.name, self.start, self.end, self.direction)

    def __repr__(self):
        return self.__str__()

class MultiBlock:
    def __init__(self, segments=None):
        self.segments = []
        if segments != None:
            self.segments.extend(segments)
            
    
    def addSegments(self, segments):
        self.segments.extend(segments)   
    
    
    def __str__(self):
        text = []
        for seg in self.segments:
            text.append(str(seg) + "\n")
        return "".join(text)

    def __repr__(self):
        return self.__str__()


class SyntenyBlock:
    blockid = -1
    genome1 = None
    chrom1  = None
    genome2 = None
    chrom2  = None
    start1  = 0
    start2  = 0
    end1    = 0
    end2    = 0
    direction = None
    
    
    """
    widths = None
    slope = None
    """
    
    genes = None
    
    def __init__(self):
        self.matches = []
        #self.genes = {}
    
    
    def init(self, match):
        (gene1, gene2) = match.genes
        
        self.genome1 = gene1.chrom.genome
        self.chrom1  = gene1.chrom
        self.genome2 = gene2.chrom.genome
        self.chrom2  = gene2.chrom
        self.start1  = gene1.start
        self.start2  = gene2.start
        self.end1    = gene1.end
        self.end2    = gene2.end
        self.matches = [match]
        
        match.block = self
    
    def merge(self, block):
        self.start1 = min(self.start1, block.start1)
        self.start2 = min(self.start2, block.start2)
        self.end1 = max(self.end1, block.end1)
        self.end2 = max(self.end2, block.end2)
        for match in block.matches:
            self.matches.append(match)
        block.matches = []
    
    def add(self, match, init = True):
        match.block = self
        if init:
            if len(self.matches) == 0:
                self.init(match)
            else:
                self.start1 = min(self.start1, match.genes[0].start)
                self.start2 = min(self.start2, match.genes[1].start)
                self.end1 = max(self.end1, match.genes[0].end)
                self.end2 = max(self.end2, match.genes[1].end)
        self.matches.append(match)
    
    def otherGenome(self, genome):
        if self.genome1 == genome:
            return self.genome2
        else:
            return self.genome1
    
    def otherChrom(self, chrom):
        if self.chrom1 == chrom:
            return self.chrom2
        else:
            return self.chrom1
    
    
    
    def hasGene(self, gene):
        for match in self.matches:
            if gene in match.genes:
                return True
        return False
    
    def getGenes(self):
        genes = {}
        for match in self.matches:
            for gene in match.genes:
                genes[gene] = True
        return genes
    
    def length1(self):
        return self.end1 - self.start1 + 1
    
    def length2(self):
        return self.end2 - self.start2 + 1
    
    

    def getDirection(self):
        if self.matches[0].direction():
            return 1
        else:
            return -1
        
        #vote = 0
        #for match in self.matches:
        #    if match.direction():
        #        vote += 1
        #    else:
        #        vote -= 1
        #if vote > 0:
        #    self.direction = 1
        #else:
        #    self.direction = -1
        #return self.direction
    
    def forget(self):
        self.widths = None
        self.slope = None
        self.genes = None




class Match:
    genes = []
    score = None
    block = None
    #pruned = False
    #resurrected = False
    
    def __init__(self, gene1 = None, gene2 = None, score = None, block = None):
        self.genes = [gene1, gene2]
        self.score = score
        self.block = block

    def otherGene(self, gene):
        if self.genes[0] != gene:
            return self.genes[0]
        else:
            return self.genes[1]
    
    def direction(self):
        return self.genes[0].direction == self.genes[1].direction
    
    """
    def prune(self):
        self.pruned = True
        for gene in self.genes:
            gene.prune(self)

    def unprune(self):
        self.pruned = False
        for gene in self.genes:
            gene.unprune(self)
        self.resurrected = True
    """
    
    def attach(self):
        for gene in self.genes:
            gene.attach(self)
    
    def isSelfMatch(self):
        return self.genes[0].chrom.genome == self.genes[1].chrom.genome
            
"""
class HomologyGroup:
    groupid = 0
    genes = []
    
    def __init__(self, groupid):
        self.groupid = groupid
        self.genes = []
"""


class Matching:   
    def __init__(self):
        self.genomes = {}
        self.matches = []
        self.blocks = []
        self.genes = {}
        self.comps = []
        self.complookup = []
        
        
        """
        self.homology = []
        """
    
    def readGenomes(self, filename, gene2species=gene2species):
        if filename.endswith(".coord"):
            self.readCoordFile(filename, gene2species)
        elif filename.endswith(".gff"):
            self.readGffFile(filename, gene2species)
    
    
    def readGffFile(self, filename, gene2species, feature="gene"):
        """assume one transcript per gene"""
        
        for region in gff.iterGff(filename,
                            regionFilter=lambda x: x.feature == feature):
            
            geneName = region.data["ID"]

            
            # don't add duplicate genes
            if geneName not in self.genes:
                genomeName = gene2species(geneName)
            
                self.addGene(geneName, genomeName, region.seqname, 
                             region.start, region.end, region.strand)
    
    def addRegions(self, regions, gene2species):
        for region in regions:
            if "ID" in region.data:
                geneName = region.data["ID"]
                genomeName = gene2species(geneName)
                self.addGene(geneName, genomeName, region.seqname, 
                             region.start, region.end, region.strand)
    
    
    def readCoordFile(self, filename, gene2species):
        infile = util.openStream(filename)
        for line in infile:
            geneName, chromName, start, end, strand = line.rstrip().split("\t")
            genomeName = gene2species(geneName)
            
            start = int(start)
            end = int(end)
            
            if strand == "+" or strand == "1":
                strand = 1
            else:
                strand = -1
            
            self.addGene(geneName, genomeName, chromName, start, end, strand)
            
    
    def addGene(self, geneName, genomeName, chromName, start, end, strand):
        # add genome if needed
        if genomeName not in self.genomes:
            genome = Genome(genomeName)
            self.genomes[genomeName] = genome
        else:
            genome = self.genomes[genomeName]

        # add chromosome if needed
        if chromName not in genome.chroms:
            chrom = Chromosome(chromName, genome)
            genome.chroms[chromName] = chrom

        # create gene
        gene = Gene()
        gene.name  = geneName
        gene.chrom = genome.chroms[chromName]
        gene.start = start
        gene.end   = end
        gene.direction = strand

        genome.genes[geneName] = gene
        self.genes[geneName] = gene
    
    
    def autoconf(self, genomes = None):
        if genomes == None:
            genomes = self.genomes.keys()
            genomes.sort()
    
        # configure genomes
        for i in xrange(len(genomes)):
            self.genomes[genomes[i]].autoconf()
            self.genomes[genomes[i]].order = i
    
    
    def getGenomeOrder(self):
        genomes = self.genomes.keys()
        genomes.sort(lambda a,b: cmp(self.genomes[a].order,
                                     self.genomes[b].order))
        return genomes
    
    
    def getGenes(self):
        # TODO: make sure caller does not change dict
        return self.genes
    
    
    
    def readMatches(self, infile, cols=[0,1,2], minscore=0):
        genes = self.getGenes()
        
        for line in infile:
            tokens = line.split()
            geneName1 = tokens[cols[0]]
            geneName2 = tokens[cols[1]]
            score     = float(tokens[cols[2]])
            
            # skip matches with unknown genes
            if (geneName1 == geneName2) or \
               (score < minscore) or \
               (geneName1 not in genes) or \
               (geneName2 not in genes):
                continue
            
            self.addMatch(genes[geneName1], genes[geneName2], score)
    
    def writeSynteny(self, filename):
        out = file(filename, "w")
        for match in self.matches:
            if match.block != None:
                print >>out, match.block.blockid
            else:
                print >>out, "-1"
                    
    def writeGeneMatches(self, outfile):
        out = file(outfile, "w")
        
        for match in self.matches:
            if match.pruned:
                continue
            
            print >>out, "%s\t%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d" % ( 
                match.genes[0].chrom.genome.name, match.genes[0].chrom.name,
                match.genes[0].start, match.genes[0].end, 
                match.genes[0].direction,
                match.genes[1].chrom.genome.name, match.genes[1].chrom.name,
                match.genes[1].start, match.genes[1].end, 
                match.genes[1].direction)
   
   
    def clearMatches(self):
        """clear all matches from Matching and genes"""
        
        self.matches = []
        self.blocks = []
        genes = self.getGenes()
        
        for gene in genes.values():
            gene.matches = []
        
    
    def addMatch(self, gene1, gene2, score):
        match = Match(gene1, gene2, score)
        self.matches.append(match)
        match.attach()
    
    
    def setGeneComponents(self, comps):
        self.comps = []
        for comp in comps:
            self.comps.append(map(lambda x: self.genes[x], comp))
        
        self.complookup = {}
        for i in xrange(len(self.comps)):
            for gene in self.comps[i]:
                self.complookup[gene] = i
    
    def getComp(self, gene):
        if gene in self.complookup:
            return self.comps[self.complookup[gene]]
        else:
            return [gene]
    

    def interpolateGene(self, gene1, chrom2):
        if gene1 not in self.complookup:
            return []
        else:
            comp = self.comps[self.complookup[gene1]]
            genes2 = filter(lambda gene2: gene2.chrom == chrom2, comp)
            return genes2

    def findNearGenes(self, chrom, pos):
        genes = chrom.genes
        low, top = algorithms.binsearch(genes, pos, lambda a, b: cmp(a.start, b))
        if low != None:
            gene1 = genes[low]
        else:
            gene1 = None
        if top != None:
            gene2 = genes[top]
        else:
            gene2 = None
        return gene1, gene2
    
    
    def interpolateCoord(self, chrom1, pos, chrom2, direction):
        # find neighboring genes
        before, after = self.findNearGenes(chrom1, pos)
        #print "query", pos, before, after                
        
        # determine gene to start walking from
        if direction == 1:
            if after != None:
                i = after.order
            else:
                i = before.order
        else:
            if before != None:
                i = before.order
            else:
                i = after.order            
        
        # walk until a match is found between chrom1 and chrom2
        while 0 <= i < len(chrom1.genes):
            genes2 = self.interpolateGene(chrom1.genes[i], chrom2)
            #print "hit", genes2
            if len(genes2) > 0:
                return genes2
            i += direction
        
        return []
    

    def interpolateSegment(self, chrom1, start1, end1, chrom2):
        genes2 = self.interpolateCoord(chrom1, start1, chrom2, 1)
        genes3 = self.interpolateCoord(chrom1, end1, chrom2, -1)
        genes4 = genes2 + genes3
        
        if len(genes2) == 0 or len(genes3) == 0:
            return None, None, None
        
        if genes2[0].start < genes3[0].start:
            direction = 1
        else:
            direction = -1

        start2 = min(map(lambda x: x.start, genes4))
        end2 = max(map(lambda x: x.end, genes4))

        return start2, end2, direction
        



#--------------------------------------------------------------------------------
# MultiBlocks
#--------------------------------------------------------------------------------


#
# TODO: multiblocks can be done with my Regions in regionlib!
#


def cutBlocks(matching, blocks, refChrom, refStart, refEnd):
    segments = [Segment(refChrom.genome, refChrom, refStart, refEnd, 1)]
    
    for block in blocks:
        if block.chrom1 == refChrom:
            # ensure cut is bigger than ref
            #assert refStart <= block.start1 and \
            #       refEnd >= block.end1
            otherGenome = block.genome2
            otherChrom = block.chrom2
        else:
            # ensure cut is bigger than ref
            #assert refStart <= block.start2 and \
            #       refEnd >= block.end2
            otherGenome = block.genome1
            otherChrom = block.chrom1
    
        # add other block to segments
        start, end, direction = \
            matching.interpolateSegment(refChrom, refStart, refEnd, otherChrom)
        if start != None:
            segments.append(Segment(otherGenome, otherChrom, start, end, direction))
    
    
    return MultiBlock(segments)



class BlockRef:   
    def __init__(self, kind, block, side):
        self.kind = kind    
        self.block = block    
        self.side = side
        if side == 1:
            if kind == "start":
                self.pos = block.start1
            else:
                self.pos = block.end1
        else:
            if kind == "start":
                self.pos = block.start2
            else:
                self.pos = block.end2
    
    def otherGenome(self):
        if self.side == 1:
            return self.block.genome2
        else:
            return self.block.genome1

    def otherChrom(self):
        if self.side == 1:
            return self.block.chrom2
        else:
            return self.block.chrom1


def makeChromMultiBlocks(conf, matching, refGenome, refChrom, blocks=None):
    chroms = {}
    
    refs = []
    
    if blocks == None:
        blocks = matching.blocks
    
    # assign blocks to chroms of ref genome
    for block in blocks:
        if block.start1 > block.end1 or \
           block.start2 > block.end2:
            # maybe raise error
            continue
    
        if block.genome1 == refGenome and block.chrom1 == refChrom:
            refs.append(BlockRef("start", block, 1))
            refs.append(BlockRef("end", block, 1))
        elif block.genome2 == refGenome and block.chrom2 == refChrom:
            refs.append(BlockRef("start", block, 2))
            refs.append(BlockRef("end", block, 2))            
    
    # sort reference in order
    refs.sort(lambda a,b: cmp(a.pos, b.pos))
    
    # produce multiblocks
    lastPos = 0
    blocks = []
    multiblocks = []
    
    for ref in refs:
        # add or remove from current blocks
        if ref.kind == "start":
            if ref.pos > lastPos:
                multiblocks.append(cutBlocks(matching, blocks, refChrom, 
                                             lastPos, ref.pos-1))
            lastPos = ref.pos
            blocks.append(ref.block)
        else:
            if ref.pos > lastPos:
                multiblocks.append(cutBlocks(matching, blocks, refChrom, 
                                             lastPos, ref.pos))
                lastPos = ref.pos+1
            blocks.remove(ref.block)
        
        

    if refChrom.size > lastPos:
        multiblocks.append(cutBlocks(matching, blocks, refChrom, 
                                     lastPos, refChrom.size))

    return multiblocks


def makeGenomeMultiBlocks(conf, matching, refGenome):
    # split blocks by chromosome
    blocks = util.Dict(1, [])
    for block in matching.blocks:
        if block.genome1 == refGenome:
            blocks[block.chrom1].append(block)
        elif block.genome2 == refGenome:
            blocks[block.chrom2].append(block)
    
    # find multiblocks for each chromosome
    multiblocks = []
    for chrom in refGenome.chroms.values():
        if chrom in blocks:
            multiblocks.extend(
                makeChromMultiBlocks(conf, matching, refGenome, 
                                     chrom, blocks[chrom]))
    return multiblocks


def writeMultiBlocks(filename, multiblocks):
    out = util.openStream(filename, "w")
    
    for multiblock in multiblocks:
        if len(multiblock.segments) > 0:
            out.write("\t".join([multiblock.segments[0].genome.name, 
                                 multiblock.segments[0].chrom.name, 
                                 str(multiblock.segments[0].start), 
                                 str(multiblock.segments[0].end),
                                 str(multiblock.segments[0].direction)]))
            
        for segment in multiblock.segments[1:]:
            out.write("\t")
            out.write("\t".join([segment.genome.name, segment.chrom.name, 
                                 str(segment.start), str(segment.end),
                                 str(segment.direction)]))
        out.write("\n")
  

def iterMultiBlocks(filename):
    """
    iterates over multiblocks in a file
    
    Warning: this keeps genomes and chroms as names (strings) NOT objects
    I think I will eventually convert everything to names.
    But currently this read is incompatiable with writeMultiBlocks
    """
    
    infile = util.openStream(filename, "r")
    
    for line in infile:
        tokens = line.split("\t")
        
        segments = []
        for i in range(0, len(tokens), 5):
            genome, chrom, start, end, direction = tokens[i:i+5]
            segments.append(Segment(genome, chrom, 
                                    int(start), int(end), int(direction)))
        
        yield MultiBlock(segments)
    


def genesInMultiBlocks(genes, blocks):
    """
    Returns which (if any) multiblock a set of genes are in.
    
    genes must hit unique segments in multiblock.
    
    MultiBlocks must refer to genomes and chroms as strings.
    
    """

    for block in blocks:
        found = set()
        for gene in genes:
            for seg in block.segments:
                if gene.chrom.genome.name == seg.genome and \
                   gene.chrom.name == seg.chrom and \
                   gene.start >= seg.start and \
                   gene.end <= seg.end:
                    found.add(seg)
                    break
            else:
                # cannot find segment for gene, no match
                break

        # did we find a unique segment for each gene?
        if len(found) == len(genes):
            return block
    return None

        
