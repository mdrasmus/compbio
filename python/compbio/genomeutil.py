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
from rasmus import stats
from rasmus import treelib
from rasmus import util

from . import fasta
from . import gff, phylo


#
# Wed May  9 21:25:13 EDT 2007
# someday, I hope to delete all this code; its very DEPRECATED
#



#-------------------------------------------------------------------------
# Common options parsing
#-------------------------------------------------------------------------

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
    
    # read species map
    if "smap" in conf:
        conf["gene2species"] = readGene2species(*conf["smap"])
    else:
        conf["gene2species"] = gene2species
    
    if "stree" in conf:
        conf["stree"] = treelib.readTree(conf["stree"])




#============================================================================
# gene2species mappings and countings
#

gene2species = phylo.gene2species
makeGene2species = phylo.make_gene2species
readGene2species = phylo.read_gene2species


#----------------------------------------------------------------------------
# MultiBlocks
#----------------------------------------------------------------------------


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
    out = util.open_stream(filename, "w")
    
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
    
    infile = util.open_stream(filename, "r")
    
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

        

