from genomeutil import *

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


class MultiBlock:
    def __init__(self):
        self.segments = []
    
    def addSegments(self, segments):
        self.segments.append(segments)



def cutBlocks(blocks, refChrom, refStart, refEnd):
    segments = []
    
    for block in blocks:
        if block.chrom1 == refChrom:
            # ensure cut is bigger than ref
            assert refStart <= block.start1 and \
                   refEnd >= block.end1
            
            # add other block to segments
        else:
            # ensure cut is bigger than ref
            assert refStart <= block.start2 and \
                   refEnd >= block.end2
            
            # add other block to segments
            
    
    return segments
        



def buildForChrom(matching, refGenome, refChrom):
    chroms = {}
    
    refs = []
    
    # assign blocks to chroms of ref genome
    for block in matching.blocks:
        if block.start1 > block.end1 or \
           block.start2 > block.end2:
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
    blocks = {}
    multiBlocks = []
    
    
    for ref in refs:
        #print " ", ref.kind, ref.pos, ref.otherGenome().name, ref.otherChrom().name
        multiBlocks.append(MultiBlock(refChrom, lastPos, ref.pos, segments))

        lastPos = ref.pos

        # add or remove from current blocks
        if ref.kind == "start":
            blocks[ref.block] = True
        else:
            del blocks[ref.block]
    multiBlocks.append(MultiBlock(refChrom, lastPos, refChrom.size, blocks))

    return multiBlocks


def buildForGenome(matching, refGenome):
    for chrom in refGenome.chroms.values():
        chrom.multiBlocks = buildForChrom(matching, refGenome, chrom)

        
        
