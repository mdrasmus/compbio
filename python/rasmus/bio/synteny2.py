# python libs
import copy
import math
import os
import sys

# rasmus libs
from rasmus import algorithms
from rasmus import env
from rasmus import graph
from rasmus import stats
from rasmus import util
from rasmus import treelib


# rasmus.bio libs
from rasmus.bio.genomeutil import *
from rasmus.bio import bionj
from rasmus.bio import genomeio
from rasmus.bio import muscle
from rasmus.bio import phylip
from rasmus.bio import phylo



def initConf(conf = None):
    if conf == None:
        conf = {}
    conf['blockExtend'] = 200000
    conf['blockCloseMerge'] = 400000
    conf['blockInterupt'] = 3
    
    conf['minscore'] = 60
    conf['relScore'] = .8
    conf['rangeGenes'] = 20
    conf['range1'] = None
    conf['range2'] = None
    conf['sigScore'] = 500
    conf['sigBbhScore'] = 500
    conf['minBlockSize'] = 2
    conf["nearRange"] = 5
    conf["pullFrac"] =  .3
    conf["minPullLen"] = .3
    conf["minGeneSize"] = 100
    conf["synteny_merge_only"] = False
    conf["synteny_resume"] = False
    
    return conf
    


def separateBlocks(blocks, method="bychrom"):
    blockLists = {}
    
    if method == "bychrom":
        for block in blocks:
            lst = blockLists.setdefault((block.chrom1, block.chrom2), [])
            lst.append(block)
    else:
        for block in blocks:
            lst = blockLists.setdefault((block.genome1, block.genome2), [])
            lst.append(block)
    
    return blockLists
    

def bestBidir(conf, genomes, matches):
    "find best bidirectional matches"
    
    bestMatches = []
    best = {}
    
    for match in matches:
        best[match] = 0
    
    # mark best matches
    for genome in genomes:
        for gene in genome.genes.values():
            if gene.length() > conf["minGeneSize"]:
                if len(gene.matches) > 0:
                    besti = util.argmax(gene.matches,
                                        key=lambda x: x.score)
                    best[gene.matches[besti]] += 1
    
    for match in matches:
        if best[match] == 2:
            bestMatches.append(match)
    return bestMatches


def bestUnidir(genomes, matches):
    zeroMatch = Match()
    zeroMatch.score = 0
    best = util.Dict(1, zeroMatch)
    for match in matches:
        if match.score > best[match.genes[0]].score:
            best[match.genes[0]] = match
        if match.score > best[match.genes[1]].score:
            best[match.genes[1]] = match
    
    return best.data

        

def blockDistance(conf, direction, block1, block2, side):
    if block1 == block2:
        return 1e1000
    
    b1 = block1
    
    # swap if finding backwards distance
    if side == -1:
        tmp = block1
        block1 = block2
        block2 = tmp
        
    if direction[b1] == 1:
        dist1 = block2.start1 - block1.end1
        dist2 = block2.start2 - block1.end2
    else:
        dist1 = block2.start1 - block1.end1
        dist2 = block1.start2 - block2.end2
    
    # blocks that are "behind" you are infinitly far away
    if dist1 < -conf['blockCloseMerge'] or \
       dist2 < -conf['blockCloseMerge'] or \
       (dist1 < 0  and dist2 < 0):
        return 1e1000
    else:
        return abs(dist1) + abs(dist2)



def getNearBlocks(conf, direction, block1, blocks, side, quadtrees):
    "optimization to speed up finding near blocks"
    
    radius = conf['blockExtend']
    step = conf['blockExtend']
    #maxRadiurs = max(block1.chrom1.size, block1.chrom2.size)
    radius1 = conf['range1']
    radius2 = conf['range2']
    nearBlocks = []
    
    # setup search box
    if direction[block1] == 1:
        if side == 1:
            x1 = block1.end1
            y1 = block1.end2
        else:
            x1 = block1.start1 - radius1
            y1 = block1.start2 - radius2
    else:
        if side == 1:        
            x1 = block1.end1
            y1 = block1.start2 - radius2
        else:
            x1 = block1.start1 - radius1
            y1 = block1.end2
    x2 = x1 + radius1
    y2 = y1 + radius2        
    
    nearBlocks = quadtrees[block1.chrom1][block1.chrom2].query(
                    algorithms.Rect(x1, y1, x2, y2))

    return nearBlocks


def findNearestBlock(conf, direction, block1, blocks, side, quadtrees):
    nearBlocks = getNearBlocks(conf, direction, block1, blocks, side, quadtrees)
    
    if len(nearBlocks) == 0:
        return None
    
    # search for nearest block
    i = util.argmin(nearBlocks, key=lambda x: 
                    blockDistance(conf, direction, block1, x, side))
    block2 = nearBlocks[i]

    if blockDistance(conf, direction, block1, block2, side) < 1e1000:
        # only connect if nearest block is connectable
        if direction[block1] == direction[block2]:
            #util.note("CONNECT %s %s" % (block1.matches[0].genes[0], 
            #                            (block2.matches[0].genes[0])))
            return block2
    #    else:
    #    
    #        #util.note("BLOCK %s %s" % (block1.matches[0].genes[0], 
    #                                  (block2.matches[0].genes[0])))
    #else:
    #    util.note("NOCONNECT %s %s" % (block1.matches[0].genes[0], 
    #                                  (block2.matches[0].genes[0])))
    return None


def findConnects(conf, blocks, quadtrees):
    blockLists = separateBlocks(blocks)
    connects = {}
    direction = {}
    
    # find all block directions in advance
    for block in blocks:
        direction[block] = block.getDirection()
    
    # find closest manhattan distance blocks with same direction
    for lst in blockLists.values():
        for block1 in lst:
            # find nearest block in the forward direction
            nearest1 = findNearestBlock(conf, direction, block1, \
                                        lst, 1, quadtrees)
            
            # find nearest block in the backward direction
            nearest2 = findNearestBlock(conf, direction, block1, \
                                        lst, -1, quadtrees)
            
            connects[block1] = []
            if nearest1 != None:
                connects[block1].append(nearest1)
            if nearest2 != None:
                connects[block1].append(nearest2)
    
    # now only keep connections that both blocks agree are best
    # (keep only bidirectional connections)
    connects2 = []
    for block1 in connects:
        for block2 in connects[block1]:
            if block1 in connects[block2]:
                connects2.append((block1, block2))
    
    return connects2


def mergeBlocks(conf, blocks, connects):
    # find resultant blocks by connected components
    import graph
    
    connections = {}
    for block in blocks:
        connections[block] = []
    for connect in connects:
        connections[connect[0]].append(connect[1])
        connections[connect[1]].append(connect[0])        
    
    def getNeighbor(block):
        return connections[block]
    
    components = graph.connectedComponents(blocks, getNeighbor)

    util.tic("merge connected blocks")
    dels = {}
    for comp in components:
        root = comp[0]
        for other in comp[1:]:
            root.merge(other)
            dels[other] = 1
    blocks[:] = filter(lambda x: x not in dels, blocks)
    util.toc()
    

def connectBlocks(conf, blocks, quadtrees):
    connects = findConnects(conf, blocks, quadtrees)
    mergeBlocks(conf, blocks, connects)
    ensureColinarity(conf, blocks)
    
    
def mergeCloseBlocks(conf, blocks, quadtrees):
    direction = {}
    connects = []
    
    # find all block directions in advance
    for block in blocks:
        direction[block] = block.getDirection()
        block.matches.sort(lambda a,b: cmp(a.genes[0].start, b.genes[0].start))

    # find all near blocks
    for block in blocks:
        nearBlocks = []
        if direction[block] == 1:
            nearBlocks.extend(quadtrees[block.chrom1][block.chrom2].query(
                        algorithms.Rect(
                            block.end1 - conf["blockCloseMerge1"],
                            block.end2 - conf["blockCloseMerge2"],
                            block.end1 + conf["blockCloseMerge1"],
                            block.end2 + conf["blockCloseMerge2"])))
            nearBlocks.extend(quadtrees[block.chrom1][block.chrom2].query(
                        algorithms.Rect(
                            block.start1 - conf["blockCloseMerge1"],
                            block.start2 - conf["blockCloseMerge2"],
                            block.start1 + conf["blockCloseMerge1"],
                            block.start2 + conf["blockCloseMerge2"])))
        else:
            nearBlocks.extend(quadtrees[block.chrom1][block.chrom2].query(
                        algorithms.Rect(
                            block.start1 - conf["blockCloseMerge1"],
                            block.end2 - conf["blockCloseMerge2"],
                            block.start1 + conf["blockCloseMerge1"],
                            block.end2 + conf["blockCloseMerge2"])))
            nearBlocks.extend(quadtrees[block.chrom1][block.chrom2].query(
                        algorithms.Rect(
                            block.end1 - conf["blockCloseMerge1"],
                            block.start2 - conf["blockCloseMerge2"],
                            block.end1 + conf["blockCloseMerge1"],
                            block.start2 + conf["blockCloseMerge2"])))
        
        # only blocks in the same direction
        for near in nearBlocks:
            if direction[near] == direction[block] and \
               (util.overlap(block.matches[0].genes[0].start, 
                             block.matches[0].genes[0].end, 
                             near.matches[-1].genes[0].start, 
                             near.matches[-1].genes[0].end) or
                util.overlap(block.matches[0].genes[1].start, 
                             block.matches[0].genes[1].end, 
                             near.matches[-1].genes[1].start, 
                             near.matches[-1].genes[1].end) or
                util.overlap(block.matches[-1].genes[0].start, 
                             block.matches[-1].genes[0].end, 
                             near.matches[0].genes[0].start, 
                             near.matches[0].genes[0].end) or
                util.overlap(block.matches[-1].genes[1].start, 
                             block.matches[-1].genes[1].end, 
                             near.matches[0].genes[1].start, 
                             near.matches[0].genes[1].end)):
                connects.append([block, near])
    mergeBlocks(conf, blocks, connects)
    
    ensureColinarity(conf, blocks)
    


def ensureColinarity(conf, blocks):
    for block in blocks:
        badmatches = []
        bdir = block.getDirection()
        block.matches.sort(lambda a,b: matchSort(a,b,bdir))
        
        lastx = 0
        if bdir == 1:
            lasty = 0
        else:
            lasty = 1e1000
        
        for match in block.matches:           
            if matchDirection(match) != bdir or \
               match.genes[0].start < lastx or \
               (bdir == 1 and match.genes[1].start < lasty) or \
               (bdir == -1 and match.genes[1].start > lasty):
                badmatches.append(match)
            else:
                lastx = match.genes[0].start
                lasty = match.genes[1].start
        
        if len(badmatches) > 0:
            util.log("found %d badmatches" % len(badmatches))
        block.matches = util.remove(block.matches, *badmatches)



def makeBlockLookup(blocks):
    # build match to block lookup
    lookup = {}
    for block in blocks:
        for match in block.matches:
            lookup[match] = block
    return lookup


def matchDirection(match):
    return match.genes[0].direction * match.genes[1].direction


def isMonotonic(x):
    if len(x) < 3:
        return True

    dx = []
    for i in xrange(len(x)-1):
        dx.append(x[i+1] - x[i])
    return min(dx) >= 0 or max(dx) <= 0
    

def fillBlocks(conf, blocks, matches, quadtrees, lookup):
    # sort matches by score
    matches.sort(lambda a,b: cmp(b.score, a.score))
    
    # create block coords
    coords1 = {}
    coords2 = {}
    orders = {}
    inblock = {}
    for block in blocks:
        bdir = block.getDirection()
        block.matches.sort(lambda a,b: matchSort(a,b,bdir))
        
        c1 = []
        c2 = []
        for match in block.matches:
            c1.append(match.genes[0].start)
            c2.append(match.genes[1].start)
        coords1[block] = c1
        coords2[block] = c2
        orders[block] = bdir
    
    
    # find only matches near blocks
    for match in matches:
        chrom1 = match.genes[0].chrom
        chrom2 = match.genes[1].chrom
        
        # skip matches that have no blocks in same chrom pair
        # and matches that are already in blocks
        if not chrom1 in quadtrees or not chrom2 in quadtrees[chrom1] or \
           match.block != None:
            continue
        
        # get near blocks
        nearBlocks = quadtrees[chrom1][chrom2]. \
                        query(algorithms.Rect(
                            match.genes[0].start - conf['range1'],
                            match.genes[1].start - conf['range2'],
                            match.genes[0].end + conf['range1'],
                            match.genes[1].end + conf['range2']))
        
        # try to add match to near blocks
        for block in nearBlocks:
            # require that the match has the same direction as the block
            if matchDirection(match) != block.getDirection():
                continue
            
            # use binary search to find where match goes in block
            low1, top1 = algorithms.binsearch(coords1[block], 
                                              match.genes[0].start,
                                              order=1)
            low2, top2 = algorithms.binsearch(coords2[block], 
                                              match.genes[1].start,
                                              order=orders[block])
            
            # if match fits colinearly in block
            # match fits in block if low and tops are the same
            if low1 == low2 and top1 == top2:
                x1 = min(match.genes[0].start, block.end1)
                y1 = min(match.genes[1].start, block.end2)
                x2 = max(match.genes[0].end, block.start1)
                y2 = max(match.genes[1].end, block.start2)
                
                # ensure we are not passing through another block to make
                # this connection
                if len(quadtrees[chrom1][chrom2].query(
                          algorithms.Rect(x1, y1, x2, y2))) > 1:
                    continue
                
                # insert match
                match.filled = True
                block.add(match)
                lookup[match] = block
                if top1 != None:
                    coords1[block].insert(top1, match.genes[0].start)
                    coords2[block].insert(top1, match.genes[1].start)
                else:
                    coords1[block].append(match.genes[0].start)
                    coords2[block].append(match.genes[1].start)
                
                assert isMonotonic(coords1[block]), \
                    (low1, top1, coords1[block], orders[block])
                assert isMonotonic(coords2[block]), \
                    (low2, top2, coords2[block], orders[block])
                
                # stop searching blocks
                break



def filterRedundantBlocks(conf, blocks):
    blockList1 = util.groupby(lambda x: x.chrom1, blocks)
    blockList2 = util.groupby(lambda x: x.chrom2, blocks)
    blocks2 = []
    
    
    for block1 in blocks:
        within1 = False
        within2 = False

        for block2 in blockList1[block1.chrom1]:
            if block1.start1 >= block2.start1 and \
               block1.end1 <= block2.end1 and \
               block1.length1() < block2.length1() and \
               block1 != block2:
                within1 = True
                break

        for block2 in blockList2[block1.chrom2]:
            if block1.start2 >= block2.start2 and \
               block1.end2 <= block2.end2 and \
               block1.length2() < block2.length2() and \
               block1 != block2:
                within2 = True
                break

        if not within1 or not within2:
            blocks2.append(block1)
                
    
    return blocks2


                    
def makeMatchQuadtrees(matches):
    quadtrees = {}
    
    for match in matches:
        # create chrom-chrom key
        chrom1 = match.genes[0].chrom
        chrom2 = match.genes[1].chrom
        
        if not chrom1 in quadtrees:
            quadtrees[chrom1] = {}
        if not chrom2 in quadtrees[chrom1]:
            quadtrees[chrom1][chrom2] = \
                algorithms.QuadTree(
                    chrom1.size/2.0,
                    chrom2.size/2.0, 
                    max(chrom1.size/2.0,
                        chrom2.size/2.0))
        
        quadtrees[chrom1][chrom2].insert(match, 
            algorithms.Rect(match.genes[0].start, 
                            match.genes[1].start, 
                            match.genes[0].end, 
                            match.genes[1].end))
    return quadtrees


def makeBlockQuadtrees(blocks):
    quadtrees = {}
    
    for block in blocks:
        # create chrom-chrom key
        chrom1 = block.chrom1
        chrom2 = block.chrom2
        
        if not chrom1 in quadtrees:
            quadtrees[chrom1] = {}
        if not chrom2 in quadtrees[chrom1]:
            quadtrees[chrom1][chrom2] = \
                algorithms.QuadTree(
                    chrom1.size/2.0,
                    chrom2.size/2.0, 
                    max(chrom1.size/2.0,
                        chrom2.size/2.0))
        
        quadtrees[chrom1][chrom2].insert(block, 
            algorithms.Rect(block.start1, 
                            block.start2, 
                            block.end1, 
                            block.end2))
    return quadtrees

def installSynteny(matching, blocks):
    matching.blocks = blocks
    matching.matches = []
    
    i = 0
    for block in blocks:
        block.blockid = i
        i += 1
        for match in block.matches:
            match.block = block
            matching.matches.append(match)
    
    

def findSyntenyComponents(syngenes, useMatches=False):
    # read whole graph
    mat = util.Dict(2)
    
    if useMatches:
        for match in syngenes:
            mat[match.genes[0].name][match.genes[1].name] = 1
            mat[match.genes[1].name][match.genes[0].name] = 1
    else:
        for match in syngenes:
            mat[match[0]][match[1]] = 1
            mat[match[1]][match[0]] = 1
    
    # find connected components in graph
    def getNeighbor(vertex):
        return mat[vertex].keys()

    components = graph.connectedComponents(mat.keys(), getNeighbor)

    return components


def matchesFromComponents(comps):
    matches = []
    for comp in comps:
        for i in range(len(comp)):
            for j in range(i+1, len(comp)):
                matches.append([comp[i], comp[j]])
    return matches

    
def writeMatches(f, matches):
    for match in matches:
        print >>f, match.genes[0].name, match.genes[1].name



def avgGeneSpan(genome, ngenes = 5):
    dists = []
    for chrom in genome.chroms.itervalues():
        end = None
        
        i = 0; j = 0
        while True:
            # skip overlapping genes and little orfs
            while j < len(chrom.genes) and \
                chrom.genes[j].start - chrom.genes[i].end < 1000:
                    j += 1
            
            if j >= len(chrom.genes):
                break
            
            dists.append(chrom.genes[j].start - chrom.genes[i].start)
            i = j
    
    avgspan = stats.mean(dists) * ngenes
    util.log("%d bases per %d genes" % (avgspan, ngenes))
    
    return avgspan





def matchSort(a, b, bdir):
    val = cmp(a.genes[0].start, b.genes[0].start)
    
    if val == 0:
        if bdir == 1:
            return cmp(a.genes[1].start, b.genes[1].start)
        else:
            return cmp(b.genes[1].start, a.genes[1].start)
    else:
        return val





def verifySynteny(conf, blocks):
    for block in blocks:
        bdir = block.getDirection()
        block.matches.sort(lambda a,b: matchSort(a,b,bdir))
        
        
        coords1 = map(lambda x: x.genes[0].start, block.matches)
        coords2 = map(lambda x: x.genes[1].start, block.matches)
        
        assert isMonotonic(coords1)
        assert isMonotonic(coords2)
        
        for match in block.matches:           
            assert matchDirection(match) == bdir, "%s %s" % \
                    (matchDirection(match), bdir)

            
def displayQuickStats(conf, matching):
    keys = util.sort(matching.genomes.values(), lambda a,b: cmp(a.name, b.name))

    # count syntenic genes
    genecount = util.Dict(1, {})
    for match in matching.matches:
        for gene in match.genes:
            genecount[gene.chrom.genome.name][gene.name] = 1
    ngenes = util.mapdict(genecount, valfunc=len)
    
    
    # calc block coverage
    util.tic("calculating coverage")
    util.log("genome   #genes #matched  percent")
    for genome in keys:
        util.log("%8s %6s %8d %8.1f" % (
                 genome.name, 
                 len(genome.genes), 
                 ngenes[genome.name], 
                 100 * ngenes[genome.name] / float(len(genome.genes))))
    util.toc()


# TODO: not complete
def displayFullStats(conf, matching):
    keys = util.sort(matching.genomes.values(), lambda a,b: cmp(a.name, b.name))

    # count syntenic genes
    ngenes = {}
    cov = {}
    for genome in matching.genomes.values():
        ngenes[genome.name] = 0
        cov[genome.name] = 0
        for gene in genome.genes.values():
            if len(gene.matches) > 0:
                ngenes[genome.name] += 1
    
    comps = findSyntenyComponents(matching.matches, True)
    matching.setGeneComponents(comps)
    for genome in keys:        
        multiblocks = makeGenomeMultiBlocks(matching, genome)
        for multiblock in multiblocks:
            if len(multiblock.segments) > 1:
                cov[genome.name] += multiblock.segments[0].end - \
                                    multiblock.segments[0].start + 1 
    
    # calc block coverage
    util.log("genome   #genes #matched  percent coverage")
    util.log("%8s %6s %8d %8.1f %8.1f" % (
                 genome.name, 
                 len(genome.genes), 
                 ngenes[genome.name], 
                 100 * ngenes[genome.name] / float(len(genome.genes)),
                 100 * float(cov[genome.name]) / genome.size))
     

def findAnchors(conf, matching, outprefix):
    util.tic("find anchor matches (best bidiretional)")
    bbh = bestBidir(conf, matching.genomes.values(), matching.matches)
    bbh = filter(lambda match: match.score > conf["sigBbhScore"], bbh)
    genomeio.writeMatches(outprefix + ".bbh.match", bbh, True)
    matches = bbh
    util.log("anchor matches", len(bbh))
    util.toc()
    
    
    if "justBbh" in conf:
        util.toc()
        matching.matches = bbh
        return
    
    util.tic("find significant matches")
    best = bestUnidir(matching.genomes.values(), matching.matches)
    matches2 = filter(lambda match: match.score > conf["sigScore"], 
                      best.itervalues())
    
    
    # merge matches
    bestset = set(matches)
    for match in matches2:
        if match not in bestset:
            matches.append(match)
    
    util.log("anchor matches", len(matches))
    util.toc()
    
    return matches, best
    
    

def findSynteny(conf, matching, anchors, outprefix, dofilter=True):
    """
    Finds pairwise synteny
    """
    
    util.tic("find syteny")
    util.log("total input anchors", len(anchors))
    util.log("total input matches", len(matching.matches))
    
    # set configuration
    conf['range1'] = avgGeneSpan(anchors[0].genes[0].chrom.genome, 
                                 conf['rangeGenes'])
    conf['range2'] = avgGeneSpan(anchors[0].genes[1].chrom.genome, 
                                 conf['rangeGenes'])
    conf['range1'] =  max(conf['range1'], conf['range2'])
    conf['range2'] =  max(conf['range1'], conf['range2'])
    conf["blockCloseMerge1"] =  conf['range1']
    conf["blockCloseMerge2"] =  conf['range2']

        
    # make every match its own block    
    blocks = []
    for match in anchors:
        block = SyntenyBlock()
        block.init(match)
        blocks.append(block)
    
    
    util.tic("connect almost linear blocks")
    quadtrees = makeBlockQuadtrees(blocks)    
    mergeCloseBlocks(conf, blocks, quadtrees)
    util.log("number of blocks", len(blocks))
    util.toc()
    
    util.tic("connect blocks")
    quadtrees = makeBlockQuadtrees(blocks)
    connectBlocks(conf, blocks, quadtrees)
    util.log("number of blocks", len(blocks))
    util.toc()
    
    if dofilter:
        util.tic("fill blocks")
        lookup = makeBlockLookup(blocks)
        quadtrees = makeBlockQuadtrees(blocks)
        before = len(lookup)
        fillBlocks(conf, blocks, matching.matches, quadtrees, lookup)
        util.log("syntenic matches (after)", len(lookup))
        util.log("matches found", len(lookup) - before)
        util.toc()
        
        util.tic("filter blocks: minsize = %d" % conf["minBlockSize"])
        blocks = filter(lambda x: len(x.matches) >= conf["minBlockSize"], 
                        blocks)
        util.log("number of blocks", len(blocks))
        util.toc()
        
        util.tic("reconnect blocks")
        quadtrees = makeBlockQuadtrees(blocks)
        connectBlocks(conf, blocks, quadtrees)
        util.log("number of blocks", len(blocks))
        util.toc()
    
    # filter redundant blocks
    util.tic("filter redundant blocks")
    util.log("blocks before: ", len(blocks))
    util.log("matches before:", sum(map(lambda x: len(x.matches), blocks)))
    blocks = filterRedundantBlocks(conf, blocks)
    util.log("blocks after:", len(blocks))
    util.log("matches after:", sum(map(lambda x: len(x.matches), blocks)))
    util.toc()
    
    verifySynteny(conf, blocks)
    
    util.tic("install answer into data structure")
    installSynteny(matching, blocks)
    util.toc()
    
    util.toc()



def findAllPairsSynteny(hitdata, conf, outputDir, outprefix, outsuffix=""):
    util.tic("finding all pair-wise synteny")
    for genome1, genome2, genomefile1, genomefile2, \
        fasta1, fasta2, matchfile in hitdata:
        
        util.tic("processing %s, %s" % (genome1, genome2))

        prefix = outputDir + genome1 + "_" + genome2 + outsuffix

        # allow resuming a previously halted run
        if conf["synteny_resume"] and \
           os.path.exists(prefix + ".syncomps") and \
           os.path.exists(prefix + ".synteny") and \
           os.path.exists(prefix + ".blocks"):
            util.log("skipping... output exists")
            util.toc()
            continue

        matching = Matching()

        # read genomes
        util.tic("read genomes")
        matching.readGenomes(genomefile1, conf["gene2species"])
        matching.readGenomes(genomefile2, conf["gene2species"])
        matching.autoconf()
        util.toc()

        # read matches
        # filter potential matches by relative score
        util.tic("read matches")
        matching.readMatches(file(matchfile), conf["cols"], conf["minscore"])                     
        util.toc()
        
        anchors, best = findAnchors(conf, matching, prefix)
        matching.matches = filter(lambda match:
            match.score > best[match.genes[0]].score * conf["relScore"] or 
            match.score > best[match.genes[1]].score * conf["relScore"],
            matching.matches)
        
        # find synteny
        findSynteny(conf, matching, anchors, prefix)
        displayQuickStats(conf, matching)
        
        # save result to files
        util.tic("write output")
        genomeio.writeMatches(prefix + ".syncomps", matching.matches)
        genomeio.writeSynteny(prefix + ".synteny", matching.blocks)
        matching.writeSynteny(prefix + ".blocks")
        util.toc()

        util.toc()
    util.toc()



def refineSynteny(hitdata, conf, comps, outputDir, outprefix, outsuffix=""):
    util.tic("refining synteny")
    for genome1, genome2, genomefile1, genomefile2, \
        fasta1, fasta2, matchfile in hitdata:
        
        util.tic("refining %s, %s" % (genome1, genome2))

        prefix = outputDir + genome1 + "_" + genome2 + outsuffix

        # DEBUG: force no skipping with False
        # allow resuming a previously halted run
        if conf["synteny_resume"] and \
           os.path.exists(prefix + ".syncomps") and \
           os.path.exists(prefix + ".synteny") and \
           os.path.exists(prefix + ".blocks") and False:
            util.log("skipping... output exists")
            util.toc()
            continue

        matching = Matching()

        # read genomes
        util.tic("read genomes")
        matching.readGenomes(genomefile1, conf["gene2species"])
        matching.readGenomes(genomefile2, conf["gene2species"])
        matching.autoconf()
        util.toc()

        # determine anchors
        util.tic("determine anchors from components")
        gene2species = conf["gene2species"]
        for comp in comps:
            groups = util.groupby(gene2species, comp)
            if genome1 not in groups or \
               genome2 not in groups:
                continue
            
            for gene1 in groups[genome1]:
                for gene2 in groups[genome2]:
                    matching.addMatch(matching.genes[gene1], 
                                      matching.genes[gene2], 0)
        util.log("total anchors", len(matching.matches))
        util.toc()
        anchors = matching.matches
        
        # find synteny
        findSynteny(conf, matching, anchors, prefix, dofilter=False)
        
        displayQuickStats(conf, matching)
        
        # save result to files
        util.tic("write output")
        genomeio.writeMatches(prefix + ".syncomps", matching.matches)
        genomeio.writeSynteny(prefix + ".synteny", matching.blocks)
        matching.writeSynteny(prefix + ".blocks")
        util.toc()

        util.toc()
    util.toc()


def partOrths(conf, seqs, stree):
    util.tic("partition component of size %d" % len(seqs))
    
    #aln = muscle.muscleFast(seqs)
    #trees = phylip.bootNeighbor(seqs, conf["bootiters"])
    #
    #for tree in trees:
    #    parttrees = phylo.partitionTree(tree, stree, gene2species)
    
    #tree = phylip.proml(aln) #bionj.bionj(aln)
    
    
    tree = muscle.buildTree(seqs, verbose=False)
    tree = phylo.reconRoot(tree, stree, conf["gene2species"])
    tree.write(util.globalTimer())
    
    parttrees = phylo.partitionTree(tree, stree, conf["gene2species"])
    
    comps = []
    for tree2 in parttrees:
        comps.append(tree.leafNames())
    
    util.log("broken into %d parts" % len(comps))
    if len(comps) > 1:
        util.log("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        
    util.toc()
    
    return comps



def partOrthoComponents(hitdata, conf, comps, outputDir, outprefix):
    util.tic("partition components into sets of othologs")

    comps2 = []
    
    # get genomes
    genomes = util.unique(util.flatten(util.cget(hitdata, 0, 1)))
    
    # read all sequence data
    util.tic("read sequences")
    seqs = fasta.FastaDict()
    fastas = util.unique(util.flatten(util.cget(hitdata, 4, 5)))
    for f in fastas:
        seqs.read(f)
    util.toc()
    
    # load sprecies tree (use only the subtree needed for these genomes)
    stree = treelib.Tree()
    stree.readNewick(env.findFile(conf["species_tree"]))
    stree = treelib.subtree(stree, stree.lca(genomes))
    
    # filter components
    for i in xrange(len(comps)):
        comp = comps[i]
        
        genomes = map(conf["gene2species"], comp)
        counts = util.histDict(genomes).values()
        
        if len(filter(lambda x: x>1, counts)) >= 2:
            util.log("partitioning component %d of %d" % (i, len(comps)))
            
            # partition componet
            seqs2 = seqs.get(comp)
            
            if len(seqs2) > 2:
                comps2.extend(partOrths(conf, seqs2, stree))
            else:
                comps2.append(comp)
        else:
            # accept component without changes
            comps2.append(comp)
    
    util.toc()
    
    return comps2


def mergePairsSynteny(hitdata, conf, outputDir, outprefix, insuffix="",
                      outsuffix=""):
    util.tic("merging pair-wise synteny into multi-way")
    
    # remember output files
    syntenyfiles = []
    syngenefiles = []
    genomefiles = []
    for genome1, genome2, genomefile1, genomefile2, \
        fasta1, fasta2, matchfile in hitdata:
        
        genomefiles.append(genomefile1)
        genomefiles.append(genomefile2)
        prefix = outputDir + genome1 + "_" + genome2 + insuffix
        syntenyfiles.append(prefix + ".synteny")
        syngenefiles.append(prefix + ".syncomps")
    genomefiles = util.unique(genomefiles)
    
    # complete syntenic gene components
    util.tic("merge synteny components")
    syngenes = []
    for f in syngenefiles:
        syngenes.extend(util.readDelim(f))
    comps = findSyntenyComponents(syngenes)
    del syngenes # free-up unneeded memory 
    util.toc()

    # output stats
    util.log("total # of genes: %d " % sum(map(len, comps)))
    util.log("largest component: %d" % max(map(len, comps)))
    
    
    # write mutli-way output
    util.writeDelim(file(outputDir + outprefix + outsuffix + ".syncomps", "w"), 
                    comps)
    # cat all synteny into one file
    out = file(outputDir + outprefix + outsuffix + ".synteny", "w")
    for f in syntenyfiles:
        infile = file(f)
        for line in infile:
            out.write(line)
    out.close()
    
    util.toc()


def generateMultiBlocks(hitdata, conf, comps, outputDir, outprefix):
    # get genome files
    genomefiles = []
    for genome1, genome2, genomefile1, genomefile2, \
        fasta1, fasta2, matchfile in hitdata:
        
        genomefiles.append(genomefile1)
        genomefiles.append(genomefile2)
    genomefiles = util.unique(genomefiles)
    
    
    # multiblocks
    util.tic("create multiblocks")
    util.tic("reading genomes")
    matching = Matching()
    for genomefile in genomefiles:
        util.tic("reading "+ os.path.basename(genomefile))
        matching.readGenomes(genomefile, conf["gene2species"])
        util.toc()
    matching.autoconf()
    util.toc()
    
    matching.setGeneComponents(comps)
    genomeio.readBlockDimensions(matching, outputDir + outprefix + ".synteny")
    multiblocks = makeGenomeMultiBlocks({}, matching, 
                                        matching.genomes[conf["refgenome"]])
    writeMultiBlocks(outputDir + outprefix + ".multiblocks", multiblocks)
    util.toc()

    

def findMultiSynteny(hitdata, conf = initConf(), outputDir = "./", 
                     outprefix = "all"):
    """Find multigenome synteny
    
     example hit data:
     hitdata = [("dog", "human", "dog.coord", "human.coord", "dog_human.blastp"),
                ...]
    """
    
    if False:
        # find synteny for each pair-wise speciation comparison    
        findAllPairsSynteny(hitdata, conf, outputDir, outprefix, "")

        # merge pair-wise synteny into draft multi-way
        mergePairsSynteny(hitdata, conf, outputDir, outprefix, "",
                          outsuffix="")
        comps = util.readDelim(outputDir + outprefix + ".syncomps")

        # generate multiblocks
        generateMultiBlocks(hitdata, conf, comps, outputDir, outprefix)
    
    
    
    if True:
        # find synteny for each pair-wise speciation comparison    
        findAllPairsSynteny(hitdata, conf, outputDir, outprefix, ".draft")

        # merge pair-wise synteny into draft multi-way
        mergePairsSynteny(hitdata, conf, outputDir, outprefix, 
                          insuffix=".draft", outsuffix=".draft")
        comps = util.readDelim(outputDir + outprefix + ".draft.syncomps")
        util.log("%d components" % len(comps))
        
        #comps = partOrthoComponents(hitdata, conf, comps, outputDir, outprefix)
        #util.log("after partitioning %d components" % len(comps))
        #util.writeDelim(outputDir + outprefix + ".parted.syncomps", comps)
        
        # refine synteny
        refineSynteny(hitdata, conf, comps, outputDir, outprefix, "")

        # merge pair-wise synteny into refined multi-way
        mergePairsSynteny(hitdata, conf, outputDir, outprefix, 
                          insuffix="", outsuffix="")

        # generate multiblocks
        generateMultiBlocks(hitdata, conf, comps, outputDir, outprefix)
        

    return comps


"""

def splitBlock(conf, block, lookup):
    direction = block.matches[0].direction()
    
    def helper(side):
        breaks = []
        
        # walk down the genes in the range of a block
        if side == 1:
            geneIter = GeneIter(block.chrom1.genes, block.start1, block.end1)
        else:
            geneIter = GeneIter(block.chrom2.genes, block.start2, block.end2)

        interupt = {}
        broken = False
        last = 0
    
        # find breaks
        while geneIter.more():
            gene = geneIter.get()
        
            for match in gene.matches:
                if match in lookup:
                    otherBlock = lookup[match]

                    if otherBlock != block:
                        if not otherBlock in interupt:
                            interupt[otherBlock] = 1
                        else:
                            interupt[otherBlock] += 1

                        if interupt[otherBlock] > conf["blockInterupt"] and \
                           not broken:
                            if not direction and side == 2:
                                breaks.append(last-1)
                            else:
                                breaks.append(last)
                            broken = True
                    else:
                        # reset interupts
                        interupt = {}
                        broken = False
                        last = match.genes[0].start
        return breaks
    
    # break in both directions
    breaks = helper(1)
    breaks += helper(2)
    breaks.sort()
    
    # create new blocks
    block.matches.sort(lambda a,b: cmp(a.genes[0].start, b.genes[0].start))
    
    blocks2 = [SyntenyBlock()]
    blocks2[0].add(block.matches[0])
    i = 0
    for match in block.matches[1:]:
        if i >= len(breaks):
            # no more breaks
            blocks2[-1].add(match)
        else:
            if breaks[i] < match.genes[0].start:
                # break block
                blocks2.append(SyntenyBlock())
                blocks2[-1].init(match)
                
                while i < len(breaks) and breaks[i] < match.genes[0].start:
                    i += 1
            else:
                # grow block
                blocks2[-1].add(match)
    
    return blocks2
                

def splitBlocks(conf, blocks, matches, lookup):    
    blocks2 = []
    
    for block in blocks:
        blocks2 += splitBlock(conf, block, lookup)
    
    return blocks2

def findHoles(matching, genes, comps, radius, logfile):
    # make gene to components lookup
    lookup = {}
    for i in xrange(len(comps)):
        for gene in comps[i]:
            lookup[genes[gene]] = i
        

    # scan genes for holes
    util.tic("scan for holes")
    holes = {}
    for genome in matching.genomes.values():
        for chrom in genome.chroms.values():
            for i in xrange(len(chrom.genes)):
                gene = chrom.genes[i]
            
                # only process genes without synteny
                if gene in lookup:
                    continue
                
                # get genes nearby in chrom
                left = chrom.genes[i-radius:i]
                right = chrom.genes[i+1:i+radius+1]
                
                # find all genes nearby in the alignment
                nearGenes = []
                for j in util.subdict(lookup, left + right).values():
                    nearGenes.extend([genes[x] for x in comps[j]])
                
                # get a histogram of the chromosomes nearby
                chromHist = util.histDict([x.chrom for x in nearGenes])
                if gene.chrom in chromHist:
                    del chromHist[gene.chrom]                
                
                # record locations of neighbors, keyed by chrom
                locs = {}
                for gene2 in nearGenes:
                    locs.setdefault(gene2.chrom, []).append(gene2.start)
                
                # gene is hole if another chrom occurs atleast half the time
                #if len(filter(lambda x: x >= radius, chromHist.values())) > 1:
                if True:
                    holes[gene.name] = locs
                    print >>logfile, "possible hole %s" % gene.name
    util.log("%d holes" % len(holes))
    util.toc()

    return holes


def findNearSynteny(matching, comps, genes, data, conf, logfile):
    # find holes (unmatched genes)
    holes = findHoles(matching, genes, comps, conf["nearRange"], logfile)
        
    util.tic("filter matches")
    
    keep = {}
    best = util.Dict(1, [None, 0])    
    spans = [avgGeneSpan(matching.genomes[x[0]], conf["nearRange"]) for x in data]
    fillRadius = stats.mean(spans)
    
    def testGene(gene1, gene2, score):
        # only process genes that we know about and that are holes with nearby
        # synteny
        if gene1 in genes and gene2 in genes and \
           gene1 in holes and \
           genes[gene2].chrom in holes[gene1]:
           
            # see if other gene is near one of the matched near genes
            for chrom, locs in holes[gene1].iteritems():
                if genes[gene2].chrom == chrom:
                    for loc in locs:
                        if abs(genes[gene2].start - loc) < fillRadius:
                            print >>logfile, "near match %s %s %f" % \
                                (gene1, gene2, score)
                            keep.setdefault(gene1, []).append((gene2, score))
                            return
    
    cols = conf["cols"]
    # loop through match files to find near-syntenic matches
    for genome1, genome2, genomefile1, genomefile2, matchfile in data:
        for line in file(matchfile):
            tokens = line.split()
            gene1 = tokens[cols[0]]
            gene2 = tokens[cols[1]]
            score = float(tokens[cols[2]])
            
            # test whether this match is near synteny
            testGene(gene1, gene2, score)
            testGene(gene2, gene1, score)
            
            # keep track of whether this match is a BUD
            if score > best[gene1][1]:
                best[gene1] = [gene2, score]
            if score > best[gene2][1]:
                best[gene2] = [gene1, score]                
    util.toc()
    
    
    # keep gene's whose best match is near syntentic    
    near = []
    for gene in keep:
        if best[gene][0] in [x[0] for x in keep[gene]]:
            near.append([gene, best[gene][0], best[gene][1]])

    return near


def filterComponents(conf, genes, comps):
    for i in xrange(len(comps)):
        maxsizes = util.Dict(1, 0)
        for gene in comps[i]:
            if genes[gene].length() > maxsizes[genes[gene].chrom.genome.name]:
                maxsizes[genes[gene].chrom.genome.name] = genes[gene].length()
    
        comps[i] = filter(lambda gene:
            genes[gene].length() >= conf["minGeneSize"] and 
            genes[gene].length() >= maxsizes[genes[gene].chrom.genome.name] * 
                                    conf["minPullLen"],
            comps[i])
            
    return filter(lambda x: len(x) > 1, comps)
            
    
def refineSynteny(matching, comps, data, conf, outputDir = "./", 
    outprefix="all", logfile = sys.stderr):
    util.tic("refine synteny")
    
    util.log("%d matches" % len(matching.matches))

    genes = matching.getGenes()
    
    # find nearly syntenic genes
    near = findNearSynteny(matching, comps, genes, data, conf, logfile)
    
    # add near matches to matching
    for gene1, gene2, score in near:
        matching.addMatch(genes[gene1], genes[gene2], score)
    
    
    # print out genes whose best match is near syntentic
    out = file(outputDir + outprefix + ".synnear", "w")
    for gene1, gene2, score in near:
        print >>out, "\t".join([gene1, gene2])
    out.close()
    util.log("wrote %s.synnear" % outprefix)
    near = [] # free memory
    
    # print out matches
    genomeio.writeMatches(outputDir + outprefix + ".refine.syngenes", 
                          matching.matches)
    util.log("wrote %s.refine.syngenes" % outprefix)  
    matching.matches = [] # free memory
    
    
    # find synteny components and add implied matches
    util.tic("find new components")
    os.system("unionpart.py -p %s.refine.syngenes > %s.refine.syncomps" % \
        (outputDir + outprefix, outputDir + outprefix))
    comps = util.readDelim(outputDir + outprefix + ".refine.syncomps")
    comps = filterComponents(conf, genes, comps)
    util.toc()
    
    util.writeDelim("%s.refine.syncomps" % (outputDir + outprefix), comps)
    

    # rebuild synteny blocks using additional genes
    util.tic("rebuild synteny blocks")    
    
    # setup some configuration
    # TODO: need to do something smarter.  Should find ranges for each pair of 
    # genomes
    ranges = map(lambda x: avgGeneSpan(x, conf['rangeGenes']),
                 matching.genomes.values())
    conf['range1'] = max(ranges)
    conf['range2'] = max(ranges)
    conf["blockCloseMerge1"] =  conf['range1']
    conf["blockCloseMerge2"] =  conf['range2']


    # rebuild synteny for each pair of genomes separately
    # because memory is too limited to do them all at once    
    genomes = matching.genomes.values()
    syntenyfiles = []
    for g1 in xrange(len(genomes)):
        for g2 in xrange(g1+1, len(genomes)):
            genome1 = genomes[g1]
            genome2 = genomes[g2]
            genomePair = {genome1:1, genome2:1}
            
            
            # set configuration
            conf['range1'] = avgGeneSpan(genome1, conf['rangeGenes'])
            conf['range2'] = avgGeneSpan(genome2, conf['rangeGenes'])
            conf["blockCloseMerge1"] =  conf['range1']
            conf["blockCloseMerge2"] =  conf['range2']
            
            
            util.tic("rebuild %s, %s" % (genome1.name, genome2.name))
            util.tic("convert components into matches")
            matches = []
            for comp in comps:
                for i in xrange(len(comp)):
                    for j in xrange(i+1, len(comp)):
                        if genes[comp[i]].chrom.genome in genomePair and \
                           genes[comp[j]].chrom.genome in genomePair and \
                           genes[comp[i]].chrom.genome != \
                           genes[comp[j]].chrom.genome:
                            matches.append(Match(genes[comp[i]], genes[comp[j]]))
            util.toc()

            # init every match into a block
            blocks = []
            for match in matches:
                match.genes.sort(lambda a,b: cmp(a.chrom.name, b.chrom.name))
                block = SyntenyBlock()
                block.init(match)
                blocks.append(block)

            # don't keep self-blocks (meaningless)
            blocks = filter(lambda x: x.genome1 != x.genome2, blocks)

            util.log("%d matches" % len(matches))

            util.tic("make quadtrees of blocks")
            quadtrees = makeBlockQuadtrees(blocks)
            util.toc()    

            util.tic("connect blocks")
            connectBlocks(conf, blocks, quadtrees)
            util.log("number of blocks", len(blocks))
            util.toc()


            util.tic("make quadtrees of new blocks")
            quadtrees = makeBlockQuadtrees(blocks)
            util.toc()
            util.tic("connect almost linear blocks")
            mergeCloseBlocks(conf, blocks, quadtrees)
            util.log("number of blocks", len(blocks))
            util.toc()    

            util.log("%d blocks" % len(blocks))

            util.tic("install answer into data structure")
            installSynteny(matching, blocks)
            util.toc()
            
            filename = outputDir + genome1.name+"_"+genome2.name + ".refine.synteny"
            syntenyfiles.append(filename)
            genomeio.writeSynteny(filename, matching.blocks)
            util.log("wrote %s" % filename)
            
            util.toc()
    util.toc()
    
    
    # cat all synteny into one file
    out = file(outputDir + outprefix + ".refine.synteny", "w")
    for f in syntenyfiles:
        infile = file(f)
        for line in infile:
            out.write(line)
    out.close()
    
    util.toc()



def pullGenes(genes, comps, data, conf, logfile = sys.stderr):
    util.tic("pull genes")

    # make component lookup
    lookup = {}
    for i in xrange(len(comps)):
        for gene in comps[i]:
            lookup[gene] = i

    
    # find best matches and inner-component matches
    best = {}
    compscores = util.Dict(1, [])
    cols = conf["cols"]
    for genome1, genome2, genomefile1, genomefile2, matchfile in data:
        for line in file(matchfile):
            tokens = line.split()
            gene1 = tokens[cols[0]]
            gene2 = tokens[cols[1]]
            score = float(tokens[cols[2]])
            
            if gene1 not in genes or gene2 not in genes:
                continue
            
            if gene1 in lookup and gene2 in lookup and \
               lookup[gene1] == lookup[gene2]:
                compscores[lookup[gene1]].append(score)
                continue
            
            if gene1 in lookup:
                gene1, gene2 = gene2, gene1
            
            if gene1 in lookup or gene2 not in lookup:
                continue
            
            if gene1 not in best or score > best[gene1][0]:
                best[gene1] = [score, gene2]
    
    
    # only pull best match if it is within the distribution
    # find median pair-wise score for each partition
    medians = []
    for i in xrange(len(comps)):
        if len(compscores[i]) == 0:
            compscores[i] = [conf["sigScore"]]
        medians.append(compscores[i][len(compscores[i])/2])
    
    def genelen(gene):
        return genes[gene].end - genes[gene].start + 1
    
    #medianlens = []
    #for comp in comps:
    #    lens = [genelen(gene) for gene in comp]
    #    medianlens.append(lens[len(lens) / 2])
    
    
    for gene1, match in best.iteritems():
        score, gene2 = match
        if score > medians[lookup[gene2]] * conf["pullFrac"]: # and \
           #genelen(gene1) > medianlens[lookup[gene2]] * conf["minPullLen"]:
            comps[lookup[gene2]].append(gene1)
    
    comps = filterComponents(conf, genes, comps)
    
    util.toc()
"""


