from genomeutil import *
import copy
import util
import sys
import graph

def bus(conf, matching):
    alg = Bus()
    alg.run(conf, matching)

def initConf(conf):
    conf['prune'] = .4
    conf['blockExtend'] = 3000
    conf['minBlockSize'] = 6
    conf['noprune'] = False

class Bus:
    timer = None

    def __init__(self):
        self.timer = util.Timer()

    def run(self, conf, matching):
        self.timer.start("BUS")
        if not conf['noprune']:
            self.pruneRelative(conf, matching)
        self.buildSynteny2(conf, matching)
        self.pruneSynteny(conf, matching)
        self.findGeneSets(conf, matching)
        self.pruneHomology(conf, matching)
        self.completeHomology(conf, matching)
        self.timer.stop()

    def pruneRelative(self, conf, matching):
        cutoff = conf['prune']
        
        self.timer.start("pruning")
        
        pruned = {}
        def setPrune(match):
            if not pruned.has_key(match):
                pruned[match] = True
            else:
                match.prune()

        # loop through all genes
        for genome in matching.genomes.values():
            for gene in genome.genes.values():
                # find best match for each otherGenome
                matches = {}
                for match in gene.matches:
                    otherGenome = match.otherGene(gene).chrom.genome
                    if (not matches.has_key(otherGenome)) or \
                       matches[otherGenome] < match.score:
                        matches[otherGenome] = match.score

                # prune matches if they are not above cutoff * best
                for match in gene.matches:
                    otherGenome = match.otherGene(gene).chrom.genome
                    if match.score < matches[otherGenome] * cutoff:
                        setPrune(match)

                if len(matches.values()) > 0:
                    # prune self matches with a stricter rule
                    selfCutoff = max(matches.values())
                    for match in gene.selfMatches:
                        otherGenome = match.otherGene(gene).chrom.genome
                        if match.score < selfCutoff:
                            setPrune(match)
        self.timer.stop()
    
    def rangeCollide(self, start1, end1, start2, end2):
        return end1 > start2 and start1 < end2
    
    def mergableBlocks(self, conf, block1, block2):
        size = conf['blockExtend']
        return self.rangeCollide(block1.start1 - size, block1.end1 + size, 
                                 block2.start1, block2.end1) and \
               self.rangeCollide(block1.start2 - size, block1.end2 + size, 
                                 block2.start2, block2.end2)
    
    def buildSynteny(self, conf, matching):
        self.timer.start("building synteny")
        
        blocks = []
        blockLists = {}
        
        # create one block for each match
        self.timer.log("create block lists")
        for match in matching.matches:
            if match.pruned or match.isSelfMatch():
                continue
        
            block = SyntenyBlock()
            block.init(match)
            blist = (block.chrom1, block.chrom2)
            if not blockLists.has_key(blist):
                blockLists[blist] = [block]
            else:
                blockLists[blist].append(block)
        
        # perform block merging
        for blockList in blockLists.values():
            self.timer.log("merge block list %d" % (len(blockList),))
        
            for i in [0, 1]:
                # sort blocks along either chrom1 or chrom2
                if i == 0:
                    blockList.sort(lambda a,b: b.start1 - a.start1)
                else:
                    blockList.sort(lambda a,b: b.start2 - a.start2)
                
                # walk down block list and compare to past blocks
                for i in range(len(blockList)):
                    for j in range(i-1, -1, -1):
                        if blockList[j] == None:
                            continue
                        
                        if self.mergableBlocks(conf, blockList[i], blockList[j]):
                            blockList[i].merge(blockList[j])
                            blockList[j] = None
                
                # remove merged blocks from blockList
                blockList = filter(lambda x: x != None, blockList)
        
            # filter blocks and add to master list
            for block in blockList:
                if len(block.matches) > conf['minBlockSize']:
                    blocks.append(block)
        
        # finalize blocks
        for i in range(len(blocks)):
            blocks[i].blockid = i
            for match in blocks[i].matches:
                match.block = blocks[i]
        matching.blocks = blocks
        
        # loop through syntenic matches and mark syntenic genes
        for block in matching.blocks:
            for match in block.matches:
                for gene in match.genes:
                    gene.isSyntenic = True
            
        self.timer.stop()
    
    

    # synteny block is connected components of match-box touch graph
    def buildSynteny2(self, conf, matching):
        self.timer.start("building synteny")                
        
        def getNeighbors(tree, match):
            n = tree.query(
                Rect(match.genes[0].start - conf['blockExtend']/2.0,
                     match.genes[1].start - conf['blockExtend']/2.0,
                     match.genes[0].end   + conf['blockExtend']/2.0,
                     match.genes[1].end   + conf['blockExtend']/2.0))
            return n

        def newComponent(component):
            blocks.append(SyntenyBlock())

        def visitMatch(match):
            blocks[-1].add(match)
        
        # mark best matches
        self.markBest(conf, matching)
        
        matchLists = {}
        
        # separate genes by chromosome pairs
        self.timer.log("separate matches")
        for match in matching.matches:
            if match.pruned or match.isSelfMatch() or not match.best:
                continue
            
            chrom1 = match.genes[0].chrom
            chrom2 = match.genes[1].chrom
            match.block = None
            
            if chrom1.order > chrom2.order:
                tmp = chrom1; chrom1 = chrom2; chrom2 = tmp 
            
            key = (chrom1, chrom2)
            if not matchLists.has_key(key):
                matchLists[key] = [match]
            else:
                matchLists[key].append(match)
        
        
        # perform block merging
        blocks = []
        for lstKey in matchLists:
            lst = matchLists[lstKey]
            self.timer.log("merge block list %d" % len(lst))
            
            # create quad tree
            tree = QuadTree(lstKey[0].size/2, lstKey[1].size/2,
                            max(lstKey[0].size, lstKey[1].size))
            
            for match in lst:
                tree.insert(match, Rect(match.genes[0].start,
                                        match.genes[1].start,
                                        match.genes[0].end,
                                        match.genes[1].end))
            graph.dfs(lst, visitMatch, 
                      lambda x: getNeighbors(tree, x), newComponent)
        
        # keep only blocks that are big enough
        blocks = filter(lambda x: len(x.matches) >= conf['minBlockSize'], blocks)
        
        # finalize blocks
        for i in range(len(blocks)):
            blocks[i].blockid = i
            for match in blocks[i].matches:
                match.block = blocks[i]
        matching.blocks = blocks
        
        # loop through syntenic matches and mark syntenic genes
        for block in matching.blocks:
            for match in block.matches:
                for gene in match.genes:
                    gene.isSyntenic = True
         
        self.timer.stop()
    
    
    def pruneSynteny(self, conf, matching):
        self.timer.start("prune based on synteny")
        
        # loop through genes, and if syntenic, remove non-syntenic matches
        for genome in matching.genomes.values():
            for gene in genome.genes.values():
                if gene.isSyntenic:
                    for match in gene.matches:
                        if match.block == None:
                            match.prune()
                
        self.timer.stop()
    
    
    def markBest(self, conf, matching):
        self.timer.start("find best subgraph")
        
        for match in matching.matches:
            match.best = False
        
        # mark matches as best or not
        for genome in matching.genomes.values():
            for gene in genome.genes.values():
                bestScore = {}
                best = {}
            
                # find match best for this gene with respect to each otherGenome
                for match in gene.matches:
                    otherGenome = match.otherGene(gene).chrom.genome
                    if not best.has_key(otherGenome) or \
                       bestScore[otherGenome] < match.score:
                        bestScore[otherGenome] = match.score
                        best[otherGenome] = match
                
                # mark matches
                for match in best.values():
                    match.best = True        
        
        self.timer.stop()
    
    
    def findGeneSets(self, conf, matching):
        self.timer.start("find gene sets")
    
        self.markBest(conf, matching)
        
        # find connected components for best matches (by DFS)
        self.timer.start("find connected components")
        
        # store all genes in one list
        matching.homology = []
        allGenes = []
        for genome in matching.genomes.values():
            allGenes += genome.genes.values()

        def getNeighbors(gene):
            lst = []
            for match in gene.allMatches():
                if match.best or match.block != None:
                    lst.append(match.otherGene(gene))
            return lst
        
        def newComponent(component):
            matching.homology.append(HomologyGroup(component))
        
        def visitGene(gene):
            gene.homology = matching.homology[-1]
            matching.homology[-1].genes.append(gene)
        
        graph.dfs(allGenes, visitGene, getNeighbors, newComponent)
                
        self.timer.stop()
        self.timer.stop()

        
    def pruneHomology(self, conf, matching):        
        # prune matches that span different components
        self.timer.start("prune based on homology")
        for match in matching.matches:
            if match.pruned:
                continue
            
            if match.genes[0].homology.groupid != \
               match.genes[1].homology.groupid:
                match.prune()
        self.timer.stop()
        
        # reattach orphans
        tmp = 0
        for genome in matching.genomes.values():
            for gene in genome.genes.values():
                if len(gene.allMatches()) == 0:
                    bestMatch = None
                    bestScore = 0
                    for match in gene.pruned:
                        if match.score > bestScore:
                            bestMatch = match
                            bestScore = match.score
                            
                    if bestMatch != None:
                        bestMatch.unprune()
                        gene.homology.genes.remove(gene)
                        gene.homology = bestMatch.otherGene(gene).homology
                        gene.homology.genes.append(gene)
                    else:
                        tmp += 1
                        print tmp

    
    def completeHomology(self, conf, matching):    
        # resurrect pruned matches that obey homology group
        def list2dict(lst):
            d = {}
            for item in lst:
                d[item] = True
            return d
        
        self.timer.start("complete homology groups")
        for group in matching.homology:
            hgenes = list2dict(group.genes)
            
            self.timer.log("homology group %d size %d" % \
                (group.groupid, len(group.genes)))
            
            if len(group.genes) > 50:
                continue
            
            for gene in group.genes:
                hgenes2 = copy.copy(hgenes)
                
                # remove all existing matches
                for gene2 in gene.neighbors():
                    try:
                        del(hgenes2[gene2])
                    except KeyError:
                        print gene.name, gene2.name, \
                            (gene2 in hgenes)
                
                # now insert any matches that don't exist yet
                for other in hgenes2:
                    # skip loops
                    if gene == other:
                        continue
                    
                    # see if pruned match can be resurrected
                    match = None
                    for match2 in gene.pruned:
                        if other == match2.otherGene(gene):
                            match = match2
                            match.unprune()
                            break
                    
                    # if no existing match then create new one
                    #if match == None:
                    #    match = Match(gene, other)
                    #    match.attach()
                    #    matching.matches.append(match)
        
        # remove empty homology groups, and renumber groups
        matching.homology = filter(lambda x: len(x.genes) > 0, matching.homology)
        i = 0
        for group in matching.homology:
            group.blockid = i
            i += 1
        
        self.timer.stop()
        
class Rect:
    x1 = 0
    y1 = 0
    x2 = 0
    y2 = 0        
    
    def __init__(self, x1, y1, x2, y2):
        if x1 < x2:
            self.x1 = x1
            self.x2 = x2
        else:
            self.x1 = x2
            self.x2 = x1
        if y1 < y2:
            self.y1 = y1
            self.y2 = y2
        else:
            self.y1 = y2
            self.y2 = y1

class QuadNode:
    item = None
    rect = None
    
    def __init__(self, item, rect):
        self.item = item
        self.rect = rect
        
        
class QuadTree:
    nodes = []
    children = []
    center = [0,0]
    size = 100
    depth = 0
    MAX = 10
    MAX_DEPTH = 10
    
    def __init__(self, x, y, size, depth = 0):
        self.nodes = []
        self.children = []
        self.center = [x, y]
        self.size = size
        self.depth = depth
    
    def insert(self, item, rect):
        if len(self.children) == 0:
            self.nodes.append(QuadNode(item, rect))
            
            if len(self.nodes) > self.MAX and self.depth < self.MAX_DEPTH:
                self.split()
        else:
            self.insertIntoChildren(item, rect)
    
    def insertIntoChildren(self, item, rect):
        if rect.x1 < self.center[0]:
            if rect.y1 < self.center[1]:
                self.children[0].insert(item, rect)
            if rect.y2 > self.center[1]:
                self.children[1].insert(item, rect)
        if rect.x2 > self.center[0]:
            if rect.y1 < self.center[1]:
                self.children[2].insert(item, rect)
            if rect.y2 > self.center[1]:
                self.children[3].insert(item, rect)
                   
    def split(self):
        self.children = [QuadTree(self.center[0] - self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] - self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1)]

        for node in self.nodes:
            self.insertIntoChildren(node.item, node.rect)
        self.nodes = []

    def query(self, rect, results = None):
        if results == None:
            results = []
    
        if len(self.children) > 0:
            if rect.x1 < self.center[0]:
                if rect.y1 < self.center[1]:
                    self.children[0].query(rect, results)
                if rect.y2 > self.center[1]:
                    self.children[1].query(rect, results)
            if rect.x2 > self.center[0]:
                if rect.y1 < self.center[1]:
                    self.children[2].query(rect, results)
                if rect.y2 > self.center[1]:
                    self.children[3].query(rect, results)
        else:
            for node in self.nodes:
                if node.rect.x2 > rect.x1 and node.rect.x1 < rect.x2 and \
                   node.rect.y2 > rect.y1 and node.rect.y1 < rect.y2:
                    results.append(node.item)
        return results
                    
    def getSize(self):
        size = 0
        for child in self.children:
            size += child.getSize()
        size += len(self.nodes)
        return size

def test():
    tree = QuadTree(0,0, 100)

    import random

    def rand():
        return random.random() * 200 - 100

    pts2 = []
    for i in range(200):
        x = rand()
        y = rand()
        pts2.append((i, x, y))
        tree.insert(i, Rect(x, y, x+1, y+1))

    pts = tree.query(Rect(100, 10, -100, -100))
    for i in pts:
        print pts2[i]



def combineHomologyGroups(files):
    def readHomology(vertices, labelFilename, groupsFilename):
        print "read", groupsFilename
        inlabels = file(labelFilename)
        ingroups = file(groupsFilename)    

        for label in inlabels:
            label = label.rstrip()
            group = "GROUP:" + groupsFilename + ingroups.next().rstrip()
            util.case(vertices, label, []).append(group)
            util.case(vertices, group, []).append(label)
        
    vertices = {}
    
    # read in all data
    for partFile, labelFile in files:
        readHomology(vertices, labelFile, partFile)
    
    # find connected components
    import graph
    def neighborFunc(vertex):
        return vertices[vertex]
    components = graph.connectedComponents(vertices, neighborFunc)
    
    # create one master label part file
    outLabels = []
    outPart = []
    for i in range(len(components)):
        for vertex in components[i]:
            if not vertex.startswith("GROUP"):
                outLabels.append(vertex)
                outPart.append(i)
    return (outPart, outLabels)
