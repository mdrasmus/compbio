
# python libs
import copy
import math
import sys

# rasmus libs
from rasmus import clustalw
from rasmus import matrix
from rasmus import stats
from rasmus import util
from rasmus.genomeutil import *



if __name__ != "__main__":
    from summon.core import *
    from summon import shapes


def setupTree(tree, genes):
    def walk(node):
        if node.isLeaf():
            if node.name in genes:
                node.gene = genes[node.name]
                return True
            else:
                node.gene = None
                return False
        else:
            node.gene = None
            removes = []
            for child in node.children:
                if not walk(child):
                    removes.append(child)
            
            # clean up tree
            for remove in removes:
                tree.remove(remove)
            return len(node.children) > 0
            
    walk(tree.root)
    
    tree.setSizes()


def setMix(node, genome1):
    if node.gene != None:
        if node.gene.chrom.genome == genome1:
            node.mix = 1
        else:
            node.mix = 0
    else:
        node.mix = 0
        for child in node.children:
            setMix(child, genome1)
            node.mix += child.mix * child.size / float(node.size)


class TreeVis:
    def __init__(self):
        self.mode = "desc"
        self.selgenes = []
        self.selalign = {}
        self.selnode = None
    

    def installBindings(self):
        def press(mode):
            def func():
                print "mode is '"+ mode + "'"
                self.mode = mode
            return func
        
        set_binding(input_key("a"), press("align"))
        set_binding(input_key("d"), press("desc"))

    def findGene(self, name):
        if name in self.tree.nodes:
            x = self.tree.nodes[name].x
            y = self.tree.nodes[name].y
            w = 5
            
            set_visible(x-w, y-w, x+w, y+w)
        else:
            print "cannot find gene '%s'" % name
            

    def nodeClick(self, node):
        self.selnode = node
    
        if self.mode == "desc":
            def walk(node):
                if node.gene != None:
                    return self.selgenes.append(node.gene)
                else:
                    for child in node.children:
                        walk(child)
            self.selgenes = []
            walk(node)

            # print descriptions
            for gene in self.selgenes:
                print gene.chrom.genome.name, gene.chrom.name, \
                      gene.name, gene.description
            

        elif self.mode == "align":
            # print alignment

            # get all protein sequences of sub tree
            def walk(node):
                if node.isLeaf():
                    return [(node.gene.name, node.gene.protein())]
                else:
                    seqs = []
                    for c in node.children:
                        seqs.extend(walk(c))
                    return seqs
            seqs = walk(node)
            seqs2 = {}

            for seq in seqs:
                seqs2[seq[0]] = seq[1]

            self.selalign = clustalw.clustalw(seqs2, verbose = True)
            clustalw.printAlign(self.selalign)
        

    def drawNode(self, conf, node, sx, sy):
        vis = []
        if "cut" in dir(node) and node.cut:
            vis.append(color(0,1,0))
        else:
            vis.append(color(node.mix, 0, 1 - node.mix))
        vis.append(lines(vertices(sx, sy, sx, sy-conf['height'])))
        
        #if node.cut:
        #    vis.append(color(0,1,0))
        #    vis.append(shapes.box(sx-1, sy-1, sx+1, sy+1))
        
        # record node position in tree
        node.x = sx
        node.y = sy - conf['height']

        def func():
            print "----------------"
            self.nodeClick(node)

        vis.append(hotspot("click", 
                           sx - node.size/2.0, sy,
                           sx + node.size/2.0, sy-conf['height'],
                           func))

        # draw horizontal line
        if len(node.children) > 0:
            left = sx - node.size/2.0 + node.children[0].size/2.0
            right = sx + node.size/2.0 - node.children[-1].size/2.0

            vis.append(lines(#color(node.mix, 0, 1 - node.mix), 
                             vertices(left, sy - conf['height'], 
                                      right, sy - conf['height'])))
        return list2group(vis)


    def drawTree(self, conf, node, drawNode = None, sx=0, sy=0):
        if drawNode == None:
            drawNode = self.drawNode
        
        vis = []

        # draw root of tree
        vis.append(self.drawNode(conf, node, sx, sy))

        # draw children
        if len(node.children) > 0:
            x = -node.size/2.0
            for child in node.children:
                x += child.size/2.0
                vis.append(self.drawTree(conf, child, drawNode, 
                                         sx+x, sy-conf['height']))
                x += child.size/2.0
        
        return list2group(vis)

    def draw(self, conf, matching, tree, drawNode = None):
        self.installBindings()
        setMix(tree.root, matching.genomes.values()[0])
        self.tree = tree
        
        return self.drawTree(conf, tree.root, drawNode)

    def drawOrthologs(self, conf, tree, orths):
        vis = [conf["color-orths"]]

        for orth in orths:
            if orth[0] in tree.nodes and orth[1] in tree.nodes:
                x1 = tree.nodes[orth[0]].x
                y1 = tree.nodes[orth[0]].y
                x2 = tree.nodes[orth[1]].x
                y2 = tree.nodes[orth[1]].y

                vis.append(lines(vertices(x1, y1, x2, y2)))

        return list2group(vis)

def drawLabels(conf, node):
    def same(bags):
        total = {}
        sum = 0
        for bag in bags:
            if len(bag) == 0:
                return True
            for key in bag:
                sum += 1            
                if not key in total:
                    total[key] = 1
                else:
                    total[key] += 1
        return len(total) > 0 and len(total) < sum

    
    def show(bag):
        print bag

    badwords = ['Probable', 
                'type', 
                'protein', 
                'containing', 
                'transcription',
                'domain',
                'gene']

    def walk(node, bag):
        if node.gene != None and node.isLeaf():
            words = node.gene.description.split()
            words2 = []
            for word in words:
                word = word.replace("[", "")
                word = word.replace("]", "")                
                word = word.replace("(", "")
                word = word.replace(")", "")                
                word = word.replace(".", "")
                word = word.lower()
                words2.append(word)
                
            words = filter(lambda x: 
                not x in badwords and len(x) > 3, words2)
            
            for word in words:
                if not word in bag:
                    bag[word] = 1
                else:
                    bag[word] += 1
            return True
        else:
            cont = True
            bags = []
            for child in node.children:
                bags.append({})
                cont = cont and walk(child, bags[-1])
            
            print "***************"
            print bags
                
            if same(bags) and cont:
                
                for b in bags:
                    for key in b:
                        if not key in bag:
                            bag[key] = b[key]
                        else:
                            bag[key] += b[key]
                print bag
                print "====="
                return True
            else:
                print "-----"
                print node.size
                for i in range(len(bags)):
                    if len(bags[i]) > 0 and \
                       max(bags[i].values()) > 3:
                        show(bags[i]) 
                return False
    walk(node, {})
    
    

# -----------------------------------------------------------------------------


def findOrths(genomes, node):
    orths = []
    
    def helper(node):
        if node.isLeaf():
            return ({node.gene.chrom.genome:True}, [node.gene], False)
        else:
            species = {}
            genes = []
            complete = False
            childGroups = []
            
            for child in node.children:
                (s, g, c) = helper(child)
                
                # merge data
                for i in s:
                    species[i] = True
                
                genes.extend(g)
                
                if c:
                    complete = True
                else:
                    childGroups.append(g)

            # determine if we have all genomes
            if complete:
                # make incomplete childGroups into homology groups
                for group in childGroups:
                    orths.append(group)
                
                # propogate complete
                return ({}, [], True)
            else:
                if len(species) == len(genomes):
                    orths.append(genes)
                    return ({}, [], True)
                else:
                    return (species, genes, False)

    helper(node)
    return orths


def findGenomeSet(genes):
    gset = {}
    for gene in genes:
        gset[gene.chrom.genome] = True
    return gset
    

def isOneToOne(genomes, orth):
    gset = findGenomeSet(orth)        
    return len(gset) == len(genomes) and len(orth) == len(genomes)
        
    
def isIncomplete(genomes, orth):
    gset = findGenomeSet(orth)
    return len(gset) < len(genomes)

def findOneToOnes(genomes, orths):
    return filter(lambda x: isOneToOne(genomes, x), orths)

def findIncompletes(genomes, orths):
    return filter(lambda x: isIncomplete(genomes, x), orths)

if __name__ == "__main__":
    import algorithms
    import genomeio
    
    m = Matching()
    genomeio.readGenomes(m, ["human", "dog"])
    
    tree = algorithms.Tree()
    tree.readNewick("test/small2.tree")
    genes = m.getGenes()
    setupTree(tree, genes)
    
    orths = findOrths(m.genomes.values(), tree.root)
    ones = findOneToOnes(m.genomes.values(), orths)
    incomps = findIncompletes(m.genomes.values(), orths)
