# python libs
import copy
import math
import random
import sys

# rasmus libs
import graph
import textdraw
import util

import pyparsing

sys.setrecursionlimit(4000)


#
# Newick parsing
#

def makeNewickParser():
    # pyparsing
    from pyparsing import Combine, Optional, Literal, CaselessLiteral, \
                       Word, alphanums, \
                       nums, oneOf, Group, Dict, Forward, \
                       ParseResults, CharsNotIn, ZeroOrMore


    # literals
    lparen    = Literal("(").suppress()
    rparen    = Literal(")").suppress()
    colon     = Literal(":").suppress()
    semicolon = Literal(":").suppress()
    comma     = Literal(",").suppress()
    point     = Literal(".")
    e         = CaselessLiteral("E")


    # terminal rules
    name    = Word(alphanums + "_" + "-" + ".")
    fnumber = Combine(Word("+-"+nums, nums) + 
                      Optional(point + Optional(Word(nums))) +
                      Optional(e + Word("+-"+nums, nums)))
    dist      = fnumber
    bootstrap = fnumber
    
    
    # recursive rules
    subtree = Forward()
    subtreelist = Forward()

    subtree << \
        Group(
            (
                (lparen + subtreelist + rparen).setResultsName("subtree") |
                name.setResultsName("name")
            ) +
            Optional(
                CharsNotIn(",);").setResultsName("data")
            )
        )
    subtreelist << subtree + Optional(comma + subtreelist)
    
    
    # top level rule
    tree = subtree + Word(";").suppress()
    
    
    return tree.parseString

newickParser = makeNewickParser()




class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.parent = None
        self.dist = 0
        self.data = {}
    
    
    def copy(self, parent=None):
        """ Returns a copy of a TreeNode and all of its children"""
        
        node = TreeNode(self.name)
        node.name = self.name
        node.dist = self.dist
        node.parent = parent
        node.data = copy.copy(self.data)
        
        for child in self.children:
            node.children.append(child.copy(node))
        
        return node
    
    
    def isLeaf(self):
        return len(self.children) == 0
    
    def recurse(self, func, *args):
        for child in self.children:
            func(child, *args)
    
    def leaves(self):
        leaves = []

        def walk(node):
            if node.isLeaf():
                leaves.append(node)
            node.recurse(walk)
        walk(self)
          
        return leaves

    def leaveNames(self):
        return map(lambda x: x.name, self.leaves())
    
    
    def writeData(self, out):
        out.write(str(self.dist))        


    


class Tree:
    def __init__(self, nextname = 1):
        self.nodes = {}
        self.root = None
        self.nextname = nextname
        self.defaultData = {}
        self.data = {}
    

    def hasData(self, dataname):
        return dataname in self.defaultData

    
    def copy(self):
        tree = Tree(nextname = self.nextname)
        if self.root != None:
            # copy all nodes
            tree.root = self.root.copy()
            
            # set all names
            def walk(node):
                tree.nodes[node.name] = node
                for child in node.children:
                    walk(child)
            walk(tree.root)
        
        # copy extra data
        tree.defaultData = copy.copy(self.defaultData)
        tree.data = copy.copy(self.data)
        tree.copyData(self)
        
        return tree
    
    
    def copyData(self, tree):
        self.defaultData = copy.copy(tree.defaultData)
        for name, node in self.nodes.iteritems():
            if name in tree.nodes:
                node.data = copy.copy(tree.nodes[name].data)
        self.setDefaultData()
    
    
    def setDefaultData(self):
        for node in self.nodes.itervalues():
            for key, val in self.defaultData.iteritems():
                node.data.setdefault(key, val)
    
    
    def clearData(self, *keys):
        for node in self.nodes.itervalues():
            if len(keys) == 0:
                node.data = {}
            else:
                for key in keys:
                    if key in node.data:
                        del node.data[key]    
    

    def makeRoot(self, name = None):
        if name == None:
            name = self.newName()
        self.root = TreeNode(name)
        self.add(self.root)


    def add(self, node):
        self.nodes[node.name] = node


    def addChild(self, parent, child):
        assert parent != child
        self.nodes[child.name] = child
        self.nodes[parent.name] = parent
        child.parent = parent
        parent.children.append(child)


    def remove(self, node):
        if node.parent:
            node.parent.children.remove(node)
        del self.nodes[node.name]
    
    
    def removeTree(self, node):
        def walk(node):
            if node.name in self.nodes:
                del self.nodes[node.name]
            for child in node.children:
                walk(child)
        walk(node)
        
        if node.parent:
            node.parent.children.remove(node)
    
    
    def rename(self, oldname, newname):
        node = self.nodes[oldname]
        del self.nodes[oldname]
        self.nodes[newname] = node
        node.name = newname
    
    
    def newName(self):
        name = self.nextname
        self.nextname += 1
        return name
    
    
    def addTree(self, parent, childTree):
        # merge nodes and change the names of childTree names if they conflict
        # with existing names
        self.mergeNames(childTree)        
        self.addChild(parent, childTree.root)
    
    
    def replaceTree(self, node, childTree):
        # merge nodes and change the names of childTree names if they conflict
        # with existing names
        self.mergeNames(childTree)
        parent = node.parent
        if parent:
            index = parent.children.index(node)
            parent.children[index] = childTree.root
            childTree.root.parent = parent
            del self.nodes[node.name]
    
    
    def mergeNames(self, tree2):
        names = tree2.nodes.keys()
        for name in names:
            if name in self.nodes:
                name2 = self.newName()
                self.nodes[name2] = tree2.nodes[name]
                self.nodes[name2].name = name2
            else:
                if type(name) == int:
                    if name >= self.nextname:
                        self.nextname = name + 1
                self.nodes[name] = tree2.nodes[name]
        
    
    def clear(self):
        self.nodes = {}
        self.root = None
    
    
    def leaves(self, node = None):
        if node == None:
            node = self.root                   
        return node.leaves()
    
    
    def leaveNames(self, node = None):
        return map(lambda x: x.name, self.leaves(node))



    #
    # input and output
    #
    def readData(self, node, data):
        """Default data reader: reads optional bootstrap and branch length"""
        
        if ":" in data:
            boot, dist = data.split(":")
            node.dist = float(dist)
            
            if boot.isalnum():
                node.data["boot"] = int(boot)
    
    def writeData(self, node):
        """Default data writer: writes optional bootstrap and branch length"""
        
        string = ""
        if "boot" in node.data and \
           not node.isLeaf() and \
           self.root != node:
            string += "%d" % node.data["boot"]
        string += ":%f" % node.dist
        return string
    
    
    def writeNewick(self, out = sys.stdout, writeData=None):
        self.writeNewickNode(self.root, util.openStream(out, "w"), 
                             writeData=writeData)
    
    
    def write(self, out = sys.stdout, writeData=None):
        self.writeNewick(util.openStream(out, "w"), writeData=writeData)

    
    def writeNewickNode(self, node, out = sys.stdout, depth = 0, writeData=None):
        # default data writer
        if writeData == None:
            writeData = self.writeData
        
        print >>out, (" " * depth),

        if len(node.children) == 0:
            # leaf
            print >>out, node.name,
        else:
            # internal node
            print >>out, "("
            for child in node.children[:-1]:
                self.writeNewickNode(child, out, depth+1)
                print >>out, ","
            self.writeNewickNode(node.children[-1], out, depth+1)            
            print >>out
            print >>out, (" " * depth),
            print >>out, ")",

        # don't print data for root node
        if depth == 0:
            print >>out, ";"
        else:
            out.write(writeData(node))
    
    
    def readNewick(self, filename, readData=None):
        try:

            # default data reader
            if readData == None:
                readData = self.readData

            # get parse tree
            text = util.readUntil(util.openStream(filename), ";")[0] + ";"
            expr = newickParser(text)[0]


            # walk the parse tree and build the tree
            self.clear()

            def walk(expr):
                if isinstance(expr, pyparsing.ParseResults):
                    # parse name
                    if "name" in expr:
                        node = TreeNode(expr["name"])
                    else:
                        node = TreeNode(self.newName())

                    if "data" in expr:
                        readData(node, expr["data"])

                    # recurse
                    for child in expr:
                        ret = walk(child)
                        if ret:
                            self.addChild(node, ret)
                    return node

            self.root = walk(expr)

            # test for boot strap presence
            for node in self.nodes.itervalues():
                if "boot" in node.data:
                    self.defaultData["boot"] = 0
                    break
            self.setDefaultData()
        
        except RuntimeError:
            # try simpler parser
            return self.readBigNewick(filename)
    
    
    def readBigNewick(self, filename):
        infile = file(filename)    
        closure = {"opens": 0}

        def readchar():
            while True:
                char = infile.read(1)
                if char != " " and char != "\n": break
            if char == "(": closure["opens"] += 1
            if char == ")": closure["opens"] -= 1
            return char
        
        def readUntil(chars):
            token = ""
            while True:
                char = readchar()
                if char in chars or char == "":
                    return [token, char]
                token += char
        
        def readDist():
            word = ""
            while True:
                char = readchar()
                if not char in "-0123456789.e":
                    return float(word)
                else:
                    word += char

        def readItem():
            try:
                char1 = readchar()

                if char1 == "(":
                    node = TreeNode(self.newName())
                    depth = closure["opens"]
                    while closure["opens"] == depth:
                        self.addChild(node, readItem())
                    
                    token, char = readUntil("):,")
                    if char == ":":
                        node.dist = readDist()
                    return node
                else:                   
                    word, char = readUntil(":),")
                    word = char1 + word.rstrip()
                    node = TreeNode(word)
                    if char == ":":
                        node.dist = readDist()
                    return node
            except:
                print sys.exc_type, ": Tree too deep to read"
                return TreeNode("TOO_DEEP")
        

        def readRoot():
            word, char = readUntil("(")
            
            assert char == "("
            
            node = TreeNode(self.newName())
            depth = closure["opens"]
            while closure["opens"] == depth:
                self.addChild(node, readItem())
            return node

        self.root = readRoot()
        self.add(self.root)
    

    def readParentTree(self, treeFile, labelFile):
        labels = util.readStrings(labelFile)
        self.makeRoot()

        i = 0
        for line in file(treeFile):
            parentid = int(line.split(" ")[0])
            
            # determine current child
            if i < len(labels):
                child = TreeNode(labels[i])
            else:
                if i in self.nodes:
                    child = self.nodes[i]
                else:
                    child = TreeNode(i)
            
            if parentid == -1:
                # keep track of all roots
                self.addChild(self.root, child)
            else:                
                if not parentid in self.nodes:
                    parent = TreeNode(parentid)
                    self.add(parent)
                else:
                    parent = self.nodes[parentid]

                try:
                    self.addChild(parent, child)
                except:
                    print i, parentid
            i += 1
    
    
    #
    # should I make these external?
    #
    
    def setSizes(self, node = None):
        if node == None:
            node = self.root
        if len(node.children) > 0:
            node.size = 0
            for child in node.children:
                node.size += self.setSizes(child)
        else:
            node.size = 1
        return node.size
    
    
    def findDepths(self, node = None):
        if not node:
            node = self.root
        
        depths = {}

        def walk(node, d):
            depths[node.name] = d
            for child in node.children:
                walk(child, d+1)
        walk(node, 0)
        return depths


    def findDist(self, name1, name2):
        if not name1 in self.nodes or \
           not name2 in self.nodes:
            return None
    
        # find root path for node1
        path1 = []
        node1 = self.nodes[name1]
        while node1 != self.root:
            node1 = node1.parent
            path1.append(node1)
        
        # find root path for node2
        path2 = []
        node2 = self.nodes[name2]
        while node2 != self.root:
            node2 = node2.parent
            path2.append(node2)
        
        # find when paths diverge
        i = 0
        while i < len(path1) and i < len(path2) and (path1[i] == path2[i]):
            i += 1
        
        return len(path1) + len(path2) - 2 * i + 1
    
    
    def findPath(self, name1, name2):
        if not name1 in self.nodes or \
           not name2 in self.nodes:
            return None
    
        # find root path for node1
        node1 = self.nodes[name1]        
        path1 = [node1]
        while node1 != self.root:
            node1 = node1.parent
            path1.append(node1)
        
        # find root path for node2
        node2 = self.nodes[name2]        
        path2 = [node2]
        while node2 != self.root:
            node2 = node2.parent
            path2.append(node2)
        
        # find when paths diverge
        i = -1
        while i < len(path1) and i < len(path2) and (path1[i] == path2[i]):
            i -= 1
        
        return path1[i+1:] + path2[i+1:]
    

    def lca(self, names, depths = None):
        """Least Common Ancestor"""
        
        if not depths:
            depths = self.findDepths(self.root)

        markers = util.list2dict(names)

        while len(markers) > 1:
            names = markers.keys()
            ind = util.argmaxfunc(lambda x: depths[x], names)
            deepest = names[ind]

            del markers[deepest]
            markers[self.nodes[deepest].parent.name] = 1

        return self.nodes[markers.keys()[0]]


    
    def subtree(self, node):
        tree = Tree(nextname = self.newName())
        tree.root = node
        
        def walk(node):
            tree.add(node)
            for child in node.children:
                walk(child)

        walk(node)
        return tree
    
    



def readTree(filename):
    tree = Tree()
    tree.readNewick(filename)
    return tree


def descendents(node, lst=None):
    if lst == None:
        lst = []
    for child in node.children:
        lst.append(child)
        descendents(child)
    return lst


def smallSubtrees(tree, maxsize):
    trees = []
    tree.setSizes()
    
    def walk(node):
        if node.size <= maxsize:
            trees.append(tree.subtree(node))
        else:
            # if too big keep digging
            for child in node.children:
                walk(child)
    walk(tree.root)
    
    return trees


def tree2graph(tree):
    mat = {}
    
    # init all rows of adjacency matrix to 
    for name, node in tree.nodes.items():
        mat[name] = {}
    
    for name, node in tree.nodes.items():
        for child in node.children:
            mat[name][child.name] = child.dist
        
        if node.parent:
            mat[name][node.parent.name] = node.dist
            
    return mat


def graph2tree(mat, root, closedset=None):
    if closedset == None:
        closedset = {}
    tree = Tree()

    def walk(name):
        node = TreeNode(name)
        node.dist = 0
        closedset[name] = 1
        for child in mat[name]:
            if child not in closedset:
                childNode = walk(child)
                childNode.dist = mat[name][child]
                tree.addChild(node, childNode)
        return node            
    tree.root = walk(root)
    
    tree.nextname = max(filter(lambda x: type(x) == int, tree.nodes.keys()))
    
    return tree


def maxDisjointSubtrees(tree, subroots):
    marks = {}

    # mark the path from each subroot to the root
    for subroot in subroots:
        ptr = subroot
        while True:
            lst = marks.setdefault(ptr, [])
            lst.append(subroot)
            if not ptr.parent: break
            ptr = ptr.parent

    # subtrees are those trees with nodes that have at most one mark
    subroots2 = []
    def walk(node):
        marks.setdefault(node, [])
        if len(marks[node]) < 2 and \
           (not node.parent or len(marks[node.parent]) >= 2):
            subroots2.append(node)
        node.recurse(walk)
    walk(tree.root)
    
    return subroots2


def removeSingleChildren(tree):
    removed = []
    
    # find single children
    def walk(node):
        if len(node.children) == 1 and node.parent:
            removed.append(node)
        node.recurse(walk)
    walk(tree.root)
    
    # actually remove children
    for node in removed:
        dist = node.dist        
        subtree = tree.subtree(node.children[0])
        node.children[0].dist += node.dist
        tree.removeTree(node.children[0])
        tree.replaceTree(node, subtree)

    # remove singleton from root
    if len(tree.root.children) == 1:
        oldroot = tree.root
        tree.root = tree.root.children[0]
        oldroot.children = []
        tree.remove(oldroot)
        tree.root.dist += oldroot.dist
    
    return map(lambda x: x.name, removed)


def unroot(tree):
    """Return an unrooted copy of tree"""
    
    tree = tree.copy()
    if len(tree.root.children) == 2:
        nodes = tree.root.children
        dist = nodes[0].dist + nodes[1].dist
        if len(nodes[0].children) < 2:
            nodes.reverse()
        tree.addChild(nodes[0], nodes[1])
        nodes[1].dist = dist
        nodes[0].dist = 0
        nodes[0].parent = None
        
        # replace root
        del tree.nodes[tree.root.name]
        tree.root = nodes[0]
    return tree


def isRooted(tree):
    return len(tree.root.children) <= 2
    #return len(tree.root.children) == 3 or len(tree.leaves()) <= 2


def assertTree(tree):
    visited = {}
    def walk(node):
        assert node.name in tree.nodes
        assert node.name not in visited
        visited[node.name] = 1
        if node.parent:
            node in node.parent.children
        for child in node.children:
            child.parent == node
        node.recurse(walk)
    walk(tree.root)
    
    assert len(tree.nodes) == len(visited)
    

def reroot(tree, newroot, mat = None, onBranch=True):
    # handle trivial case
    if tree.root.name == newroot or \
       (newroot in map(lambda x: x.name, tree.root.children) and \
        len(tree.root.children) == 2):
        return tree        
    
    # convert tree to a graph
    if mat == None:
        mat = tree2graph(tree)
        assert len(mat.keys()) == len(tree.nodes.keys())
        for i in mat:
            for j in mat[i]:
                assert mat[j][i] == mat[i][j]
        #print mat, newroot
        
    # intialize new tree
    tree2 = Tree(nextname = tree.newName())
    
    # generate new node if rooting on a branch
    # branch designated as branch above newroot
    if onBranch:
        name1 = newroot
        name2 = tree.nodes[newroot].parent.name
        dist = mat[name1][name2]
        del mat[name1][name2]
        del mat[name2][name1]
        
        # create new name for newroot
        newroot = tree2.newName()
        mat[newroot] = {}
        mat[newroot][name1] = dist / 2.0
        mat[name1][newroot] = dist / 2.0
        mat[newroot][name2] = dist / 2.0
        mat[name2][newroot] = dist / 2.0
        
        #print name1, name2, newroot
    
    # build new tree   
    tree2.makeRoot(newroot)
    closedset = {newroot:1}
    def walk(node):
        for child in mat[node]:
            if not child in closedset:
                childNode = TreeNode(child)
                childNode.dist = mat[child][node]
                tree2.addChild(tree2.nodes[node], childNode)
                closedset[child] = 1
                walk(child)
    walk(newroot)
    
    # copy over data
    # TODO: need to fix how extra data is copied, should be copied with dist above
    # branches have changed identies
    tree2.copyData(tree)
    
    assert len(tree2.nodes.keys()) >= len(tree.nodes.keys())
    
    # clean up tree and update mat
    removed = removeSingleChildren(tree2)
    
    # undo changes to mat
    if onBranch:
        mat[name1][name2] = dist
        mat[name2][name1] = dist
        del mat[newroot]
        del mat[name1][newroot]
        del mat[name2][newroot]
    
    assert len(tree2.nodes.keys()) >= len(tree.nodes.keys())
    
    return tree2



def outgroupRoot(tree, outgroup, closedset = None):
    # find root of outgroup
    mat = tree2graph(tree)
    neighbors = []
    if closedset == None:
        closedset = {}
    roots = {}
    
    # remove a vertex from a graph
    def remove(mat, v):
        for u in mat[v]:
            del mat[u][v]
        del mat[v]
    
    # start a special bfs
    openset = outgroup
    while len(openset) > 0:
        vertex = openset[0]
        openset = openset[1:]
        
        # skip closed vertices
        if vertex in closedset:
            continue

        # visit vertex
        if len(mat[vertex]) == 1:
            # add neighbors to openset
            openset.extend(mat[vertex].keys())
            
            # close and remove this vertex
            closedset[vertex] = 1
            remove(mat, vertex)
            if vertex in roots:
                del roots[vertex]
        else:
            roots[vertex] = 1
        
    return roots.keys()


def removeOutgroup(tree, outgroup):
    removed = {}
    roots = outgroupRoot(tree, outgroup, removed)

    if len(roots) == 1:
        tree2 = reroot(tree, roots[0])
    else:
        tree2 = tree
        
    for root in roots:
        for child in tree2.nodes[root].children:
            if child.name in removed:
                tree2.removeTree(child)
    
    return (tree2, len(roots) == 1)
    


def drawTree(tree, labels={}, scale=40, spacing=2, out=sys.stdout,
             canvas=None, x=0, y=0, display=True, labelOffset=-1,
             minlen=1,maxlen=10000):
    if canvas == None:
        canvas = textdraw.TextCanvas()
    
    xscale = scale
    yscale = spacing

    
    # determine node sizes
    sizes = {}
    nodept = {}
    def walk(node):
        if node.isLeaf():
            sizes[node] = 1
            nodept[node] = yscale - 1 
        else:
            sizes[node] = 0
        for child in node.children:
            sizes[node] += walk(child)
        if not node.isLeaf():
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2
        return sizes[node]
    walk(tree.root)
    
    
    def walk(node, x, y):
        # calc coords
        xchildren = int(x+min(max(node.dist*xscale,minlen),maxlen))
        
        # draw branch
        canvas.line(x, y+nodept[node], xchildren, y+nodept[node], '-')
        if node.name in labels:
            branchlen = xchildren - x
            lines = str(labels[node.name]).split("\n")
            labelwidth = max(map(len, lines))
            
            labellen = min(labelwidth, 
                           max(int(branchlen-1),0))
            canvas.text(x + 1 + (branchlen - labellen)/2., 
                        y+nodept[node]+labelOffset, 
                        labels[node.name], width=labellen)
        
        if node.isLeaf():
            canvas.text(xchildren +1, y+yscale-1, node.name)
        else:
            top = y + nodept[node.children[0]]
            bot = y + (sizes[node]-sizes[node.children[-1]]) * yscale + \
                      nodept[node.children[-1]]
        
            # draw children
            canvas.line(xchildren, top, xchildren, bot, '|')
            
            ychild = y
            for child in node.children:
                walk(child, xchildren, ychild)
                ychild += sizes[child] * yscale

            
            canvas.set(xchildren, y+nodept[node], '+')
            canvas.set(xchildren, top, '/')
            canvas.set(xchildren, bot, '\\')
        canvas.set(x, y+nodept[node], '+')
    walk(tree.root, x+0, 0)
    
    if display:
        canvas.display(out)


def drawTreeLens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        labels[node.name] = "%f" % node.dist
    
    drawTree(tree, labels, *args, **kargs)


def drawTreeBootLens(tree, *args, **kargs):
    if not tree.hasData("boot"):
        drawTreeLens(tree, *args, **kargs)
        return

    labels = {}
    for node in tree.nodes.values():
        if node.isLeaf():
            labels[node.name] = "%f" % node.dist
        else:
            labels[node.name] = "(%d) %f" % (node.data["boot"], node.dist)
    
    drawTree(tree, labels, *args, **kargs)


def drawTreeNames(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.isLeaf():
            labels[node.name] = "%s" % node.name
    
    drawTree(tree, labels, *args, **kargs)

