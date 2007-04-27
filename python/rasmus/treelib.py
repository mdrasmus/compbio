#
# Tree data structures 
#
# Contains special features for representing phylogeny.  
# See rasmus.phylo for more.
#
#


# python libs
import copy
import math
import random
import sys
import os
import StringIO

# rasmus libs
from rasmus import textdraw
from rasmus import util

# newick parsing support
import pyparsing



# allow tree reading extra recursion levels
sys.setrecursionlimit(4000)

#============================================================================
# Newick parsing
#


#
# TODO: try D-parser, it might be faster
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
    name    = Word(alphanums + "_" + "-" + "." + "+")
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



#============================================================================
# Tree data structures
#


class TreeNode:
    def __init__(self, name):
        self.name = name
        self.children = []
        self.parent = None
        self.dist = 0
        self.data = {}
    
    
    def __iter__(self):
        """Iterate through child nodes"""
        return iter(self.children)
    
    
    def copy(self, parent=None, copyChildren=True):
        """ Returns a copy of a TreeNode and all of its children"""
        
        node = TreeNode(self.name)
        node.name = self.name
        node.dist = self.dist
        node.parent = parent
        node.data = copy.copy(self.data)
        
        if copyChildren:
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
    
    def leafNames(self):
        return map(lambda x: x.name, self.leaves())
    
    def writeData(self, out):
        out.write(str(self.dist))        


    


class Tree:
    """Basic rooted tree"""

    def __init__(self, nextname = 1):
        self.nodes = {}
        self.root = None
        self.nextname = nextname
        self.defaultData = {}
        self.data = {}
    
    
    def __iter__(self):
        """Iterate through nodes of tree"""
        return self.nodes.itervalues()
        

    def hasData(self, dataname):
        """Does the tree contain 'dataname' in its extra data"""
        return dataname in self.defaultData

    
    def copy(self):
        """Returns a copy of the tree"""
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
        tree.copyData(self)
        tree.copyNodeData(self)
        
        return tree
    
    
    def copyData(self, tree):
        """Copy tree data to another"""
        self.defaultData = copy.copy(tree.defaultData)
        self.data = copy.copy(tree.data)
    
    
    def copyNodeData(self, tree):
        """Copy node data to another tree"""
        for name, node in self.nodes.iteritems():
            if name in tree.nodes:
                node.data = copy.copy(tree.nodes[name].data)
        self.setDefaultData()
    
    
    def setDefaultData(self):
        """Set default values in each node's data"""
        for node in self.nodes.itervalues():
            for key, val in self.defaultData.iteritems():
                node.data.setdefault(key, val)
    
    
    def clearData(self, *keys):
        """Clear tree data"""
        for node in self.nodes.itervalues():
            if len(keys) == 0:
                node.data = {}
            else:
                for key in keys:
                    if key in node.data:
                        del node.data[key]    
    

    def makeRoot(self, name = None):
        """Create a new root node"""
        if name == None:
            name = self.newName()
        self.root = TreeNode(name)
        self.add(self.root)


    def add(self, node):
        """Add a node to the tree
           Does not add node to any specific location (use addChild instead).
        """
        self.nodes[node.name] = node


    def addChild(self, parent, child):
        """Add a child node to an existing node 'parent' in the tree"""
        assert parent != child
        self.nodes[child.name] = child
        self.nodes[parent.name] = parent
        child.parent = parent
        parent.children.append(child)


    def remove(self, node):
        """
        Removes a node from a tree.
        Notifies parent (if it exists) that node has been removed.
        """
        
        if node.parent:
            node.parent.children.remove(node)
        del self.nodes[node.name]
    
    
    def removeTree(self, node):
        """
        Removes subtree rooted at 'node' from tree.
        Notifies parent (if it exists) that node has been removed.
        """
        
        def walk(node):
            if node.name in self.nodes:
                del self.nodes[node.name]
            for child in node.children:
                walk(child)
        walk(node)
        
        if node.parent:
            node.parent.children.remove(node)
    
    
    def rename(self, oldname, newname):
        """Rename a node in the tree"""
        node = self.nodes[oldname]
        del self.nodes[oldname]
        self.nodes[newname] = node
        node.name = newname
    
    
    def newName(self):
        """Returns a new node name that should be unique in the tree"""
        name = self.nextname
        self.nextname += 1
        return name
    
    
    def addTree(self, parent, childTree):
        """Add a subtree to the tree."""
        
        # Merge nodes and change the names of childTree names if they conflict
        # with existing names
        self.mergeNames(childTree)        
        self.addChild(parent, childTree.root)
    
    
    def replaceTree(self, node, childTree):
        """Remove node and replace it with the root of childTree"""
    
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
        """Merge the node names from tree2 into this tree.
           Change any names that conflict"""
    
        for name in tree2.nodes:
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
        """Clear all nodes from tree"""
        self.nodes = {}
        self.root = None
    
    
    def leaves(self, node=None):
        """Return the leaves of the tree in order"""
        if node == None:
            node = self.root                   
        return node.leaves()
    
    
    def leafNames(self, node = None):
        """Returns the leaf names of the tree in order"""
        return map(lambda x: x.name, self.leaves(node))
    
    
    #
    # input and output
    #
    def write(self, out = sys.stdout, writeData=None, oneline=False):
        """Write the tree in newick notation"""
        self.writeNewick(util.openStream(out, "w"), writeData=writeData, 
                         oneline=oneline)        
    
    def readData(self, node, data):
        """Default data reader: reads optional bootstrap and branch length"""
        
        if ":" in data:
            boot, dist = data.split(":")
            node.dist = float(dist)
            
            if len(boot) > 0:
                if boot.isdigit():
                    node.data["boot"] = int(boot)
                else:
                    node.data["boot"] = float(boot)
    
    
    
    def writeData(self, node):
        """Default data writer: writes optional bootstrap and branch length"""
        
        string = ""
        if "boot" in node.data and \
           not node.isLeaf() and \
           self.root != node:
            if isinstance(node.data["boot"], int):
                string += "%d" % node.data["boot"]
            else:
                string += "%f" % node.data["boot"]
        string += ":%f" % node.dist
        return string
    
    
    def writeNewick(self, out = sys.stdout, writeData=None, oneline=False):
        """Write the tree in newick notation"""
        self.writeNewickNode(self.root, util.openStream(out, "w"), 
                             writeData=writeData, oneline=oneline)
   
    
    def writeNewickNode(self, node, out = sys.stdout, 
                              depth = 0, writeData=None, oneline=False):
        # default data writer
        if writeData == None:
            writeData = self.writeData
        
        if not oneline:
            out.write(" " * depth)

        if len(node.children) == 0:
            # leaf
            out.write(node.name)
        else:
            # internal node
            if oneline:
                out.write("(")
            else:
                out.write("(\n")
            for child in node.children[:-1]:
                self.writeNewickNode(child, out, depth+1, 
                                     writeData=writeData, oneline=oneline)
                if oneline:
                    out.write(",")
                else:
                    out.write(",\n")
            self.writeNewickNode(node.children[-1], out, depth+1,
                                 writeData=writeData, oneline=oneline)
            if oneline:
                out.write(")")
            else:            
                out.write("\n" + (" " * depth) + ")")

        # don't print data for root node
        if depth == 0:
            if oneline:
                out.write(";")
            else:
                out.write(";\n")
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
    

    def readParentTree(self, treeFile, labelFile=None):
        if labelFile:
            labels = util.readStrings(labelFile)
        else:
            nitems = (len(file(treeFile).readlines()) + 1)/ 2
            labels = map(str, range(nitems))
        
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
        
        # remove redunant root
        if len(self.root.children) == 1:
            self.root = self.root.children[0]
            self.remove(self.root.parent)
            self.root.parent = None
    
    
    def writeParentTree(self, treeFile, labels):        
        ids = {}
        
        # assign ids to leaves
        for leafname in labels:
            ids[self.nodes[leafname]] = len(ids)
        
        # assign ids to internal nodes
        def walk(node):
            node.recurse(walk)
            if not node.isLeaf():
                ids[node] = len(ids)
        walk(self.root)
        
        # build ptree array
        ptree = [0] * len(ids)
        for node, idname in ids.iteritems():
            if node.parent != None:
                ptree[idname] = ids[node.parent]
            else:
                ptree[idname] = -1
        
        util.writeVector(treeFile, ptree)
    
    
    def getOnelineNewick(self):
        stream = StringIO.StringIO()
        self.write(stream, oneline=True)
        return stream.getvalue()
    
    

    
    
    #
    # should I make these external?
    #
    
    
    def findDepths(self, node = None):
        """DEPRECATED"""
        
        if not node:
            node = self.root
        
        depths = {}

        def walk(node, d):
            depths[node.name] = d
            for child in node.children:
                walk(child, d+1)
        walk(node, 0)
        return depths



#============================================================================
# Misc. functions for manipulating trees
#

def readTree(filename):
    """Read a tree from a newick file"""
    tree = Tree()
    tree.readNewick(filename)
    return tree


def parseNewick(newick):
    """Read a tree from newick notation stored in a string"""
    tree = Tree()
    stream = StringIO.StringIO(newick)
    tree.readNewick(stream)
    return tree


def assertTree(tree):
    """Assert that the tree is internally consistent"""
    
    visited = {}
    def walk(node):
        assert node.name in tree.nodes
        assert node.name not in visited
        visited[node.name] = 1
        if node.parent:
            assert node in node.parent.children
        for child in node.children:
            assert child.parent == node
        node.recurse(walk)
    walk(tree.root)
    
    assert tree.root.parent == None
    assert len(tree.nodes) == len(visited), "%d %d" % (len(tree.nodes), len(visited))




def lca(nodes):
    """Returns the Least Common Ancestor (LCA) of a list of nodes"""
    
    if len(nodes) == 1:
        return nodes[0]
    elif len(nodes) > 2:
        return lca([lca(nodes[:2])] + nodes[2:])
    elif len(nodes) == 2:
        node1, node2 = nodes
        set1 = set([node1])
        set2 = set([node2])
        
        while True:
            if node1 in set2:
                return node1
            if node2 in set1:
                return node2
            if node1.parent != None:
                node1 = node1.parent
            if node2.parent != None:
                node2 = node2.parent
            
            set1.add(node1)
            set2.add(node2)
    else:
        raise Exception("No nodes given")
        

def countDescendents(node, sizes=None):
    if sizes == None:
        sizes = {}
    
    if len(node.children) > 0:
        sizes[node] = 0
        for child in node.children:
            countDescendents(child, sizes)
            sizes[node] += sizes[child]
    else:
        sizes[node] = 1
    
    return sizes


def descendents(node, lst=None):
    if lst == None:
        lst = []
    for child in node.children:
        lst.append(child)
        descendents(child, lst=lst)
    return lst


def subtree(tree, node):
    """Return a copy of a subtree of 'tree' rooted at 'node'"""
    
    # make new tree
    tree2 = Tree(nextname = tree.newName())
    
    # copy nodes and data
    tree2.root = node.copy()    
    tree2.copyData(tree)
    
    # add nodes
    def walk(node):
        tree2.add(node)
        node.recurse(walk)
    walk(tree2.root)
    
    return tree2


def findDist(tree, name1, name2):
    """Returns the branch distance between two nodes in a tree"""

    if not name1 in tree.nodes or \
       not name2 in tree.nodes:
        raise Exception("nodes '%s' and '%s' are not in tree" %
                        (name1, name2))
    
    # find root path for node1
    node1 = tree.nodes[name1]
    path1 = [node1]    
    while node1 != tree.root:
        node1 = node1.parent
        path1.append(node1)
    
    # find root path for node2
    node2 = tree.nodes[name2]
    path2 = [node2]
    while node2 != tree.root:
        node2 = node2.parent
        path2.append(node2)
    
    # find when paths diverge
    i = 1
    while i <= len(path1) and i <= len(path2) and (path1[-i] == path2[-i]):
        i += 1
    
    dist = 0
    for j in range(i, len(path1)+1):
        dist += path1[-j].dist
    for j in range(i, len(path2)+1):
        dist += path2[-j].dist
    
    return dist

"""
def findPath(tree, name1, name2):
    '''
    TODO: double check this code before use
    '''
    if not name1 in tree.nodes or \
       not name2 in tree.nodes:
        return None

    # find root path for node1
    node1 = tree.nodes[name1]        
    path1 = [node1]
    while node1 != tree.root:
        node1 = node1.parent
        path1.append(node1)

    # find root path for node2
    node2 = tree.nodes[name2]        
    path2 = [node2]
    while node2 != tree.root:
        node2 = node2.parent
        path2.append(node2)

    # find when paths diverge
    i = -1
    while i < len(path1) and i < len(path2) and (path1[i] == path2[i]):
        i -= 1

    return path1[i+1:] + path2[i+1:]
"""

def smallSubtrees(tree, maxsize):
    trees = []
    sizes = countDescendents(tree.root)
    
    def walk(node):
        if sizes[node] <= maxsize:
            trees.append(subtree(tree, node))
        else:
            # if too big keep digging
            for child in node.children:
                walk(child)
    walk(tree.root)
    
    return trees


def tree2graph(tree):
    """Convert a tree to a graph data structure (sparse matrix)"""
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
    """Convert a graph to a tree data structure"""
    
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
    """Returns a list of rooted subtrees with atmost one node from 
       the list 'subroots'
    """
    
    marks = {}

    # mark the path from each subroot to the root
    for subroot in subroots:
        ptr = subroot
        while ptr != None:
            lst = marks.setdefault(ptr, [])
            lst.append(subroot)
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
    """
    Remove all nodes from the tree that have exactly one child
    
    Branch lengths are added together when node is removed.
    
    """
    removed = []
    
    # find single children
    def walk(node):
        if len(node.children) == 1 and node.parent:
            removed.append(node)
        node.recurse(walk)
    walk(tree.root)
    
    # actually remove children
    for node in removed:
        newnode = node.children[0]
        
        # add distance
        newnode.dist += node.dist
        
        # change parent and child pointers
        newnode.parent = node.parent
        index = node.parent.children.index(node)
        node.parent.children[index] = newnode
        
        # remove old node
        del tree.nodes[node.name]

    # remove singleton from root
    if len(tree.root.children) == 1:
        oldroot = tree.root
        tree.root = tree.root.children[0]
        oldroot.children = []
        tree.remove(oldroot)
        tree.root.dist += oldroot.dist
    
    return map(lambda x: x.name, removed)


def removeExposedInternalNodes(tree):
    """
    Remove all leaves that were internal nodes (e.g. node.name is an int)
    """
    
    def walk(node):
        children = copy.copy(node.children)
        for child in children:
            walk(child)

        if len(node.children) == 0 and \
           isinstance(node.name, int):
            tree.remove(node)
    walk(tree.root)
    



#=============================================================================
# Rerooting functions
#


def isRooted(tree):
    return len(tree.root.children) <= 2
    #return len(tree.root.children) == 3 or len(tree.leaves()) <= 2


def unroot(tree, newCopy = True):
    """Return an unrooted copy of tree"""
    
    if newCopy:
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


def reroot(tree, newroot, mat=None, onBranch=True, newCopy=True):
    """
    Change the rooting of a tree
    """
    
    # mat is a DEPRECATED parameter.  But if it is used then call old function
    # TODO: remove newCopy (or assert newCopy=False)
    if mat != None:
        return reroot_old(tree, newroot, mat, onBranch=onBranch)
    
    if newCopy:
        tree = tree.copy()
    

    # handle trivial case
    if tree.root.name == newroot or \
       (newroot in [x.name for x in tree.root.children] and \
        len(tree.root.children) == 2):
        return tree        
    
    
    unroot(tree, newCopy=False)
    
    if onBranch:
        # add new root in middle of branch
        newNode = TreeNode(tree.newName())
        node1 = tree.nodes[newroot]
        rootdist = node1.dist
        node1.dist = rootdist / 2.0
        newNode.dist = rootdist / 2.0
        node2 = node1.parent
        node2.children.remove(node1)
        tree.addChild(newNode, node1)
        tree.addChild(node2, newNode)
        
        ptr = node2
        ptr2 = newNode
        newRoot = newNode
    else:
        # root directly on root
        ptr2 = tree.nodes[newroot]
        ptr = ptr2.parent
        newRoot = ptr2
    
    newRoot.parent = None
    
    # reverse parent child relationship of all nodes on path node1 to root
    oldroot = tree.root    
    nextDist = ptr2.dist
    ptr2.dist = 0
    while True:
        nextPtr = ptr.parent
        ptr.children.remove(ptr2)
        tree.addChild(ptr2, ptr)
        
        tmp = ptr.dist
        ptr.dist = nextDist
        nextDist = tmp
        
        ptr2 = ptr
        ptr = nextPtr
        
        if nextPtr is None:
            break
    tree.root = newRoot
    
    return tree
    


def reroot_old(tree, newroot, mat = None, onBranch=True):
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
                childNode.data = copy.copy(tree.nodes[child].data)
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
    


#=========================================================================
# Draw Tree ASCII art 
#

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
            canvas.text(xchildren +1, y+yscale-1, str(node.name))
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
            if isinstance(node.data["boot"], int):
                labels[node.name] = "(%d) %f" % (node.data["boot"], node.dist)
            else:
                labels[node.name] = "(%.2f) %f" % (node.data["boot"], node.dist)
    
    drawTree(tree, labels, *args, **kargs)


def drawTreeNames(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.isLeaf():
            labels[node.name] = "%s" % node.name
    
    drawTree(tree, labels, *args, **kargs)


def drawTreeNameLens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.isLeaf():
            labels[node.name] = "%s " % node.name
        else:
            labels[node.name] =""
        labels[node.name] += "%f" % node.dist
    
    drawTree(tree, labels, *args, **kargs)
