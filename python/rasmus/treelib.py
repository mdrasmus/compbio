#
# Tree data structures 
#
# Contains special features for representing phylogeny.  
# See rasmus.bio.phylo for more.
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
from rasmus import util


# ply parsing support
from rasmus import treelib_parser
    
try:
    from rasmus import treelib_parser
except ImportError:
    try:
        import treelib_parser
    except ImportError:
        treelib_parser = None




# allow tree reading extra recursion levels
sys.setrecursionlimit(4000)



#============================================================================
# Tree data structures
#


class TreeNode:
    """A class for nodes in a rooted Tree
    
    Contains fields for branch length 'dist' and custom data 'data'
    """

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
        """Returns a copy of a TreeNode and all of its children"""
        
        node = TreeNode(self.name)
        node.name = self.name
        node.dist = self.dist
        node.parent = parent
        node.data = copy.copy(self.data)
        
        if copyChildren:
            for child in self.children:
                node.children.append(child.copy(node))
        
        return node
        
    def is_leaf(self):
        """Returns True if the node is a leaf (no children)"""
        return len(self.children) == 0
    isLeaf = is_leaf
    
    def recurse(self, func, *args):
        """Applies a function 'func' to the children of a node"""
        for child in self.children:
            func(child, *args)
    
    def leaves(self):
        """Returns the leaves beneath the node in traversal order"""
        leaves = []

        def walk(node):
            if node.is_leaf():
                leaves.append(node)
            for child in node.children:
                walk(child)
        walk(self)
          
        return leaves
    
    def leaf_names(self):
        """Returns the leaf names beneath the node in traversal order"""
        return [x.name for x in self.leaves()]
    leafNames = leaf_names
    
    def write_data(self, out):
        """Writes the data of the node to the file stream 'out'"""
        out.write(str(self.dist))
    writeData = write_data


# class for managing branch data
class BranchData (object):
    """A class for managing branch specific data for a Tree
    
    By default, this class implements bootstrap data for TreeNode's.  
    
    To incorporate new kinds of branch data, do the following.  Subclass this
    class (say, MyBranchData).  Create Tree's with
    Tree(branch_data=MyBranchData()).  This will ensure your new branch data
    manager is used for manipulations to the tree.  Any tree's copied from the
    tree (via tree.copy()) will also use the same branch manager.
    """

    def __init__(self):
        pass

    def get_branch_data(self, node):
        """Returns branch specific data from a node"""
        if "boot" in node.data:
            return {"boot": node.data["boot"]}
        else:
            return {}
    
    def set_branch_data(self, node, data):
        """Set the branch specific data from 'data' to node.data"""
        if "boot" in data:
            node.data["boot"] = data["boot"]
    
    def split_branch_data(self, node):
        """Split a branch's data into two copies"""
        if "boot" in node.data:
            return {"boot": node.data["boot"]}, {"boot": node.data["boot"]}
        else:
            return {}, {}
    
    def merge_branch_data(self, data1, data2):
        """Merges the branch data from two neighboring branches into one"""
        if "boot" in data1 and "boot" in data2:
            assert data1["boot"] == data2["boot"], (data1, data2)
            return {"boot": data1["boot"]}
        else:
            return {}
    
    

class Tree:
    """
    Basic rooted tree
    
    Well suited for phylogenetic trees
    """

    def __init__(self, nextname = 1, branch_data=BranchData()):
        self.nodes = {}
        self.root = None
        self.nextname = nextname
        self.defaultData = {}
        self.data = {}
        self.branch_data = branch_data

    #=========================================
    # iterators
    
    def __iter__(self):
        """Iterate through nodes of tree"""
        return self.nodes.itervalues()


    def preorder(self, node=None):
        """Iterate through nodes in pre-order traversal"""
        
        if node is None:
            node = self.root

        queue = [node]

        while len(queue) > 0:
            node = queue.pop()
            yield node

            for child in reversed(node.children):
                queue.append(child)

                
    def postorder(self, node=None):
        """Iterate through nodes in post-order traversal"""
        
        if node is None:
            node = self.root

        stack = [[node, 0]]

        while len(stack) > 0:
            node, i = stack[-1]

            if i < len(node.children):
                stack.append([node.children[i], 0])
                stack[-2][1] += 1
            else:
                yield node
                stack.pop()
    
  
    #=============================
    # structure functions

    def make_root(self, name = None):
        """Create a new root node"""
        if name == None:
            name = self.newName()
        self.root = TreeNode(name)
        self.add(self.root)
    makeRoot = make_root


    def add(self, node):
        """Add a node to the tree
           Does not add node to any specific location (use addChild instead).
        """
        self.nodes[node.name] = node


    def add_child(self, parent, child):
        """Add a child node to an existing node 'parent' in the tree"""
        assert parent != child
        self.nodes[child.name] = child
        self.nodes[parent.name] = parent
        child.parent = parent
        parent.children.append(child)
    addChild = add_child


    def remove(self, node):
        """
        Removes a node from a tree.
        Notifies parent (if it exists) that node has been removed.
        """
        
        if node.parent:
            node.parent.children.remove(node)
        del self.nodes[node.name]
    
    
    def remove_tree(self, node):
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
    removeTree = remove_tree

    
    
    def rename(self, oldname, newname):
        """Rename a node in the tree"""
        node = self.nodes[oldname]
        del self.nodes[oldname]
        self.nodes[newname] = node
        node.name = newname
    
    
    def new_name(self):
        """Returns a new node name that should be unique in the tree"""
        name = self.nextname
        self.nextname += 1
        return name
    newName = new_name
    
    
    def add_tree(self, parent, childTree):
        """Add a subtree to the tree."""
        
        # Merge nodes and change the names of childTree names if they conflict
        # with existing names
        self.mergeNames(childTree)        
        self.addChild(parent, childTree.root)
    addTree = add_tree
    
    
    def replace_tree(self, node, childTree):
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
    replaceTree = replace_tree
    
    
    def merge_names(self, tree2):
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
    mergeNames = merge_names
        
    
    def clear(self):
        """Clear all nodes from tree"""
        self.nodes = {}
        self.root = None
    
    
    def leaves(self, node=None):
        """Return the leaves of the tree in order"""
        if node == None:
            node = self.root                   
        return node.leaves()
    
    
    def leaf_names(self, node = None):
        """Returns the leaf names of the tree in order"""
        return map(lambda x: x.name, self.leaves(node))
    leafNames = leaf_names
    
    
    #===============================
    # data functions

    def has_data(self, dataname):
        """Does the tree contain 'dataname' in its extra data"""
        return dataname in self.defaultData
    hasData = has_data

    
    def copy(self):
        """Returns a copy of the tree"""
        tree = Tree(nextname = self.nextname)
        
        # copy structure
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
    
    
    def copy_data(self, tree):
        """Copy tree data to another"""
        self.branch_data = tree.branch_data
        self.defaultData = copy.copy(tree.defaultData)
        self.data = copy.copy(tree.data)
    copyData = copy_data
    
    
    def copy_node_data(self, tree):
        """Copy node data to another tree"""
        for name, node in self.nodes.iteritems():
            if name in tree.nodes:
                node.data = copy.copy(tree.nodes[name].data)
        self.setDefaultData()
    copyNodeData = copy_node_data
    
    def set_default_data(self):
        """Set default values in each node's data"""
        for node in self.nodes.itervalues():
            for key, val in self.defaultData.iteritems():
                node.data.setdefault(key, val)
    setDefaultData = set_default_data
    
    
    def clear_data(self, *keys):
        """Clear tree data"""
        for node in self.nodes.itervalues():
            if len(keys) == 0:
                node.data = {}
            else:
                for key in keys:
                    if key in node.data:
                        del node.data[key]
    clearData = clear_data

    
    
    #======================================================================
    # branch data functions
    # forward branch data calles to branch data manager        
    
    def get_branch_data(self, node):
        """Returns branch specific data from a node"""
        return self.branch_data.get_branch_data(node)

    def set_branch_data(self, node, data):
        """Set the branch specific data from 'data' to node.data"""
        return self.branch_data.set_branch_data(node, data)
    
    def split_branch_data(self, node):
        """Split a branch's data into two copies"""
        return self.branch_data.split_branch_data(node)
    
    def merge_branch_data(self, data1, data2):
        """Merges the branch data from two neighboring branches into one"""    
        return self.branch_data.merge_branch_data(data1, data2)
    
    
    #=======================================================================
    # input and output
    #
    
    def write(self, out = sys.stdout, writeData=None, oneline=False,
              rootData=False):
        """Write the tree in newick notation"""
        self.writeNewick(util.open_stream(out, "w"), writeData=writeData, 
                         oneline=oneline, rootData=rootData)
    
    def read_data(self, node, data):
        """Default data reader: reads optional bootstrap and branch length"""
        
        if ":" in data:
            boot, dist = data.split(":")
            node.dist = float(dist)
            
            if len(boot) > 0:
                if boot.isdigit():
                    node.data["boot"] = int(boot)
                else:
                    try:
                        node.data["boot"] = float(boot)
                    except ValueError:
                        # treat as node name
                        name = boot.strip()
                        if name and not node.is_leaf():
                            node.name = name
        else:
            data = data.strip()
            
            # treat as name
            if data and not node.is_leaf():
                node.name = data
    readData = read_data
    
    
    def write_data(self, node):
        """Default data writer: writes optional bootstrap and branch length"""
        
        string = ""
        if "boot" in node.data and \
           not node.is_leaf() and \
           self.root != node:
            if isinstance(node.data["boot"], int):
                string += "%d" % node.data["boot"]
            else:
                string += "%f" % node.data["boot"]
        else:
            # see if internal node names exist
            if not node.is_leaf() and isinstance(node.name, str):
                string += node.name

        string += ":%f" % node.dist
        return string
    writeData = write_data
    
    
    def write_newick(self, out = sys.stdout, writeData=None, oneline=False,
                    rootData=False):
        """Write the tree in newick notation"""
        self.write_newick_node(self.root, util.open_stream(out, "w"), 
                               writeData=writeData, oneline=oneline,
                               rootData=rootData)
    writeNewick = write_newick
   
    
    def write_newick_node(self, node, out = sys.stdout, 
                          depth=0, writeData=None, oneline=False,
                          rootData=False):
        """Write the node in newick format to the out file stream"""
        
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
            if rootData:
                out.write(writeData(node))
            if oneline:
                out.write(";")
            else:
                out.write(";\n")
        else:
            out.write(writeData(node))
    writeNewickNode = write_newick_node
    
    
    def read_newick(self, filename, readData=None):
        """
        Reads newick tree format from a file stream
        
        You can specify a specialized node data reader with 'readData'
        """
    
        # use simple parsing if PLY is not available

        try:
            if treelib_parser is not None:
                return self._read_newick(filename, readData)
        except RuntimeError:
            pass

        # try simpler parser
        return self.readBigNewick(filename)
    readNewick = read_newick


    def _read_newick(self, filename, readData=None):
        """read with PLY"""
        
        # default data reader
        if readData == None:
            readData = self.readData

        # get parse tree
        text = util.read_until(util.open_stream(filename), ";")[0] + ";"
        print text
        expr = treelib_parser.yacc.parse(text)


        # walk the parse tree and build the tree
        self.clear()

        def walk(expr):
            children, name, data = expr
            
            # parse name
            if name == "":
                node = TreeNode(self.newName())
            else:
                node = TreeNode(name)

            # parse data
            readData(node, data)

            # recurse
            for child in children:
                ret = walk(child)
                if ret:
                    self.add_child(node, ret)
            return node

        self.root = walk(expr)

        # test for boot strap presence
        for node in self.nodes.itervalues():
            if "boot" in node.data:
                self.defaultData["boot"] = 0
                break
        self.setDefaultData()
    
    
    def read_big_newick(self, filename):
        """Reads a big newick file with a custom parser"""
    
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
    readBigNewick = read_big_newick
    

    def read_parent_tree(self, treeFile, labelFile=None, labels=None):
        """Reads a parent array from a file"""
    
        if labelFile:
            labels = util.readStrings(labelFile)
        elif labels == None:
            nitems = (len(open(treeFile).readlines()) + 1)/ 2
            labels = map(str, range(nitems))
        
        self.makeRoot()

        for i, line in enumerate(open(treeFile)):
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

        # remove unused internal nodes
        labelset = set(labels)
        for child in list(self.root.children):
            if child.is_leaf() and child.name not in labelset:
                self.remove(child)
        
        # remove redunant root
        if len(self.root.children) == 1:
            self.root = self.root.children[0]
            self.remove(self.root.parent)
            self.root.parent = None
    readParentTree = read_parent_tree
    
    
    def write_parent_tree(self, treeFile, labels):
        """Writes tree to the  parent array format"""
        ids = {}
        
        # assign ids to leaves
        for leafname in labels:
            ids[self.nodes[leafname]] = len(ids)
        
        # assign ids to internal nodes
        def walk(node):
            node.recurse(walk)
            if not node.is_leaf():
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
    writeParentTree = write_parent_tree
    
    
    def get_one_line_newick(self):
        """Get a presentation of the tree in a oneline string newick format"""
        stream = StringIO.StringIO()
        self.write(stream, oneline=True)
        return stream.getvalue()
    getOnelineNewick = get_one_line_newick
    
    #
    # should I make these external?
    #
    
    
    def find_depths(self, node = None):
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
    findDepths = find_depths



#============================================================================
# Convenient Input/Output functions
#

def read_tree(filename):
    """Read a tree from a newick file"""
    tree = Tree()
    tree.readNewick(filename)
    return tree
readTree = read_tree

def parse_newick(newick):
    """Read a tree from newick notation stored in a string"""
    tree = Tree()
    stream = StringIO.StringIO(newick)
    tree.readNewick(stream)
    return tree
parseNewick = parse_newick

#============================================================================
# Misc. functions for manipulating trees
#


def assert_tree(tree):
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
assertTree = assert_tree



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


def find_dist(tree, name1, name2):
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
findDist = find_dist
        

def count_descendants(node, sizes=None):
    """Returns a dict with number of leaves beneath each node"""
    if sizes == None:
        sizes = {}
    
    if len(node.children) > 0:
        sizes[node] = 0
        for child in node.children:
            countDescendants(child, sizes)
            sizes[node] += sizes[child]
    else:
        sizes[node] = 1
    
    return sizes
countDescendants = count_descendants


def descendants(node, lst=None):
    """Return a list of all the descendants beneath a node"""
    if lst == None:
        lst = []
    for child in node.children:
        lst.append(child)
        descendants(child, lst=lst)
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


def smallSubtrees(tree, maxsize):
    trees = []
    sizes = countDescendants(tree.root)
    
    def walk(node):
        if sizes[node] <= maxsize:
            trees.append(subtree(tree, node))
        else:
            # if too big keep digging
            for child in node.children:
                walk(child)
    walk(tree.root)
    
    return trees


def max_disjoint_subtrees(tree, subroots):
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
maxDisjointSubtrees = max_disjoint_subtrees


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


def remove_single_children(tree):
    """
    Remove all nodes from the tree that have exactly one child
    
    Branch lengths are added together when node is removed.
    
    """
    
    # find single children
    removed = [node
               for node in tree
               if len(node.children) == 1 and node.parent]
            

    """
    def walk(node):
        if len(node.children) == 1 and node.parent:
            removed.append(node)
        node.recurse(walk)
    walk(tree.root)
    """
    
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
        tree.root.parent = None
        tree.root.dist += oldroot.dist
    
    return removed
removeSingleChildren = remove_single_children


def remove_exposed_internal_nodes(tree, leaves=None):
    """
    Remove all leaves that were originally internal nodes
    
    leaves -- a list of original leaves that should stay
    
    if leaves is not specified, only leaves with strings as names will be kept
    """
    
    if leaves != None:
        stay = set(leaves)
    else:
        # use the fact that the leaf name is a string to determine
        # wether to keep it
        stay = set()
        for leaf in tree.leaves():
            if isinstance(leaf.name, str):
                stay.add(leaf)
    
    # post order traverse tree
    def walk(node):
        # keep a list of children to visit, since they may remove themselves
        children = copy.copy(node.children)
        for child in children:
            walk(child)

        if node.is_leaf() and node not in stay:
            tree.remove(node)
    walk(tree.root)
removeExposedInternalNodes = remove_exposed_internal_nodes


def subtree_by_leaf_names(tree, leafNames, newCopy=False):
    """Returns a subtree with only the leaves specified"""
    
    if newCopy:
        tree = tree.copy()
    
    remove_set = set(tree.leafNames()) - set(leafNames)
    for sp in remove_set:
    	tree.remove(tree.nodes[sp])
    removeExposedInternalNodes(tree)
    removeSingleChildren(tree)
    
    return tree
subtreeByLeafNames = subtree_by_leaf_names


def reorder_tree(tree, tree2):
    """Reorders the branches of tree to match tree2"""

    leaf_lookup = util.list2lookup(tree2.leafNames())

    def mean(lst):
        return sum(lst) / float(len(lst))
    
    def walk(node):
        if node.is_leaf():
            return set([node.name])
        else:
            leaf_sets = []

            for child in node.children:
                leaf_sets.append(walk(child))

            scores = [mean(util.mget(leaf_lookup, l)) for l in leaf_sets]
            rank = util.sortrank(scores)
            node.children = util.mget(node.children, rank)

            # return union
            ret = leaf_sets[0]
            for l in leaf_sets[1:]:
                ret = ret.union(l)
            return ret

    walk(tree.root)


#=============================================================================
# parent tables

def tree2parent_table(tree, leafNames=None):
    """Converts tree to a parent table
    
    parent table is a standard format of the Compbio Lab as of 02/01/2007.
    It is a list of triples (node_name, parent_name, dist)
    
    * If the node is a leaf node_name is the leaf name (a string)
    * If the node is internal node_name is an int representing which row
      (0-based) the node is in the table.
    
    * parent_name indicates the parent of the node.  If the parent is root, a 
      -1 is used as the parent_name.
    
    * dist is the distance between the node and its parent.
    
    Arguments:
    leafNames -- specifies that a tree with only a subset of the leaves 
                 should be used
    
    NOTE: root is not given a row, because root does not have a distance
    the nodeid of the root is -1
    """
    
    if leafNames != None:
        tree = subtree_by_leaf_names(tree, leafNames, newCopy=True)
    else:
        leafNames = tree.leafNames()

    # assign a numbering to the leaves as specified
    nodeid = 0
    nodeids = {}
    nodes = []
    for leaf in leafNames:
        nodeids[tree.nodes[leaf]] = leaf
        nodes.append(tree.nodes[leaf])
        nodeid += 1
    
    # assign a numbering to the internal nodes
    for node in tree:
        if node.is_leaf():
            continue
        if node == tree.root:
            nodeids[node] = -1
        else:
            nodeids[node] = nodeid
            nodeid += 1
            nodes.append(node)

    # make parent table
    parentTable = []
    for node in nodes:
        parentTable.append([nodeids[node], nodeids[node.parent], node.dist])
    
    return parentTable
tree2parentTable = tree2parent_table
    

def parent_table2tree(parentTable):
    """Converts a parent table to a Tree
    
    See tree2parentTable for details
    """

    # TODO: allow named internal nodes
    
    tree = Tree()
    
    # create nodes
    maxint = 0
    for name, parent_name, dist in parentTable:
        node = TreeNode(name)
        node.dist = dist
        tree.add(node)
        
        if isinstance(name, int):
            maxint = max(name, maxint)
    
    # make a root node
    tree.nextname = maxint + 1
    tree.makeRoot()

    # link up parents
    for name, parent_name, dist in parentTable:
        if parent_name == -1:
            parent = tree.root
        else:
            parent = tree.nodes[parent_name]
        tree.addChild(parent, tree.nodes[name])
    
    return tree
parentTable2tree = parent_table2tree


def write_parent_table(parentTable, out=sys.stdout):
    """Writes a parent table to out
    
    out can be a filename or file stream
    """
    
    out = util.open_stream(out, "w")
    
    for name, parent, dist in parentTable:
        print >>out, "%s\t%d\t%f" % (str(name), parent, dist)
writeParentTable = write_parent_table


def read_parent_table(filename):
    """Reads a parent table from the file 'filename'
    
    filename can also be an open file stream
    """
    
    infile = util.open_stream(filename)    
    parentTable = []
    
    for line in infile:
        name, parent, dist = line.split("\t")
        
        if name.is_digit():
            name = int(name)
        parentTable.append([name, int(parent), float(dist)])
    
    return parentTable
readParentTable = read_parent_table
    

#=============================================================================
# Rerooting functions
#


def is_rooted(tree):
    return len(tree.root.children) <= 2
    #return len(tree.root.children) == 3 or len(tree.leaves()) <= 2
isRooted = is_rooted

def unroot(tree, newCopy = True):
    """Return an unrooted copy of tree"""
    
    if newCopy:
        tree = tree.copy()
    
    if len(tree.root.children) == 2:
        nodes = tree.root.children
        dist = nodes[0].dist + nodes[1].dist
        data = tree.merge_branch_data(nodes[0].data, nodes[1].data)
        if len(nodes[0].children) < 2:
            nodes.reverse()
        tree.addChild(nodes[0], nodes[1])
        nodes[1].dist = dist
        tree.set_branch_data(nodes[1], data)
        nodes[0].dist = 0
        tree.set_branch_data(nodes[0], {})
        nodes[0].parent = None
        
        # replace root
        del tree.nodes[tree.root.name]
        tree.root = nodes[0]
    return tree


def reroot(tree, newroot, onBranch=True, newCopy=True):
    """
    Change the rooting of a tree
    """
    
    # TODO: remove newCopy (or assert newCopy=False)
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
        rootdata1, rootdata2 = tree.split_branch_data(node1)
        node1.dist = rootdist / 2.0
        tree.set_branch_data(node1, rootdata1)
        newNode.dist = rootdist / 2.0
        tree.set_branch_data(newNode, rootdata2)
        
        node2 = node1.parent
        node2.children.remove(node1)
        tree.addChild(newNode, node1)
        tree.addChild(node2, newNode)
        
        ptr = node2
        ptr2 = newNode
        newRoot = newNode
    else:
        # root directly on node
        ptr2 = tree.nodes[newroot]
        ptr = ptr2.parent
        newRoot = ptr2
    
    newRoot.parent = None
    
    # reverse parent child relationship of all nodes on path node1 to root
    oldroot = tree.root    
    nextDist = ptr2.dist
    nextData = tree.get_branch_data(ptr2)
    ptr2.dist = 0
    while True:
        nextPtr = ptr.parent
        ptr.children.remove(ptr2)
        tree.addChild(ptr2, ptr)
        
        tmp = ptr.dist
        tmpData = tree.get_branch_data(ptr)
        ptr.dist = nextDist
        tree.set_branch_data(ptr, nextData)
        nextDist = tmp
        nextData = tmpData
        
        ptr2 = ptr
        ptr = nextPtr
        
        if nextPtr is None:
            break
    tree.root = newRoot
    
    return tree


def midpoint_root(tree):

    # get maximum distance from leaves to each node
    depths = {}
    for node in tree.postorder():
        if node.is_leaf():
            depths[node] = (0.0, node)
        else:
            depths[node] = max((c.dist + depths[c][0], depths[c][1])
                               for c in node.children)

    # find maximum path
    dists = []
    for node in tree:
        if node.is_leaf():
            continue
        assert len(node.children) != 1
        tmp = sorted([(c.dist + depths[c][0], depths[c][1], c)
                      for c in node.children])
        dists.append((tmp[-1][0] + tmp[-2][0], node,
                      tmp[-1][2], tmp[-1][1],
                      tmp[-2][2], tmp[-2][1]))
    
    maxdist, top, child1, leaf1, child2, leaf2 = max(dists)
    middist = maxdist / 2.0


    # find longer part of path
    if depths[child1][0] + child1.dist >= middist:
        ptr = leaf1
    else:
        ptr = leaf2

    # find branch that contains midpoint
    dist = 0.0
    while ptr != top:        
        if ptr.dist + dist >= middist:
            # reroot tree
            reroot(tree, ptr.name, onBranch=True, newCopy=False)

            # fixup branch lengths and return
            pdist = sum(c.dist for c in tree.root.children)
            other = filter(lambda x: x != ptr, tree.root.children)[0]
            ptr.dist = middist - dist
            other.dist = pdist - ptr.dist
            return tree

        dist += ptr.dist
        ptr = ptr.parent
    
    
    assert 0 # shouldn't get here


def make_ptree(tree):
    """Make parent tree array from tree"""
    
    nodes = []
    nodelookup = {}
    ptree = []
    
    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)
    walk(tree.root)
    
    def leafsort(a, b):
        if a.is_leaf():
            if b.is_leaf():
                return 0
            else:
                return -1
        else:
            if b.is_leaf():
                return 1
            else:
                return 0
    
    # bring leaves to front
    nodes.sort(cmp=leafsort)
    nodelookup = util.list2lookup(nodes)
    
    for node in nodes:
        if node == tree.root:
            ptree.append(-1)
        else:
            ptree.append(nodelookup[node.parent])
    
    assert nodes[-1] == tree.root
    
    return ptree, nodes, nodelookup
makePtree = make_ptree

#=============================================================================
# Tree visualization
   
def layoutTree(tree, xscale, yscale, minlen=-util.INF, maxlen=util.INF,
               rootx=0, rooty=0):
    """\
    Determines the x and y coordinates for every branch in the tree.
    
    Branch lengths are determined by node.dist
    """
    
    coords = {}
    
    """
       /-----   ] 
       |        ] nodept[node]
    ---+ node   ]
       |
       |
       \---------
    """
    
    # first determine sizes and nodepts
    sizes = {}          # number of descendants (leaves have size 1)
    nodept = {}         # distance between node y-coord and top bracket y-coord
    def walk(node):
        # calculate new y-coordinate for node
                
        # compute node sizes
        sizes[node] = 0
        for child in node.children:
            sizes[node] += walk(child)
        
        if node.is_leaf():
            sizes[node] = 1
            nodept[node] = yscale - 1
        else:
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2.0
        
        return sizes[node]
    walk(tree.root)
    
    # determine x, y coordinates
    def walk(node, x, y):
        xchildren = x+min(max(node.dist*xscale, minlen), maxlen)        
        coords[node] = [xchildren, y + nodept[node]]
                
        if not node.is_leaf():
            ychild = y
            for child in node.children:
                walk(child, xchildren, ychild)
                ychild += sizes[child] * yscale
    walk(tree.root, rootx, rooty)
    
    return coords


def layoutTreeHierarchical(tree, xscale, yscale, rootx=0, rooty=0):
    """\
    Determines the x and y coordinates for every branch in the tree.
    
    Leaves are drawn to line up.  Best used for hierarchical clustering.
    """
    
    coords = {}
    
    """
       /-----   ] 
       |        ] nodept[node]
    ---+ node   ]
       |
       |
       \---------
    """
    
    # first determine sizes and nodepts
    sizes = {}          # number of descendants (leaves have size 1)
    depth = {}          # how deep in tree is node
    nodept = {}         # distance between node y-coord and top bracket y-coord
    def walk(node):
        # calculate new y-coordinate for node
        
        # compute node sizes
        sizes[node] = 0        
        for child in node.children:
            sizes[node] += walk(child)
        
        if node.is_leaf():
            sizes[node] = 1
            nodept[node] = yscale - 1
            depth[node] = 0
        else:
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2.0
            depth[node] = max(depth[child] + 1 for child in node.children)
        
        return sizes[node]
    walk(tree.root)
    
    # determine x, y coordinates
    maxdepth = depth[tree.root]
    def walk(node, x, y):
        xchildren = xscale * (maxdepth - depth[node])
        coords[node] = [xchildren, y + nodept[node]]
        
        if not node.is_leaf():
            ychild = y
            for child in node.children:
                walk(child, xchildren, ychild)
                ychild += sizes[child] * yscale
    walk(tree.root, rootx, rooty)
    
    return coords

#=============================================================================
# Tree color map

def tree_color_map(leafmap=lambda x: (0, 0, 0)):
    """Returns a simple color mixing colormap"""

    def func(tree):
        def walk(node):
            if node.is_leaf():
                node.color = leafmap(node)
            else:
                colors = []
                for child in node.children:
                    walk(child)
                    colors.append(child.color)
                node.color = color_mix(colors)
        walk(tree.root)
    return func
    


def ensemblTreeColorMap(tree):
    """Example of a color map for a tree of Enseml genes"""
    
    def leafmap(node):
        if type(node.name) == str:
            if node.name.startswith("ENSG"): return [1,0,0]
            elif node.name.startswith("ENSCAFG"): return [1,1,0]
            elif node.name.startswith("ENSMUSG"): return [0,0,1]
            elif node.name.startswith("ENSRNOG"): return [0,1,0]
        return [0,0,0]
    
    return treeColorMap(leafmap)(tree)

    
def color_mix(colors):
    """Mixes together several color vectors into one"""
    
    sumcolor = [0, 0, 0]
    for c in colors:
        sumcolor[0] += c[0]
        sumcolor[1] += c[1]
        sumcolor[2] += c[2]                
    for i in range(3):
        sumcolor[i] /= float(len(colors))
    return sumcolor


def make_expr_mapping(maps):
    """Returns a function that maps strings matching an expression to a value
    
       maps -- a list of pairs (expr, value)
    """

    # find exact matches and expressions
    exacts = {}
    exps = []
    for key, val in maps:
        if "*" not in key:
            exacts[key] = val
        else:
            exps.append((key, val))
    
    # create mapping function
    def mapping(key):
        if key in exacts:
            return exacts[key]    
        
        # return default color
        if not isinstance(key, str):
            return (0, 0, 0)
        
        # eval expressions first in order of appearance
        for exp, val in exps:
            if exp[-1] == "*":
                if key.startswith(exp[:-1]):
                    return val
            elif exp[0] == "*":
                if key.endswith(exp[1:]):
                    return val
        
        raise Exception("Cannot map key '%s' to any value" % key)
    return mapping


def read_tree_color_map(filename):
    """Reads a tree colormap from a file"""
    
    infile = util.open_stream(filename)
    maps = []
    
    for line in infile:
        expr, red, green, blue = line.rstrip().split("\t")
        maps.append([expr, map(float, (red, green, blue))])
    
    name2color = make_expr_mapping(maps)
    
    def leafmap(node):
        return name2color(node.name)

    return tree_color_map(leafmap)

#=========================================================================
# Draw Tree ASCII art 
#

try:
    from rasmus import textdraw
except:
    pass


def draw_tree(tree, labels={}, scale=40, spacing=2, out=sys.stdout,
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
        if node.is_leaf():
            sizes[node] = 1
            nodept[node] = yscale - 1 
        else:
            sizes[node] = 0
        for child in node.children:
            sizes[node] += walk(child)
        if not node.is_leaf():
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
        
        if node.is_leaf():
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
drawTree = draw_tree



def draw_tree_lens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        labels[node.name] = "%f" % node.dist
    
    draw_tree(tree, labels, *args, **kargs)
drawTreeLens = draw_tree_lens


def draw_tree_boot_lens(tree, *args, **kargs):
    if not tree.hasData("boot"):
        draw_tree_lens(tree, *args, **kargs)
        return

    labels = {}
    for node in tree.nodes.values():
        if node.is_leaf():
            labels[node.name] = "%f" % node.dist
        else:
            if isinstance(node.data["boot"], int):
                labels[node.name] = "(%d) %f" % (node.data["boot"], node.dist)
            else:
                labels[node.name] = "(%.2f) %f" % (node.data["boot"], node.dist)
    
    draw_tree(tree, labels, *args, **kargs)
drawTreeBootLens = draw_tree_boot_lens


def draw_tree_names(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s" % node.name
    
    draw_tree(tree, labels, *args, **kargs)
drawTreeNames = draw_tree_names


def draw_tree_name_lens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s " % node.name
        else:
            labels[node.name] =""
        labels[node.name] += "%f" % node.dist
    
    draw_tree(tree, labels, *args, **kargs)
drawTreeNameLens = draw_tree_name_lens



if __name__ == "__main__":
    from StringIO import StringIO
    
    infile = StringIO("((a:1,b:2)x:3,(c:4,d:5)y:6)r;")
    tree = read_tree(infile)
    tree.write(rootData=True)

    infile = StringIO("((d:1,c:2)x:3,(b:4,a:5)y:6)r;")
    tree2 = read_tree(infile)

    reorder_tree(tree, tree2)
    draw_tree(tree, maxlen=5)
