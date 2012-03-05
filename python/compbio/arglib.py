"""
    arglib.py
    
    Ancestral recombination graph (ARG)

"""


#=============================================================================
# imports

from __future__ import division

# python libs
import sys
import random
from itertools import izip, chain
from collections import defaultdict
import heapq
from math import *

# compbio libs
from . import fasta

# rasmus libs
from rasmus import treelib, util, stats
from rasmus.intervals import iter_intersections


#=============================================================================
# Ancestral Reconstruction Graph

class ArgNode (object):
    """
    A node in an ARG
    """

    def __init__(self, name="n", age=0, event="gene", pos=0):
        self.name = name
        self.parents = []
        self.children = []
        self.event = event
        self.age = age
        self.pos = pos        # recomb position
        self.data = {}

    def __repr__(self):
        return "<node %s>" % self.name

    def get_dist(self, parent_index=0):
        """Get branch length distance from node to parent_index'th parent"""
        if len(self.parents) == 0:
            return 0.0
        return self.parents[parent_index].age - self.age

    def get_dists(self):
        """Get all branch length distances from node to parents"""
        return [p.age - self.age for p in self.parents]

    def copy(self):
        """Returns a copy of this node"""
        node = ArgNode(self.name, age=self.age, event=self.event,
                        pos=self.pos)
        node.data = dict(self.data)
        return node

    def is_leaf(self):
        """Returns True if this node is a leaf"""
        return len(self.children) == 0
    



class ARG (object):
    """
    A ancestral recombination graph (ARG)
    """

    def __init__(self, start=0.0, end=1.0):
        self.root = None
        self.nodes = {}
        self.nextname = 1
        self.start = start
        self.end = end


    def __iter__(self):
        """Iterates over the nodes in the ARG"""
        return self.nodes.itervalues()


    def __len__(self):
        """Returns number of nodes in the ARG"""
        return len(self.nodes)


    def __getitem__(self, name):
        """Returns node by name"""
        return self.nodes[name]


    def __setitem__(self, name, node):
        """Adds a node to the ARG"""
        node.name = name
        self.add(node)


    def __contains__(self, name):
        """
        Returns True if node in ARG has name 'name'
        """
        return name in self.nodes


    #=================================
    # node manipulation methods

    def new_name(self):
        """
        Returns a new name for a node
        """
        name = self.nextname
        self.nextname += 1
        return name


    def new_node(self, name=None, parents=[], children=[],
                 age=0, event="gene", pos=0):
        """
        Returns a new node
        """
        if name is None:
            name = self.new_name()
        node = self.add(ArgNode(name, age=age, event=event, pos=pos))
        node.parents = list(parents)
        node.children = list(children)
        return node


    def new_root(self, age=0, event="gene", pos=0):
        """
        Returns a new root
        """
        self.root = self.new_node(age=age, event=event, pos=pos)
        return self.root
        

    def add(self, node):
        """
        Adds a node to the ARG
        """
        self.nodes[node.name] = node
        return node

    def remove(self, node):
        """
        Removes a node from the ARG
        """
        for child in node.children:
            child.parents.remove(node)
        for parent in node.parents:
            parent.children.remove(node)
        del self.nodes[node.name]


    def rename(self, oldname, newname):
        """
        Renames a node in the ARG
        """
        node = self.nodes[oldname]
        node.name = newname
        del self.nodes[oldname]
        self.nodes[newname] = node


    def leaves(self, node=None):
        """
        Iterates over the leaves of the ARG
        """
        if node is None:
            for node in self:
                if len(node.children) == 0:
                    yield node
        else:
            for node in self.preorder(node):
                if len(node.children) == 0:
                    yield node


    def leaf_names(self, node=None):
        """
        Iterates over the leaf names of the ARG
        """
        if node is None:
            for node in self:
                if len(node.children) == 0:
                    yield node.name
        else:
            for node in self.preorder(node):
                if len(node.children) == 0:
                    yield node.name


    def copy(self):
        """
        Returns a copy of this ARG
        """

        arg = ARG(start=self.start, end=self.end)
        arg.nextname = self.nextname

        # copy all nodes
        for name, node in self.nodes.iteritems():
            arg.nodes[name] = node.copy()

        # connect nodes
        for node in self.nodes.itervalues():
            node2 = arg[node.name]
            for child in node.children:
                node2.children.append(arg[child.name])
            for parent in node.parents:
                node2.parents.append(arg[parent.name])

        if self.root:
            arg.root = arg[self.root.name]

        return arg
            

    #================================
    # iterator methods
        
    def postorder(self, node=None):
        """
        Iterates through nodes in postorder traversal
        """

        visit = defaultdict(lambda: 0)
        queue = list(self.leaves(node))

        for node in queue:
            yield node
            for parent in node.parents:
                visit[parent] += 1

                # if all children of parent has been visited then queue parent
                if visit[parent] == len(parent.children):
                    queue.append(parent)

        
    def preorder(self, node=None):
        """
        Iterates through nodes in preorder traversal
        """

        visit = set()
        if node is None:
            node = self.root
        queue = [node]

        for node in queue:
            if node in visit:
                continue
            yield node
            visit.add(node)
            
            for child in node.children:
                queue.append(child)


    def postorder_marginal_tree(self, pos):
        """
        Iterate postorder over the nodes in the marginal tree at position 'pos'

        NOTE: nodes are iterated in order of age
        """
        
        # initialize heap
        heap = [(node.age, node) for node in self.leaves()]
        seen = set([None])
        visited = set([])
        visit_age = min(x[0] for x in heap) - 1
        
        def unreachable(node):
            # returns True if node is unreachable from leaves
            if node in visited or node.is_leaf():
                return False
            if node.age < visit_age:
                return True
            for child in self.get_local_children(node, pos):
                if not unreachable(child):
                    return False
            return True

        def ready(node):
            # returns True if node is ready to yield
            # node is ready if all unvisited child are unreachable
            for child in self.get_local_children(node, pos):
                if child not in visited:
                    if not unreachable(child):
                        return False
            return True
        

        # add all ancestor of lineages
        unready = []
        while len(heap) > 0:
            # yield next ready node
            del unready[:]
            while True:
                age, node = heapq.heappop(heap)
                if ready(node):
                    break
                unready.append((age, node))
            for x in unready:
                heapq.heappush(heap, x)
            yield node
            visited.add(node)
            visit_age = node.age
            if len(heap) == 0:
                # MRCA reached
                return
            
            # find correct marginal parent
            # add parent to lineages if it has not been seen before
            parent = self.get_local_parent(node, pos)
            if parent not in seen:
                heapq.heappush(heap, (parent.age, parent))
                seen.add(parent)


    def preorder_marginal_tree(self, pos, node=None):
        """
        Iterate preorder over the nodes in the marginal tree at position 'pos'

        NOTE: this might also include unreachable nodes
        """

        if node is None:
            node = arg.root

        # initialize heap
        heap = [node]
        seen = set([node])
        
        # add all ancestor of lineages
        while len(heap) > 0:
            node = heap.pop()
            yield node

            for child in node.children:
                if self.get_local_parent(child, pos) == node:
                    if child not in seen:
                        heap.append(child)
                        seen.add(child)
                        # NOTE: this prevents error when
                        # children[0] == children[1]


    def get_local_parent(self, node, pos):
        """Returns the local parent of 'node' for position 'pos'"""
        if node.event == "gene" or node.event == "coal":
            if len(node.parents) > 0:
                return node.parents[0]
            else:
                return None
        elif node.event == "recomb":
            if len(node.parents) > 0:
                return node.parents[0 if pos < node.pos else 1]
            else:
                return None
        else:
            raise Exception("unknown event '%s'" % node.event)


    def get_local_parents(self, node, start, end):
        """Returns the parents of 'node' with ancestral sequence within (start, end)"""
        if node.event == "recomb":
            parents = []
            if node.pos > start:
                parents.append(node.parents[0])
            if node.pos < end:
                parents.append(node.parents[1])
        else:
            parents = node.parents
        return parents


    def get_local_children(self, node, pos):
        """
        Returns the local children of 'node' for position 'pos'
        
        NOTE: the local children are not necessarily in the local tree
        because the children may be unreachable from the leaves
        """

        return [child for child in node.children
                if self.get_local_parent(child, pos) == node]
        

    def get_local_dist(self, node, pos):
        """Returns the local parent of 'node' for position 'pos'"""
        parent = self.get_local_parent(node, pos)
        if parent:
            return parent.age - node.age
        else:
            return 0.0


    #===============================
    # ancestral sequence methods

    def set_recomb_pos(self, start=None, end=None, descrete=False):
        """
        Randomly aample all recombination positions in the ARG
        """

        if start is not None:
            self.start = start
        if end is not None:
            self.end = end

        length = self.end - self.start

        for node in self:
            if node.event == "recomb":
                if descrete:
                    node.pos = random.randint(self.start, self.end-1) + .5
                else:
                    node.pos = random.random() * length + self.start


    def set_ancestral(self):
        """
        Set all ancestral regions for the nodes of the ARG

        NOTE: recombination positions must be set first (set_recomb_pos)
        """

        for node in self:
            node.data["ancestral"] = []

        for block, tree in iter_tree_tracks(self):
            for node in tree:
                if node.name in self.nodes:
                    ancestral = self[node.name].data["ancestral"]
                    if len(ancestral) > 0 and ancestral[-1][1] == block[0]:
                        # extend
                        ancestral[-1] = (ancestral[-1][0], block[1])
                    else:
                        ancestral.append(block)
                else:
                    # cap node?
                    pass


    def set_ancestral3(self):
        """
        Set all ancestral regions for the nodes of the ARG

        NOTE: recombination positions must be set first (set_recomb_pos)
        """

        # NOTE: block_counts is used to determine when the MRCA of a block
        # is found.
        
        # get all non-recomb blocks (identified by starting pos)
        nleaves = len(list(self.leaves()))
        block_counts = [[self.start, self.end, nleaves]]

        
        
        for node in self.postorder():
            print "node", node.age

            if node.is_leaf():
                # initialize leaves with entire extant sequence
                node.data["ancestral"] = [(self.start, self.end)]

            elif node.event == "coal":
                # union of ancestral of children
                # get all child regions
                assert len(node.children) == 2, node
                
                # walk through regions for both children and determine
                # whether they coal
                regions = list(block_counts)
                regions3 = []
                if node.children[0] == node.children[1]:
                    # special case
                    regions.extend((start, end, "coal") for (start, end) in 
                                    node.children[0].data["ancestral"])
                else:
                    regions.extend(
                        (start, end, "coal") for (start, end) 
                        in chain(self.get_ancestral(
                                node.children[0], parent=node),
                                 self.get_ancestral(
                                node.children[1], parent=node)))

                regions.sort()

                print "block_counts", block_counts
                print "regions", regions

                block_counts = []
                for start, end, regs in iter_intersections(regions):

                    ncoals = util.ilen(reg for reg in regs
                                       if reg[2] == "coal")
                    count = sum(reg[2] for reg in regs if reg[2] != "coal")
                    
                    print ncoals, count, regs

                    if ncoals == 2:
                        # region coal
                        block_counts.append((start, end, count - 1))
                        regions3.append((start, end))
                    
                    elif ncoals == 1:
                        # pass single block upwards
                        block_counts.append((start, end, count))
                        if count > 1:
                            regions3.append((start, end))

                    elif ncoals == 0:
                        block_counts.append((start, end, count))

                    else:
                        print (self.get_ancestral(
                                node.children[0], parent=node), 
                               self.get_ancestral(
                                node.children[1], parent=node))
                        print (node.children[0].data["ancestral"],
                               node.children[1].data["ancestral"])

                        raise Exception("unknown error")
                        
                
                node.data["ancestral"] = regions3
                if max([0] + util.hist_dict(regions3).values()) > 1:
                    print regions3
                    raise Exception("unknown error2")


            elif node.event == "recomb":
                # inherit all ancestral
                assert len(node.children) == 1, (node, len(node.children))
                regions3 = []
                regions = [(start, end, "node") for start, end in 
                           self.get_ancestral(node.children[0], parent=node)]
                regions.extend(block_counts)
                regions.sort()

                for start, end, regs in iter_intersections(regions):
                    nnodes = util.ilen(reg for reg in regions
                                       if reg[2] == "node")
                    count = sum(reg[2] for reg in regions
                                       if reg[2] != "node")
                    if nnodes > 0 and count > 1:
                        regions3.append((start, end))
                node.data["ancestral"] = regions3

                if max([0] + util.hist_dict(regions3).values()) > 1:
                    print regions3
                    raise Exception("unknown error2")

            else:
                raise Exception("unknown event '%s'" % node.event)


        print block_counts
    
    def set_ancestral2(self):
        """
        Set all ancestral regions for the nodes of the ARG

        NOTE: recombination positions must be set first (set_recomb_pos)
        """

        # NOTE: block_counts is used to determine when the MRCA of a block
        # is found.
        
        # get all non-recomb blocks (identified by starting pos)
        nleaves = len(list(self.leaves()))
        all_blocks = list(iter_recomb_blocks(self))
        block_counts = dict((block, nleaves) for block in all_blocks)
        
        for node in self.postorder():
            if node.is_leaf():
                # initialize leaves with entire extant sequence
                node.data["ancestral"] = list(all_blocks)
            elif node.event == "coal":
                # union of ancestral of children
                # get all child regions
                assert len(node.children) == 2, node
                
                # walk through regions for both children and determine
                # whether they coal
                if node.children[0] == node.children[1]:
                    # special case
                    regions1 = node.children[0].data["ancestral"]
                    regions2 = []
                else:
                    regions1 = self.get_ancestral(
                        node.children[0], parent=node)
                    regions2 = self.get_ancestral(
                        node.children[1], parent=node)
                regions3 = []

                i = j = 0
                while True:
                    reg1 = regions1[i] if i < len(regions1) else None
                    reg2 = regions2[j] if j < len(regions2) else None
                    if reg1 is None and reg2 is None:
                        # stop when all regions have been considered
                        break

                    if reg1 == reg2:
                        # region coal
                        block_counts[reg1] -= 1
                        regions3.append(reg1)
                        i += 1
                        j += 1
                    elif reg2 is None or (reg1 and reg1[0] < reg2[0]):
                        if block_counts[reg1] > 1:
                            regions3.append(reg1)
                        i += 1
                    else:
                        assert reg2, reg2
                        if block_counts[reg2] > 1:
                            regions3.append(reg2)
                        j += 1

                node.data["ancestral"] = regions3

            elif node.event == "recomb":
                # inherit all ancestral
                assert len(node.children) == 1, (node, len(node.children))
                node.data["ancestral"] = [
                    reg for reg in self.get_ancestral(
                    node.children[0], parent=node)
                    if block_counts[reg] > 1]

            else:
                raise Exception("unknown event '%s'" % node.event)

                    
    def get_ancestral(self, node, side=None, parent=None):
        """
        Get the ancestral sequence from an edge above a node 
        
        node -- node to get ancestral sequence from
        side -- 0 for left parent edge, 1 for right parental edge
        parent -- if given, determine side from parent node
        """

        # set side from parent
        if parent:
            side = node.parents.index(parent)

        if node.event == "recomb":
            if (parent and len(node.parents) == 2 and
                node.parents[0] == node.parents[1]):
                # special case where both children of a coal node are the same
                # recomb node.
                return node.data["ancestral"]
            
            regions = []
            for reg in node.data["ancestral"]:
                if side == 0:
                    if reg[1] <= node.pos:
                        # keep all regions fully left of recomb position
                        regions.append(reg)
                    elif reg[0] < node.pos:
                        # cut region
                        regions.append((reg[0], node.pos))
                elif side == 1:
                    if reg[0] >= node.pos:
                        # keep all regions fully right of recomb position
                        regions.append(reg)
                    elif reg[1] > node.pos:
                        # cut region
                        regions.append((node.pos, reg[1]))
                else:
                    raise Exception("side not specified")
            return regions
        
        elif node.event == "gene" or node.event == "coal":
            return node.data["ancestral"]

        else:
            raise Exception("unknown event '%s'" % node.event)                
            

    def prune(self, remove_single=True):
        """
        Prune ARG to only those nodes with ancestral sequence
        """

        # NOTE: be careful when removing nodes that you call get_ancestral
        # before changing parent/child orders
        
        # find pruned edges
        prune_edges = []
        for node in list(self):
            for parent in list(node.parents):
                if len(self.get_ancestral(node, parent=parent)) == 0:
                    prune_edges.append((node, parent))

        # remove pruneded edges
        for node, parent in prune_edges:
            parent.children.remove(node)
            node.parents.remove(parent)
        
        # remove pruned nodes
        for node in list(self):
            if len(node.data["ancestral"]) == 0:
                self.remove(node)


        for node in self:
            assert not node.is_leaf() or node.age == 0.0

        # remove single children
        if remove_single:
            remove_single_lineages(self)
            
        # set root
        # TODO: may need to actually use self.roots
        for node in list(self):
            if len(node.parents) == 0:
                dellist = []
                while len(node.children) == 1:
                    delnode = node
                    node = node.children[0]
                    self.remove(delnode)
                self.root = node


    #===========================
    # marginal tree methods
                
    def get_marginal_tree(self, pos):
        """
        Returns the marginal tree of the ARG containing position 'pos'
        """
        
        # make new ARG to contain marginal tree
        tree = ARG(self.start, self.end)
        tree.nextname = self.nextname

        # populate tree with marginal nodes
        for node in self.postorder_marginal_tree(pos):
            tree.add(node.copy())
        
        # set parent and children
        roots = []
        for node2 in tree:
            node = self[node2.name]
            parent = self.get_local_parent(node, pos)
            if parent is not None and parent.name in tree.nodes:
                parent2 = tree[parent.name]
                node2.parents = [parent2]
                parent2.children.append(node2)
            else:
                roots.append(node2)

        # make root
        if len(roots) == 1:
            tree.root = roots[0]
        elif len(roots) > 1:
            # make cap node since marginal tree does not fully coallesce
            tree.root = tree.new_node(event="coal", 
                                      name=self.new_name(),
                                      age=max(x.age for x in roots)+1)
            tree.nextname = self.nextname
            for node in roots:
                tree.root.children.append(node)
                node.parents.append(tree.root)

        assert tree.root is not None, (tree.nodes, pos)
        
        return tree
        
        
    def get_tree(self, pos=None):
        """
        Returns a treelib.Tree() object representing the ARG if it is a tree

        if 'pos' is given, return a treelib.Tree() for the marginal tree at
        position 'pos'.
        """

        # TODO: make more efficient

        # get marginal tree first
        if pos is not None:
            return self.get_marginal_tree(pos).get_tree()

        tree = treelib.Tree()

        # add all nodes
        for node in self:
            node2 = treelib.TreeNode(node.name)
            tree.add(node2)

        # set parent, children, dist
        for node in tree:
            node2 = self[node.name]
            node.parent = (tree[node2.parents[0].name]
                           if len(node2.parents) > 0 else None)
            node.children = [tree[c.name] for c in node2.children]

            if node.parent:
                node.dist = self[node.parent.name].age - node2.age

        tree.root = tree[self.root.name]
        return tree


    #=======================
    # input/output

    def read(self, filename=sys.stdin):
        read_arg(filename, arg=self)


    def write(self, filename=sys.stdout):
        write_arg(filename, self)

            

#=============================================================================

def assert_arg(arg):
    """Asserts that the arg data structure is consistent"""

    for node in arg:
        # check parent, child links
        for parent in node.parents:
            assert node in parent.children
        for child in node.children:
            assert node in child.parents

        # check ages
        for parent in node.parents:
            assert node.age <= parent.age, ((node.name, node.age),
                                            (parent.name, parent.age))

    leaves = set(arg.leaf_names())
    for tree in iter_marginal_trees(arg):
        assert set(tree.leaf_names()) == leaves


#=============================================================================
# coalescence with recombination

def sample_coal_recomb(k, n, r):
    """
    Returns a sample time for either coal or recombination

    k -- chromosomes
    n -- effective population size (haploid)
    r -- recombination rate (recombinations / chromosome / generation)

    Returns (event, time) where
    event -- 0 for coalesce event, 1 for recombination event
    time  -- time (in generations) of event
    """

    # coal rate = (k choose 2) / 2
    # recomb rate = k * r
    coal_rate = (k * (k-1) / 2) / n
    recomb_rate = k * r
    rate = coal_rate + recomb_rate

    event = ("coal", "recomb")[int(random.random() < (recomb_rate / rate))]
    
    return event, random.expovariate(rate)


def sample_coal_recomb_times(k, n, r, t=0):
    """
    Returns a sample time for either coal or recombination

    k -- chromosomes
    n -- effective population size (haploid)
    r -- recombination rate (recombinations / chromosome / generation)
    t -- initial time (default: 0)

    Returns (event, time) where
    event -- 0 for coalesce event, 1 for recombination event
    time  -- time (in generations) of event
    """

    times = []
    events = []

    while k > 1:
        event, t2 = sample_coal_recomb(k, n, r)
        t += t2
        times.append(t)
        events.append(event)
        if event == "coal":
            k -= 1
        elif event == "recomb":
            k += 1
        else:
            raise Exception("unknown event '%s'" % event)

    return times, events


def sample_arg(k, n, rho, start=0.0, end=1.0, t=0, names=None, 
               make_names=True):
    """
    Returns an ARG sampled from the coalescent with recombination (pruned)

    k   -- chromosomes
    n   -- effective population size (haploid)
    rho -- recombination rate (recombinations / site / generation)
    start -- staring chromosome coordinate
    end   -- ending chromsome coordinate
    t   -- initial time (default: 0)
    names -- names to use for leaves (default: None)
    make_names -- make names using strings (default: True)

    Returns (event, time) where
    event -- 0 for coalesce event, 1 for recombination event
    time  -- time (in generations) of event
    """

    arg = ARG(start, end)

    class Lineage (object):
        def __init__(self, node, regions, seqlen):
            self.node = node
            self.regions = regions
            self.seqlen = seqlen

    # init ancestral lineages
    # (node, region, seqlen)
    total_seqlen = k * (end - start)
    if make_names:
        names = ["n%d" % i for i in range(k)]
    if names is None:
        lineages = set(Lineage(arg.new_node(), [(start, end)], end-start)
                       for i in xrange(k))
    else:
        lineages = set(Lineage(arg.new_node(name=names[i]),
                               [(start, end)], end-start)
                       for i in xrange(k))
    for lineage in lineages:
        lineage.node.data["ancestral"] = [(start, end)]
    recomb_parent_lineages = {}
    lineage_parents = {}

    # block start -> lineage count
    block_starts = [start]
    block_counts = {start: k}

    # perform coal, recomb
    while len(lineages) > 1:
        # sample time and event
        k = len(lineages)
        coal_rate = (k * (k-1) / 2) / n  # (k choose 2) / n
        recomb_rate = rho * total_seqlen
        rate = coal_rate + recomb_rate
        t2 = random.expovariate(rate)
        event = ("coal", "recomb")[int(random.random() < (recomb_rate / rate))]
        t += t2

        
        # process event
        if event == "coal":
            node = arg.new_node(age=t, event=event)

            # choose lineages to coal
            a, b = random.sample(lineages, 2)
            lineages.remove(a)
            lineages.remove(b)
            lineage_parents[a] = node
            lineage_parents[b] = node
            total_seqlen -= a.seqlen + b.seqlen

            # set parent, child links
            node.children = [a.node, b.node]
            a.node.parents.append(node)
            b.node.parents.append(node)

            # coal each non-overlapping region
            regions = []
            lineage_regions = []
            nblocks = len(block_starts)
            i = 0
            
            for start, end, count in count_region_overlaps(
                a.regions, b.regions):                
                assert start != end, count in (0, 1, 2)
                #assert end == arg.end or end in block_starts
                i = block_starts.index(start, i)
                start2 = start
                while start2 < end:
                    end2 = block_starts[i+1] if i+1 < nblocks else arg.end

                    # region coalesces
                    if count == 2:
                        block_counts[start2] -= 1
                    if count >= 1:
                        regions.append((start2, end2)) # ancestral seq
                        if block_counts[start2] > 1:
                            # regions moves on, since not MRCA
                            lineage_regions.append((start2, end2))

                    # move to next region
                    i += 1
                    start2 = end2
            node.data["ancestral"] = regions

            # create 1 new lineage if any regions remain
            if len(lineage_regions) > 0:
                for reg in lineage_regions:
                    assert block_counts[reg[0]] > 1, (reg, block_counts)
                seqlen = lineage_regions[-1][1] - lineage_regions[0][0]
                lineages.add(Lineage(node, lineage_regions, seqlen))
                total_seqlen += seqlen

            
        elif event == "recomb":
            node = arg.new_node(age=t, event=event)

            # choose lineage and pos to recombine (weighted by seqlen)
            pick = random.random() * total_seqlen
            i = 0
            for lineage in lineages:
                i += lineage.seqlen
                if i >= pick and lineage.seqlen > 0:
                    break

            # set parent, child links
            lineage_parents[lineage] = node
            lineages.remove(lineage)
            node.children = [lineage.node]
            lineage.node.parents.append(node)
            node.data["ancestral"] = lineage.regions

            # choose recomb pos
            node.pos = random.uniform(lineage.regions[0][0],
                                      lineage.regions[-1][1])

            # does recomb pos break an existing block?
            for reg in lineage.regions:
                if reg[0] < node.pos < reg[1]:
                    # split block
                    block_starts.append(node.pos)
                    block_starts.sort()
                    prev_pos = block_starts[block_starts.index(node.pos)-1]
                    block_counts[node.pos] = block_counts[prev_pos]

            # create 2 new lineages
            regions1 = list(split_regions(node.pos, 0, lineage.regions))
            regions2 = list(split_regions(node.pos, 1, lineage.regions))
            
            regions1_len = regions1[-1][1] - regions1[0][0]
            regions2_len = regions2[-1][1] - regions2[0][0]
            total_seqlen += regions1_len + regions2_len - lineage.seqlen
            a = Lineage(node, regions1, regions1_len)
            b = Lineage(node, regions2, regions2_len)
            lineages.add(a)
            lineages.add(b)
            recomb_parent_lineages[node] = (a, b)
        else:
            raise Exception("unknown event '%s'" % event)

    assert len(lineages) == 0, (lineages, block_counts.values())

    # fix recomb parent order, so that left is before pos and right after
    for node, (a, b) in recomb_parent_lineages.iteritems():
        an = lineage_parents[a]
        bn = lineage_parents[b]
        for reg in a.regions: assert reg[1] <= node.pos
        for reg in b.regions: assert reg[0] >= node.pos
        node.parents = [an, bn]

    # TODO: set root(s)

    return arg


# TODO: make a sample_trees_smc()
# would be much faster

def sample_arg_smc(k, n, rho, start=0.0, end=0.0, t=0):

    arg = ARG(start, end)
    pos = start # position simulated so far

    # sample initial coal tree with no recomb
    times, events = sample_coal_recomb_times(k, n, 0)
    # at time[i] we go from (k - i) lineages to (k - i - 1) lineages
    
    lineages = set(arg.new_node() for i in xrange(k))
    for t in times:
        a, b = random.sample(lineages, 2)
        parent = arg.new_node(age=t, children=[a,b], event="coal")
        a.parents.append(parent)
        b.parents.append(parent)
        lineages.remove(a)
        lineages.remove(b)
        lineages.add(parent)

    '''
    # get nodes from last tree
    nodes = list(arg.postorder_marginal_tree(pos+1))
    oldest = nodes.pop()
    nodes_set = set(nodes)
    nodes_set.add(oldest)
    parent_count = defaultdict(lambda:0)
    parents = {}
    for node in nodes_set:
        parent = arg.get_local_parent(node, start+1)
        parents[node] = parent
        parent_count[parent] += 1
    '''

    while True:

        # TODO: could make faster by modifying nodes, parents, and parents_count
        # get nodes from last tree
        nodes = list(arg.postorder_marginal_tree(pos+1))
        oldest = nodes.pop()
        
        # sample next recomb
        totlen = sum(arg.get_local_dist(x, pos+1) for x in nodes)
        recomb_per_site = rho * totlen

        # sample recomb pos
        pos += random.expovariate(recomb_per_site)
        if pos > end:
            # no more recombinations
            break

        # sample recomb node
        nodes_set = set(nodes)
        nodes_set.add(oldest)
        parent_count = defaultdict(lambda:0)
        parents = {}
        for node in nodes_set:
            parent = arg.get_local_parent(node, pos+1)
            parents[node] = parent
            parent_count[parent] += 1
        times = [parent.age for parent in parent_count
                 if parent_count[parent] == 2]
        times.sort()
        lens = [arg.get_local_dist(x, pos+1) for x in nodes]
        child = nodes[stats.sample(lens)]

        # make recomb node and insert into arg
        old_parent = arg.get_local_parent(child, pos+1)
        recomb_node = arg.new_node(age=child.age + random.random() *
                                   arg.get_local_dist(child, pos+1),
                                   children=[child],
                                   parents=[old_parent],
                                   pos=pos, event="recomb")
        i = old_parent.children.index(child)
        old_parent.children[i] = recomb_node
        child.parents[child.parents.index(old_parent)] = recomb_node


        # coalesce recomb_node back into previous marginal tree
        for i in xrange(len(times)):
            if recomb_node.age < times[i]:
                break
        for i in xrange(i, len(times) - 1):
            k2 = k - i - 1
            # k2 = 0 = k - i - 1
            # i = k - 1
            coal_time = times[i] + random.expovariate(k2 / float(n))
            if coal_time < times[i+1]:
                # coal event, choose sister
                candidates = [(node, parents[node]) for node in nodes
                              if node.age < coal_time < parents[node].age]
                assert len(candidates) > 0, (coal_time, times)
                sister, sister_parent = random.sample(candidates, 1)[0]
                if sister == child:
                    sister = recomb_node
                    assert sister_parent == old_parent
                break
        else:
            # no coal at all choose oldest sister
            sister, sister_parent = oldest, parents[oldest]
            assert sister_parent is None or sister_parent in sister.parents, \
                   (sister_parent, sister.parents)
            coal_time = times[-1] + random.expovariate(1.0 / n)


        # make coal node
        coal_node = arg.new_node(age=coal_time, event="coal",
                                 children=[sister, recomb_node])
        recomb_node.parents.append(coal_node)
        if sister_parent:
            i = sister_parent.children.index(sister)
            sister_parent.children[i] = coal_node
            i = sister.parents.index(sister_parent)
            sister.parents[i] = coal_node
            coal_node.parents.append(sister_parent)
        else:
            sister.parents.append(coal_node)
            

    nodes = [(node.age, node) for node in arg]
    nodes.sort()
    arg.root = nodes[-1][1]

    return arg

    



#=============================================================================
# arg functions


def lineages_over_time(k, events):
    """
    Computes number of lineage though time using coal/recomb events
    """

    for event in events:
        if event == "coal":
            k -= 1
        elif event == "recomb":
            k += 1
        else:
            raise Exception("unknown event '%s'" % event)        
        yield k
        

def make_arg_from_times(k, times, events):
    """
    Returns an ARG given 'k' samples and a list of 'times' and 'events'

    times  -- ordered times of coalescence or recombination
    events -- list of event types (either 'coal' or 'recomb')
    """

    arg = ARG()

    # make leaves
    lineages  = set((arg.add(ArgNode(arg.new_name())), 1)
                     for i in xrange(k))

    # process events
    for t, event in izip(times, events):
        if event == "coal":
            node = arg.add(ArgNode(arg.new_name(), age=t, event=event))
            a, b = random.sample(lineages, 2)
            lineages.remove(a)
            lineages.remove(b)
            node.children = [a[0], b[0]]
            a[0].parents.append(node)
            b[0].parents.append(node)
            lineages.add((node, 1))
            
        elif event == "recomb":
            node = arg.add(ArgNode(arg.new_name(), age=t, event=event))
            a = random.sample(lineages, 1)[0]
            lineages.remove(a)
            node.children = [a[0]]
            a[0].parents.append(node)
            lineages.add((node, 1))
            lineages.add((node, 2))

        else:
            raise Exception("unknown event '%s'" % event)

    
    if len(lineages) == 1:
        arg.root = lineages.pop()[0]    

    return arg


def get_recomb_pos(arg, visible=False):
    """
    Returns a sorted list of an ARG's recombination positions

    visible -- if True only iterate recombination break points that are 
               visible to extant sequences
    """
    
    if visible:
        rpos = [node.pos for node in
                arg if node.event == "recomb" and 
                len(node.data["ancestral"]) > 0]
        rpos.sort()
        return rpos

    else:
        rpos = [node.pos for node in
                arg if node.event == "recomb"]
        rpos.sort()
        return rpos


def iter_recomb_blocks(arg, start=None, end=None, visible=False):
    """
    Iterates over the recombination blocks of an ARG

    arg     -- ARG to iterate over
    start   -- starting position in chromosome to iterate over
    end     -- ending position in chromosome to iterate over
    visible -- if True only iterate recombination break points that are 
               visible to extant sequences
    """

    if start is None:
        start = arg.start
    if end is None:
        end = arg.end

    a = start
    b = start
    for pos in get_recomb_pos(arg, visible=visible):
        if pos < start:
            continue
        if pos > end:
            pos = end
            break
        b = pos
        yield (a, b)
        a = pos

    yield (a, end)


def iter_marginal_trees(arg, start=None, end=None):
    """
    Iterate over the marginal trees of an ARG
    """
    for block, tree in iter_tree_tracks(arg, start, end):
        yield tree


def iter_tree_tracks(arg, start=None, end=None, convert=False):
    """
    Iterate over the marginal trees of an ARG

    Yeilds ((start, end), tree) for each marginal tree where (start, end)
    defines the block of the marginal tree
    """

    if start is None:
        start = arg.start
    if end is None:
        end = arg.end

    i = 0
    rpos = get_recomb_pos(arg)
    rpos.append(end)
    if len(rpos) > 0:
        r = rpos[i]
    else:
        r = end

    while start < end:
        # find next rpos
        while i < len(rpos)-1 and r <= start:
            i += 1
            r = rpos[i]

        tree = arg.get_marginal_tree((start+r) / 2.0)
        
        # find block end
        end2 = arg.end
        for node in tree:
            if node.event == "recomb" and start < node.pos < end2:
                end2 = node.pos

        if convert:
            tree = tree.get_tree()
        yield (start, end2), tree
        start = end2

    
def descendants(node, nodes=None):
    """
    Return all descendants of a node in an ARG
    """
    if nodes is None:
        nodes = set()
    nodes.add(node)
    for child in node.children:
        if child not in nodes:
            descendants(child, nodes)
    return nodes


def remove_single_lineages(arg):
    """
    Remove unnecessary nodes with single parent and single child
    """
    queue = list(arg)

    for node in queue:
        if node.name not in arg:
            continue
        
        if len(node.children) == 1:
            if len(node.parents) == 1:
                child = node.children[0]
                parent = node.parents[0]

                del arg.nodes[node.name]
                child.parents[child.parents.index(node)] = parent
                parent.children[parent.children.index(node)] = child

            elif len(node.parents) == 0:
                child = node.children[0]

                del arg.nodes[node.name]
                child.parents.remove(node)
                arg.root = node

                queue.append(child)

    # relabel events for leaves that were recombinations
    for node in arg:
        if node.is_leaf() and len(node.parents) == 1:
            node.event = "gene"

    return arg



def postorder_subarg(arg, start, end):
    """Iterates postorder over the nodes of the 'arg' that are ancestral to (start,end)"""

    # initialize heap
    heap = [(node.age, node) for node in arg.leaves()]
    seen = set([None])

    # add all ancestor of lineages
    while len(heap) > 0:
        age, node = heapq.heappop(heap)
        yield node
        if len(heap) == 0:
            # MRCA reached
            return

        # find parents within (start, end)
        # add parent to lineages if it has not been seen before
        for parent in arg.get_local_parents(node, start, end):
            if parent not in seen:
                heapq.heappush(heap, (parent.age, parent))
                seen.add(parent)


def subarg(arg, start, end):
    """Returns a new ARG that only contains recombination within (start, end)"""

    arg2 = ARG(start, end)

    # add nodes
    for node in postorder_subarg(arg, start, end):
        arg2.root = arg2.new_node(node.name, event=node.event, age=node.age,
                                  pos=node.pos)

    # add edges
    for node2 in arg2:
        node = arg[node2.name]
        for parent in arg.get_local_parents(node, start, end):
            pname = parent.name
            if pname in arg2:
                parent2 = arg2[pname]
                node2.parents.append(parent2)
                parent2.children.append(node2)

    return arg2


def subarg_by_leaves(arg, leaves, keep_single=False):
    """
    Removes any leaf from the arg that is not in leaves set
    """

    stay = set(leaves)
    remove = []

    # find nodes to remove
    for node in arg.postorder():
        nchildren = sum(1 for child in node.children if child in stay)
        if nchildren == 0 and node not in stay:
            remove.append(node)
        else:
            stay.add(node)

    # remove nodes
    for node in remove:
        arg.remove(node)
    
    if not keep_single:
        remove_single_lineages(arg)

    return arg


def iter_arg_sprs(arg, start=None, end=None, use_leaves=False):
    """
    Iterate through the SPR moves of an ARG

    Yields (recomb_pos, (rnode, rtime), (cnode, ctime))
    """

    trees = iter_tree_tracks(arg, start, end)
    block, last_tree = trees.next()
    
    for block, tree in trees:

        # find recombination node
        recomb_pos = block[0]
        node = (x for x in tree if x.pos == recomb_pos).next()
        rtime = node.age
        ptr = last_tree[node.name]
        while len(ptr.children) == 1:
            ptr = ptr.children[0]
        rnode = ptr.name

        # find recoal node
        ptr = node
        while len(ptr.children) != 2:
            last = ptr
            ptr = ptr.parents[0]
        ctime = ptr.age
        if ptr.name in last_tree:
            c = ptr.children
            child = c[0] if c[1] == last else c[1]
            child = last_tree[child.name]
            while len(child.children) == 1:
                child = child.children[0]
            cnode = child.name
        else:
            cnode = last_tree.root.name

        # yield SPR
        if use_leaves:
            rleaves = list(last_tree.leaf_names(last_tree[rnode]))
            cleaves = list(tree.leaf_names(last_tree[cnode]))
            yield (recomb_pos, (rleaves, rtime), (cleaves, ctime))
        else:
            yield (recomb_pos, (rnode, rtime), (cnode, ctime))
        last_tree = tree
        


def make_arg_from_sprs(init_tree, sprs, ignore_self=False):
    """
    Make an ARG from an initial tree 'init_tree' and a list of SPRs 'sprs'

    NOTE: sprs should indicate branches by their leaf set (use_leaves=True)
    """

    def add_node(arg, node, time, pos, event):
        node2 = arg.new_node(event=event, age=time, children=[node], pos=pos)
        if event == "coal":
            node2.pos = 0
        parent = arg.get_local_parent(node, pos)
        if parent:
            node.parents[node.parents.index(parent)] = node2
            parent.children[parent.children.index(node)] = node2
            node2.parents.append(parent)
        else:
            node.parents.append(node2)
        return node2

    arg = init_tree

    for rpos, (rleaves, rtime), (cleaves, ctime) in sprs:

        node1 = arg_lca(arg, rleaves, rpos, time=rtime)
        node2 = arg_lca(arg, cleaves, rpos, time=ctime)

        # check whether self cycles are wanted
        if ignore_self and node1 == node2:
            continue
        
        recomb = add_node(arg, node1, rtime, rpos, "recomb")
        if node1 == node2:
            node2 = recomb
        coal = add_node(arg, node2, ctime, rpos, "coal")

        recomb.parents.append(coal)
        coal.children.append(recomb)

    return arg
    


def smcify_arg(arg):
    """
    Rebuild an ARG so that is follows the SMC assumptions
    """
    
    arg2 = arg.get_marginal_tree(-.5)
    remove_single_lineages(arg2)
    sprs = iter_arg_sprs(arg, use_leaves=True)
    make_arg_from_sprs(arg2, sprs, ignore_self=True)

    return arg2


def has_self_cycles(arg):
    """
    Return True if there are lineages that coalesce with themselves

    Requires ancestral sequences set.
    """

    # Such a cycle does not contain 'local coalescent nodes' on the sides
    # but it might have non-local coalescent and recombination nodes.
    # The relative order of theses nodes from the left and right side
    # does not matter
    #
    #      |
    #     coal
    #    /   \
    #    |   |
    #    \   /
    #    recomb
    #      |
    
    # get overall postorder
    # assumes stable sort
    nodes = list(arg.postorder())
    nodes.sort(key=lambda x: x.age)
    order = dict((x, i) for i, x in enumerate(nodes))    

    # find cycles by their recombination nodes
    recombs = [x.name for x in arg if x.event == "recomb"]

    # find smallest separation
    recomb_pos = [arg[x].pos for x in recombs]
    recomb_pos.sort()
    eps = .5
    for i in xrange(1, len(recomb_pos)):
        sep = recomb_pos[i] - recomb_pos[i-1]
        if sep > 0 and sep/2.0 < eps:
            eps = sep / 2.0
    
    for recomb_name in recombs:
        if recomb_name not in arg:
            continue
        if is_self_cycle(arg, arg[recomb_name], order=order, eps=eps):
            return True

    return False


def iter_self_cycles(arg):
    """
    Return True if there are lineages that coalesce with themselves

    Requires ancestral sequences set.
    """

    # Such a cycle does not contain 'local coalescent nodes' on the sides
    # but it might have non-local coalescent and recombination nodes.
    # The relative order of theses nodes from the left and right side
    # does not matter
    #
    #      |
    #     coal
    #    /   \
    #    |   |
    #    \   /
    #    recomb
    #      |
    
    # get overall postorder
    # assumes stable sort
    nodes = list(arg.postorder())
    nodes.sort(key=lambda x: x.age)
    order = dict((x, i) for i, x in enumerate(nodes))    

    # find cycles by their recombination nodes
    recombs = [x.name for x in arg if x.event == "recomb"]

    # find smallest separation
    recomb_pos = [arg[x].pos for x in recombs]
    recomb_pos.sort()
    eps = .5
    for i in xrange(1, len(recomb_pos)):
        sep = recomb_pos[i] - recomb_pos[i-1]
        if sep > 0 and sep/2.0 < eps:
            eps = sep / 2.0
    
    for recomb_name in recombs:
        if recomb_name not in arg:
            continue
        if is_self_cycle(arg, arg[recomb_name], order=order, eps=eps):
            yield arg[recomb_name]




def is_self_cycle(arg, recomb, order=None, eps=1e-4):

    def is_local_coal(node, child, pos):
        if node.event != "coal":
            return False

        i = node.children.index(child)
        other_child = node.children[1 - i]

        for start, end in other_child.data["ancestral"]:
            if start < pos < end:
                return True

        return False

    if order is None:
        # get overall postorder
        # assumes stable sort
        nodes = list(arg.postorder())
        nodes.sort(key=lambda x: x.age)
        order = dict((x, i) for i, x in enumerate(nodes))    

    # find cycle
    # also check for local coal nodes along the way
    rpos = recomb.pos
    path1 = []
    path2 = []
    ptr1 = arg.get_local_parent(recomb, rpos-eps)
    ptr2 = arg.get_local_parent(recomb, rpos+eps)
    while ptr1 and ptr2:
        order1 = order[ptr1]
        order2 = order[ptr2]

        if order1 < order2:
            if is_local_coal(ptr1,path1[-1] if path1 else recomb,rpos-eps):
                break
            path1.append(ptr1)
            ptr1 = arg.get_local_parent(ptr1, rpos-eps)

        elif order1 > order2:
            if is_local_coal(ptr2,path2[-1] if path2 else recomb,rpos+eps):
                break
            path2.append(ptr2)
            ptr2 = arg.get_local_parent(ptr2, rpos+eps)

        else:
            # we have reached coal node
            assert ptr1 == ptr2
            coal = ptr1
            return True

    return False
    

def remove_self_cycles(arg):
    """
    Removes cycles that represent a lineage coalescing with itself

    Requires ancestral sequences set.
    """

    # Such a cycle does not contain 'local coalescent nodes' on the sides
    # but it might have non-local coalescent and recombination nodes.
    # The relative order of theses nodes from the left and right side
    # does not matter
    #
    #      |
    #     coal
    #    /   \
    #    |   |
    #    \   /
    #    recomb
    #      |

    def is_local_coal(node, child, pos):
        if node.event != "coal":
            return False

        i = node.children.index(child)
        other_child = node.children[1 - i]

        for start, end in other_child.data["ancestral"]:
            if start < pos < end:
                return True

        return False

    # get overall postorder
    # assumes stable sort
    nodes = list(arg.postorder())
    nodes.sort(key=lambda x: x.age)
    order = dict((x, i) for i, x in enumerate(nodes))    

    # find cycles by their recombination nodes
    recombs = [x.name for x in arg if x.event == "recomb"]

    # find smallest separation
    recomb_pos = [arg[x].pos for x in recombs]
    recomb_pos.sort()
    eps = .5
    for i in xrange(1, len(recomb_pos)):
        sep = recomb_pos[i] - recomb_pos[i-1]
        if sep > 0 and sep/2.0 < eps:
            eps = sep / 2.0
            
    
    for recomb_name in recombs:
        if recomb_name not in arg:
            continue
        recomb = arg[recomb_name]
        rpos = recomb.pos
        
        # find cycle
        # also check for local coal nodes along the way
        is_cycle = False
        path1 = []
        path2 = []
        ptr1 = arg.get_local_parent(recomb, rpos-eps)
        ptr2 = arg.get_local_parent(recomb, rpos+eps)
        while ptr1 and ptr2:
            order1 = order[ptr1]
            order2 = order[ptr2]
            
            if order1 < order2:
                if is_local_coal(ptr1,path1[-1] if path1 else recomb,rpos-eps):
                    break
                path1.append(ptr1)
                ptr1 = arg.get_local_parent(ptr1, rpos-eps)
                
            elif order1 > order2:
                if is_local_coal(ptr2,path2[-1] if path2 else recomb,rpos+eps):
                    break
                path2.append(ptr2)
                ptr2 = arg.get_local_parent(ptr2, rpos+eps)
                
            else:
                # we have reached coal node
                assert ptr1 == ptr2
                coal = ptr1
                is_cycle = True
                break

        if not is_cycle:
            # this recombination node is not a cycle
            # either because it contains a local coal node or never recoals
            # which can happen in SMC ARGs
            continue

        if path1:
            assert coal in path1[-1].parents
        else:
            assert coal in recomb.parents
        if path2:
            assert coal in path2[-1].parents
        else:
            assert coal in recomb.parents

        if len(set(path1) & set(path2)) != 0:
            print [(order[x], x) for x in path1]
            print [(order[x], x) for x in path2]
            assert False
        
        # remove coal node
        top = coal.parents[0] if coal.parents else None
        if top:
            util.replace(top.children, coal, None)

        # remove recomb node
        bottom = recomb.children[0]
        util.replace(bottom.parents, recomb, None)

        # unlink nodes in left path
        last = recomb
        for node in path1:
            util.replace(last.parents, node, None)
            util.replace(node.children, last, None)
            last = node
        util.replace(last.parents, coal, None)

        # unlink nodes in right path
        last = recomb
        #print "--"
        for node in path2:
            #print last.parents, node
            util.replace(last.parents, node, None)
            util.replace(node.children, last, None)
            last = node
        util.replace(last.parents, coal, None)
        
        
        # merge paths
        combine = path1 + path2
        combine.sort(key=lambda x: (x.age, order[x]))
        last = bottom
        
        for n in combine:
            util.replace(last.parents, None, n)
            util.replace(n.children, None, last)
            last = n
            
        if top:
            util.replace(last.parents, None, top)
            util.replace(top.children, None, last)
        else:
            # no top node
            if last.event == "coal" or last.event == "gene":
                last.parents = []
            elif last.event == "recomb":
                # remove last recomb node since it is a single lineage
                c = last.children[0]
                p = last.parents[1 - last.parents.index(None)]
                util.replace(c.parents, last, p)
                util.replace(p.children, last, c)
                del arg.nodes[last.name]
            else:
                raise Exception("unknown event '%s'" % last.event)

        del arg.nodes[recomb.name]
        del arg.nodes[coal.name]
            




def remove_self_cycles2(arg):
    """
    Removes cycles that represent a lineage coalescing with itself

    Requires ancestral sequences set.
    """

    # Such a cycle does not contain 'local coalescent nodes' on the sides
    # but it might have non-local coalescent and recombination nodes.
    # The relative order of theses nodes from the left and right side
    # does not matter
    #
    #      |
    #     coal
    #    /   \
    #    |   |
    #    \   /
    #    recomb
    #      |

    # get overall postorder
    # assumes stable sort
    nodes = list(arg.postorder())
    nodes.sort(key=lambda x: x.age)
    order = dict((x, i) for i, x in enumerate(nodes))    

    # find cycles by their recombination nodes
    recombs = [x.name for x in arg if x.event == "recomb"]

    # find smallest separation
    recomb_pos = [arg[x].pos for x in recombs]
    recomb_pos.sort()
    eps = .5
    for i in xrange(1, len(recomb_pos)):
        sep = recomb_pos[i] - recomb_pos[i-1]
        if sep > 0 and sep/2.0 < eps:
            eps = sep / 2.0
            
    
    for recomb_name in recombs:
        if recomb_name not in arg:
            continue
        recomb = arg[recomb_name]
        rpos = recomb.pos
        
        # find cycle
        # also check for local coal nodes along the way
        is_cycle = False
        path1 = []
        path2 = []
        ptr1 = arg.get_local_parent(recomb, rpos-eps)
        ptr2 = arg.get_local_parent(recomb, rpos+eps)
        while ptr1 and ptr2:
            order1 = order[ptr1]
            order2 = order[ptr2]
            
            if order1 < order2:
                if ptr1.event == "coal":
                    break
                path1.append(ptr1)
                ptr1 = arg.get_local_parent(ptr1, rpos-eps)
                
            elif order1 > order2:
                if ptr2.event == "coal":
                    break
                path2.append(ptr2)
                ptr2 = arg.get_local_parent(ptr2, rpos+eps)
                
            else:
                # we have reached coal node
                assert ptr1 == ptr2
                coal = ptr1
                is_cycle = True
                break

        if not is_cycle:
            # this recombination node is not a cycle
            # either because it contains a local coal node or never recoals
            # which can happen in SMC ARGs
            continue
        
        # remove coal node
        top = coal.parents[0] if coal.parents else None
        if top:
            util.replace(top.children, coal, None)

        # remove recomb node
        bottom = recomb.children[0]
        util.replace(bottom.parents, recomb, None)

        # unlink nodes in left path
        last = recomb
        for node in path1:
            util.replace(last.parents, node, None)
            util.replace(node.children, last, None)
            last = node
        util.replace(last.parents, coal, None)

        # unlink nodes in right path
        last = recomb
        for node in path2:
            util.replace(last.parents, node, None)
            util.replace(node.children, last, None)
            last = node
        util.replace(last.parents, coal, None)
        
        
        # merge paths
        combine = path1 + path2
        combine.sort(key=lambda x: (x.age, order[x]))
        last = bottom
        
        for n in combine:
            util.replace(last.parents, None, n)
            util.replace(n.children, None, last)
            last = n
            
        if top:
            util.replace(last.parents, None, top)
            util.replace(top.children, None, last)
        else:
            # no top node
            if last.event == "coal":
                last.parents = []
            elif last.event == "recomb":
                # remove last recomb node since it is a single lineage
                c = last.children[0]
                p = last.parents[1 - last.parents.index(None)]
                util.replace(c.parents, last, p)
                util.replace(p.children, last, c)
                del arg.nodes[last.name]
            else:
                raise Exception("unknown event '%s'" % node.event)

        del arg.nodes[recomb.name]
        del arg.nodes[coal.name]


    

def subarg_by_leaf_names(arg, leaf_names, keep_single=False):
    """
    Removes any leaf from the arg that is not in leaf name set
    """

    return subarg_by_leaves(arg, [arg[x] for x in leaf_names],
                            keep_single=keep_single)


def arg_lca(arg, leaves, pos, time=None, local=None):
    """
    Find the Least Common Ancestor (LCA) of a set of leaves in the ARG

    arg    -- an ARG
    leaves -- a list of nodes in arg
    pos    -- position along sequence to perform LCA
    time   -- the time ascend to (optional)
    local  -- the set of nodes considered local (optional)
    """

    def is_local_coal(arg, node, pos, local):
        return (len(node.children) == 2 and
                node.children[0] in local and
                arg.get_local_parent(node.children[0], pos) == node and
                node.children[1] in local and
                arg.get_local_parent(node.children[1], pos) == node and
                node.children[0] != node.children[1])


    order = dict((node, i) for i, node in enumerate(
        arg.postorder_marginal_tree(pos)))
    if local is None:
        local = order
    
    queue = [(order[arg[x]], arg[x]) for x in leaves]
    seen = set(x[1] for x in queue)
    heapq.heapify(queue)

    while len(queue) > 1:
        i, node = heapq.heappop(queue)
        parent = arg.get_local_parent(node, pos)
        if parent and parent not in seen:
            seen.add(parent)
            heapq.heappush(queue, (order[parent], parent))
    node = queue[0][1]
    parent = arg.get_local_parent(node, pos)


    if time is not None:
        while parent and parent.age <= time:
            if is_local_coal(arg, parent, pos, local):
                break            
            node = parent
            parent = arg.get_local_parent(node, pos)

        # safety check
        if parent and parent.age < time:
            print (pos, leaves, parent.age, time)
            tree = arg.get_marginal_tree(pos).get_tree()
            tree.write()
            treelib.draw_tree_names(tree, maxlen=8, minlen=8)
            assert False

    return node



#=============================================================================
# region functions


def split_regions(pos, side, regions):
    """
    Iterates through the regions on the left (side=0) or right (side=1) of 'pos'
    """
    
    for reg in regions:
        if side == 0:
            if reg[1] <= pos:
                # keep all regions fully left of recomb position
                yield reg
            elif reg[0] < pos:
                # cut region
                yield (reg[0], pos)
        elif side == 1:
            if reg[0] >= pos:
                # keep all regions fully right of recomb position
                yield reg
            elif reg[1] > pos:
                # cut region
                yield (pos, reg[1])
        else:
            raise Exception("side not specified")


def count_region_overlaps(*region_sets):
    """
    Count how many regions overlap each interval (start, end)
    
    Iterates through (start, end, count) sorted
    """

    # build endpoints list
    end_points = []
    for regions in region_sets:
        for reg in regions:
            end_points.append((reg[0], 0))
            end_points.append((reg[1], 1))
    end_points.sort()

    count = 0
    start = None
    end = None
    last = None
    for pos, kind in end_points:
        if last is not None and pos != last:
            yield last, pos, count
        if kind == 0:
            count += 1
        elif kind == 1:
            count -= 1
        last = pos

    if last is not None and pos != last:
        yield last, pos, count
    

        
def groupby_overlaps(regions, bygroup=True):
    """
    Group ranges into overlapping groups
    Ranges must be sorted by start positions
    """

    start = -util.INF
    end = -util.INF
    group = None
    groupnum = -1
    for reg in regions:        
        if reg[0] > end:
            # start new group
            start, end = reg
            groupnum += 1

            if bygroup:
                if group is not None:
                    yield group
                group = [reg]
            else:
                yield (groupnum, reg)

        else:
            # append to current group
            if reg[1] > end:
                end = reg[1]

            if bygroup:
                group.append(reg)
            else:
                yield (groupnum, reg)

    if bygroup and group is not None and len(group) > 0:
        yield group


#=============================================================================
# mutations and splits


def sample_arg_mutations(arg, mu):
    """
    mu -- mutation rate (mutations/site/gen)
    """

    mutations = []

    locsize = arg.end - arg.start
    u = mu * locsize

    for node in arg:
        for parent in node.parents:
            for region in arg.get_ancestral(node, parent=parent):
                # ensure node is not MRCA
                for pregion in parent.data["ancestral"]:
                    if pregion[0] <= region[0] < pregion[1]:
                        break
                else:
                    continue
                
                frac = (region[1] - region[0]) / locsize
                dist = parent.age - node.age
                t = parent.age
                while True:
                    t -= random.expovariate(u * frac)
                    if t < node.age:
                        break
                    pos = random.uniform(region[0], region[1])
                    mutations.append((node, parent, pos, t))

    return mutations


def get_marginal_leaves(arg, node, pos):
    return (x for x in arg.preorder_marginal_tree(pos, node) if x.is_leaf())


def get_mutation_split(arg, mutation):
    node, parent, pos, t = mutation
    return tuple(sorted(x.name for x in get_marginal_leaves(arg, node, pos)))


def split_to_tree_branch(tree, split):
    
    node = treelib.lca([tree[name] for name in split])

    if sorted(split) != sorted(node.leaf_names()):
        inv = [x for x in tree.leaf_names() if x not in split]
        node = treelib.lca([tree[name] for name in inv])
        if sorted(inv) != sorted(node.leaf_names()):
            # split does not conform to tree
            return None

    return node


def split_to_arg_branch(arg, pos, split):

    # TODO: make more efficient
    tree = arg.get_tree(pos)
    node = split_to_tree_branch(tree, split)
    if node is not None:
        return arg[node.name]
    else:
        None


def iter_tree_splits(tree):
    for node in tree:
        if len(node.children) != 2 or node.children[0] == node.children[1]:
            continue
        split = tuple(sorted(tree.leaf_names(node)))
        if len(split) > 1:
            yield split
    

def is_split_compatible(split1, split2):
    """Returns True if two splits are compatible"""
    i = j = 0
    intersect = 0
    while i < len(split1) and j < len(split2):
        if split1[i] == split2[j]:
            intersect += 1
            i += 1
            j += 1
        elif split1[i] < split2[j]:
            i += 1
        else:
            j += 1
    #intersect = len(split1 & split2)
    return intersect == 0 or intersect == min(len(split1), len(split2))


def is_split_compatible_unpolar2(split1, split2, leaves):

    a = set(split1)
    b = set(split2)
    x00 = False
    x01 = False
    x10 = False
    x11 = False
    for l in leaves:
        if l in a:
            if l in b:
                x11 = True
            else:
                x10 = True
        else:
            if l in b:
                x01 = True
            else:
                x00 = True

    return not (x00 and x01 and x10 and x11)

def is_split_compatible_unpolar(split1, split2, leaves):
    if is_split_compatible(split1, split2):
        return True
    else:
        split1rev = tuple(x for x in leaves if x not in split1)
        return is_split_compatible(split1rev, split2)


def split_relation(split1, split2):

    i = j = 0
    intersect = 0
    while i < len(split1) and j < len(split2):
        if split1[i] == split2[j]:
            intersect += 1
            i += 1
            j += 1
        elif split1[i] < split2[j]:
            i += 1
        else:
            j += 1

    if intersect == 0:
        return "disjoint"

    elif intersect == len(split1):
        if intersect == len(split2):
            assert split1 == split2
            return "equal"
        else:
            return "child"

    elif intersect == len(split2):
        return "parent"

    else:
        return "conflict"
    
    return intersect == 0 or intersect == min(len(split1), len(split2))


def iter_mutation_splits(arg, mutations):

    nleaves = sum(1 for x in arg.leaves())

    for node, parent, pos, t in mutations:
        split = tuple(sorted(x.name for x in get_marginal_leaves(
            arg, node, pos)))
        if len(split) != 1 and len(split) != nleaves:
            yield pos, split






#=============================================================================
# alignments

def make_alignment(arg, mutations, ancestral="A", derived="C"):
    aln = fasta.FastaDict()
    alnlen = int(arg.end - arg.start)
    leaves = list(arg.leaf_names())
    nleaves = len(leaves)

    # sort mutations by position
    mutations.sort(key=lambda x: x[2])

    # make align matrix
    mat = []
    
    pos = arg.start
    muti = 0
    for i in xrange(alnlen):
        if muti >= len(mutations) or i < int(mutations[muti][2]):
            # no mut
            mat.append(ancestral * nleaves)
        else:
            # mut
            node, parent, mpos, t = mutations[muti]
            row = []
            split = set(x.name for x in get_marginal_leaves(arg, node, mpos))
            mat.append("".join((derived if leaf in split else ancestral)
                               for leaf in leaves))
            muti += 1
    
    # make fasta
    for i, leaf in enumerate(leaves):
        aln[leaf] = "".join(x[i] for x in mat)

    return aln


#=============================================================================
# input/output


def write_arg(filename, arg):
    out = util.open_stream(filename, "w")

    # write ARG key values
    out.write("start=%s\tend=%s\n" % (arg.start, arg.end))

    # write nodes header
    out.write("\t".join(("name", "event", "age", "pos", "parents", "children"))
              + "\n")

    # write nodes
    for node in arg:
        util.print_row(
            node.name, node.event, node.age, node.pos,
            ",".join(str(x.name) for x in node.parents),
            ",".join(str(x.name) for x in node.children),
            out=out)
        
    if isinstance(filename, basestring):
        out.close()


def parse_number(text):
    if text.isdigit():
        return int(text)
    else:
        return float(text)


def parse_node_name(text):
    if text.isdigit():
        return int(text)
    else:
        return text

def parse_key_value(field):
    try:
        i = field.index("=")
        return field[:i], field[i+1:]
    except:
        raise Exception("improper key-value field '%s'" % text)
    


def read_arg(filename, arg=None):
    infile = util.DelimReader(filename)

    if arg is None:
        arg = ARG()

    # read ARG key values
    row = infile.next()
    for field in row:
        key, val = parse_key_value(field)
        if key == "start":
            arg.start = int(val)
        elif key == "end":
            arg.end = int(val)
    
    # read header
    row = infile.next()
    assert row == ["name", "event", "age", "pos", "parents", "children"]

    # read nodes
    clinks = {}
    plinks = {}
    for row in infile:
        node = arg.new_node(name=parse_node_name(row[0]), event=row[1], 
                            age=float(row[2]), 
                            pos=parse_number(row[3]))
        if len(row) > 4 and len(row[4]) > 0:
            plinks[node.name] = map(parse_node_name, row[4].split(","))
        if len(row) > 5 and len(row[5]) > 0:
            clinks[node.name] = map(parse_node_name, row[5].split(","))

    # setup parents
    for node in arg:
        for parent_name in plinks.get(node.name, ()):
            parent = arg.nodes.get(parent_name)
            if parent:
                node.parents.append(parent)
            else:
                raise Exception("node '%s' has unknown parent '%s'" % 
                                (node.name, parent_name))

        # detect root
        if parent_name not in plinks:
            arg.root = node

    # setup children
    for node in arg:
        for child_name in clinks.get(node.name, ()):
            child = arg.nodes.get(child_name)
            if child:
                node.children.append(child)
                assert node in child.parents, \
                    "node '%s' doesn't have parent '%s' (%s)" % (
                    child.name, node.name, str(child.parents))
            else:
                raise Exception("node '%s' has unknown child '%s'" % 
                                (node.name, child_name))


    # set nextname
    for name in arg.nodes:
        if isinstance(name, int):
            arg.nextname = max(arg.nextname, name+1)

    return arg


def write_tree_tracks(filename, arg, start=None, end=None, verbose=False):
    out = util.open_stream(filename, "w")
    for block, tree in iter_tree_tracks(arg, start, end):
        if verbose:
            print >>sys.stderr, "writing block", block
        remove_single_lineages(tree)
        tree = tree.get_tree()
        out.write(str(int(block[0]))+"\t"+str(int(block[1]))+"\t")
        tree.write(out, oneline=True)
        out.write("\n")
    if isinstance(filename, basestring):
        out.close()


def read_tree_tracks(filename):
    for row in util.DelimReader(filename):
        yield (int(row[0]), int(row[1])), treelib.parse_newick(row[2])



def write_mutations(filename, arg, mutations):
    out = util.open_stream(filename, "w")

    for mut in mutations:
        l = get_marginal_leaves(arg, mut[0], mut[2])
        util.print_row(mut[2], mut[3], ",".join(x.name for x in l), out=out)

    if isinstance(filename, basestring):
        out.close()


def read_mutations(filename):
    for row in util.DelimReader(filename):
        chroms = row[2].split(",") if row[2] else []
        yield int(row[0]), float(row[1]), chroms


def write_ancestral(filename, arg):
    out = util.open_stream(filename, "w")

    for node in arg:
        regions = util.flatten(node.data.get("ancestral", ()))
        util.print_row(node.name, *regions, out=out)

    if isinstance(filename, basestring):
        out.close()


def read_ancestral(filename, arg):
    for row in util.DelimReader(filename):
        node = arg[parse_node_name(row[0])]
        node.data["ancestral"] = [(int(row[i]), int(row[i+1]))
                                  for i in xrange(1, len(row), 2)]


#=============================================================================
# OLD CODE



def sample_mutations(arg, u):
    """
    u -- mutation rate (mutations/locus/gen)

    DEPRECATED: use sample_arg_mutations() instead
    """

    mutations = []

    locsize = arg.end - arg.start

    for node in arg:
        for parent in node.parents:
            for region in arg.get_ancestral(node, parent=parent):
                # ensure node is not MRCA
                for pregion in parent.data["ancestral"]:
                    if pregion[0] <= region[0] < pregion[1]:
                        break
                else:
                    continue
                
                frac = (region[1] - region[0]) / locsize
                dist = parent.age - node.age
                t = parent.age
                while True:
                    t -= random.expovariate(u * frac)
                    if t < node.age:
                        break
                    pos = random.uniform(region[0], region[1])
                    mutations.append((node, parent, pos, t))

    return mutations


'''
SLOW remove cycles

def remove_self_cycles(arg, eps=.5):
    """
    Removes cycles that represent a lineage coalescing with itself

    Requires ancestral sequences set.
    """

    # Such a cycle does not contain 'local coalescent nodes' on the sides
    # but it might have non-local coalescent and recombination nodes.
    # The relative order of theses nodes from the left and right side
    # does not matter
    #
    #      |
    #     coal
    #    /   \
    #    |   |
    #    \   /
    #    recomb
    #      |

    def is_local_coal(node, child, pos):
        if node.event != "coal":
            return False
        return True

        i = node.children.index(child)
        other_child = node.children[1 - i]

        for start, end in other_child.data["ancestral"]:
            #print node, other_child, (start, end), pos
            if start < pos < end:
                return True

        return False

    # get overall postorder
    order = dict((x, i) for i, x in enumerate(arg.postorder()))


    # find cycles by their recombination nodes
    recombs = [x.name for x in arg if x.event == "recomb"]
    for recomb_name in recombs:
        if recomb_name not in arg:
            continue
        recomb = arg[recomb_name]
        rpos = recomb.pos

        # get left path
        path1 = []
        ptr = arg.get_local_parent(recomb, rpos-eps)
        while ptr:
            path1.append(ptr)
            ptr = arg.get_local_parent(ptr, rpos-eps)

        # get right path
        path2 = []
        ptr = arg.get_local_parent(recomb, rpos+eps)
        while ptr:
            path2.append(ptr)
            ptr = arg.get_local_parent(ptr, rpos+eps)
                            

        # find recoal node
        i = -1
        length = min(len(path1), len(path2))
        while -i <= length and path1[i] == path2[i]:
            i -= 1
        if i == -1:
            # this happens with SMC ARGs
            continue
        a = len(path1) + (i + 1)
        b = len(path2) + (i + 1)
        coal = path1[a]
        

        # are there any coal nodes in left and right paths?
        is_cycle = True
        for i in range(a):
            if is_local_coal(
                path1[i], path1[i-1] if i > 0 else recomb, rpos-eps):
                is_cycle = False
                break

        for i in range(b):
            if is_local_coal(
                path2[i], path2[i-1] if i > 0 else recomb, rpos+eps):
                is_cycle = False
                break

        if not is_cycle:
            # this recombination node is not a cycle
            print "recomb", recomb, "pos=", rpos, "is not a cycle"
            continue

        
        print path1, path2, recomb, coal, i, a, b
        print path1[:a], path2[:b]


        # remove coal node
        top = coal.parents[0] if coal.parents else None
        if top:
            util.replace(top.children, coal, None)

        # remove recomb node
        bottom = recomb.children[0]
        util.replace(bottom.parents, recomb, None)

        # unlink nodes in left path
        last = recomb
        for node in path1[:a]:
            util.replace(last.parents, node, None)
            util.replace(node.children, last, None)
            last = node
        util.replace(last.parents, coal, None)

        # unlink nodes in right path
        last = recomb
        for node in path2[:b]:
            util.replace(last.parents, node, None)
            util.replace(node.children, last, None)
            last = node
        util.replace(last.parents, coal, None)
        
        
        # merge paths
        combine = path1[:a] + path2[:b]
        #print "path1", [order[x] for x in path1]
        #print "path2", [order[x] for x in path2]
        combine.sort(key=lambda x: (x.age, order[x]))
        #print [(x.name, x.age, order[x]) for x in combine]

        
        last = bottom
        for n in combine:
            util.replace(last.parents, None, n)
            util.replace(n.children, None, last)
            last = n
        if top:
            util.replace(last.parents, None, top)
            util.replace(top.children, None, last)
        else:
            print "no top for", recomb, rpos
            
            # no top node
            if last.event == "coal":
                last.parents = []
            elif last.event == "recomb":
                # remove last recomb node since it is a single lineage
                c = last.children[0]
                p = last.parents[1 - last.parents.index(None)]
                util.replace(c.parents, last, p)
                util.replace(p.children, last, c)
                del arg.nodes[last.name]
                print "remove last", last.name
            else:
                raise Exception("unknown event '%s'" % node.event)

        del arg.nodes[recomb.name]
        del arg.nodes[coal.name]
        print "remove", recomb.name, coal.name
        print "  ", order[bottom], order[recomb], order[coal], order[top]

        assert_arg(arg)
        #print "good"
            


'''


'''

def iter_marginal_trees2(arg, start=None, end=None, visible=False):
    """
    Iterate over the marginal trees of an ARG
    """
    
    for a,b in iter_recomb_blocks(arg, start, end, visible=visible):
        yield arg.get_marginal_tree((a+b) / 2.0)



def iter_tree_tracks2(arg, start=None, end=None, visible=False):
    for a, b in iter_recomb_blocks(arg, start, end, visible=visible):
        tree = arg.get_marginal_tree((a+b) / 2.0)
        remove_single_lineages(tree)
        yield (a, b), tree.get_tree()


'''
'''
# TODO: this mostly implements Hein's sampling method, just need to be
# cleaned up and should remove bgsel features since its not exactly
# acturate anyways.

def sample_arg_bgsel(k, ns, rho, start=0.0, end=0.0, t=0):
    """
    ns = [[n1, end1], [n2, end2], ...]
    """

    # TODO: convert this into simply a SMC

    arg = ARG(start, end)
    ni = 0  # index into ns (list popsizes)
    pos = 0 # position simulated so far

    # sample initial coal tree with no recomb
    times, events = sample_coal_recomb_times(k, ns[ni][0], 0)
    lineages = set(arg.new_node() for i in xrange(k))
    nlineages = []
    for t in times:
        a, b = random.sample(lineages, 2)
        parent = arg.new_node(age=t, children=[a,b], event="coal")
        a.parents.append(parent)
        b.parents.append(parent)
        lineages.remove(a)
        lineages.remove(b)
        lineages.add(parent)
        nlineages.append(len(lineages))
    oldest = parent

    nrecomb = 0

    while True:
        # find correct ni
        while ns[ni][1] < pos:
            ni += 1

        # sample next recomb
        totlen = sum(arg.get_local_dist(x, pos+1)
                     for x in arg.postorder_marginal_tree(pos+1))
        recomb_per_site = rho * totlen;

        # sample recomb pos
        pos += random.expovariate(recomb_per_site)

        if pos > end:
            # no more recombinations
            break

        # sample recomb node
        nodes = list(arg.postorder_marginal_tree(pos+1))
        lens = [arg.get_local_dist(x, pos+1) for x in nodes]
        child = nodes[stats.sample(lens)]

        # make recomb node and insert into arg
        old_parent = arg.get_local_parent(child, pos+1)
        recomb_node = arg.new_node(age=child.age + random.random() *
                                   arg.get_local_dist(child, pos+1),
                                   children=[child],
                                   parents=[old_parent],
                                   pos=pos, event="recomb")
        #old_parent = child.parents[-1]
        i = old_parent.children.index(child)
        old_parent.children[i] = recomb_node
        child.parents[child.parents.index(old_parent)] = recomb_node

        # add recomb_node age into times
        for i in xrange(len(times)):
            if recomb_node.age < times[i]:                
                times.insert(i, recomb_node.age)
                if i == 0:
                    nlineages.insert(i, k+1)
                else:
                    nlineages.insert(i, nlineages[i-1]+1)
                break
        else:
            assert False

        # coalesce recomb_node back into ARG
        for i in xrange(i, len(times) - 1):
            coal_time = times[i] + random.expovariate(
                (nlineages[i]-1) / ns[ni][0])
            if coal_time < times[i+1]:
                # coal event, choose sister
                candidates = [(node, parent) for node in arg
                              for parent in node.parents
                              if node.age < coal_time < parent.age]
                sister, sister_parent = random.sample(candidates, 1)[0]
                times.insert(i+1, coal_time)
                nlineages.insert(i+1, nlineages[i]-1)
                break
            else:
                # no coal, keep moving up
                nlineages[i+1] += 1
        else:
            # no coal at all choose oldest sister
            i += 1
            sister, sister_parent = oldest, None
            coal_time = times[i] + random.expovariate(
                (nlineages[i]-1) / 2.0 / ns[ni][0])
            times.append(coal_time)
            nlineages.append(1)

        print pos, ni, ns[ni][0]

        # make coal node
        coal_node = arg.new_node(age=coal_time, event="coal",
                                 children=[sister, recomb_node])
        recomb_node.parents.append(coal_node)
        if sister_parent:
            i = sister_parent.children.index(sister)
            sister_parent.children[i] = coal_node
            i = sister.parents.index(sister_parent)
            sister.parents[i] = coal_node
            coal_node.parents.append(sister_parent)
        else:
            sister.parents.append(coal_node)
            oldest = coal_node

        nrecomb += 1

    arg.root = oldest

    return arg
'''
    
