"""
    arglib.py
    
    Ancestral recombination graph (ARG)

"""


#=============================================================================
# imports

from __future__ import division

# python libs
import random
from itertools import izip
from collections import defaultdict
import heapq

# rasmus libs
from rasmus import treelib, util



#=============================================================================
# Ancestral Reconstruction Graph

class CoalNode (object):

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

    def get_dist(self, parent_index):
        return self.parents[parent_index].age - self.age

    def get_dists(self):
        return [p.age - self.age for p in self.parents]

    def copy(self):
        node = CoalNode(self.name, age=self.age, event=self.event,
                        pos=self.pos)
        node.data = dict(self.data)
        return node

    def is_leaf(self):
        return len(self.children) == 0
    



class ARG (object):

    def __init__(self, start=0.0, end=1.0):
        self.root = None
        self.nodes = {}
        self.nextname = 1
        self.start = start
        self.end = end


    def __iter__(self):
        return self.nodes.itervalues()


    def __len__(self):
        """Returns number of nodes in tree"""
        return len(self.nodes)


    def __getitem__(self, key):
        """Returns node by name"""
        return self.nodes[key]


    def __setitem__(self, key, node):
        """Adds a node to the tree"""
        node.name = key
        self.add(node)


    def new_name(self):
        name = self.nextname
        self.nextname += 1
        return name

    def add(self, node):
        self.nodes[node.name] = node
        return node

    def remove(self, node):
        for child in node.children:
            child.parents.remove(node)
        for parent in node.parents:
            parent.children.remove(node)
        del self.nodes[node.name]


    def leaves(self):
        """
        Iterates over the leaves of the ARG
        """
        for node in self:
            if len(node.children) == 0:
                yield node

    def postorder(self):
        """
        Iterates through nodes in postorder traversal
        """

        visit = defaultdict(lambda: 0)
        queue = list(self.leaves())

        for node in queue:
            yield node
            for parent in node.parents:
                visit[parent] += 1

                # if all children of parent has been visited then queue parent
                if visit[parent] == len(parent.children):
                    queue.append(parent)
        
                    

    #==============================

    def set_recomb_pos(self, start=None, end=None, descrete=False):
        """
        Set all recombination positions in the ARG
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

        for node in self.postorder():
            if node.is_leaf():
                # initialize leaves with entire extant sequence
                node.data["ancestral"] = [(self.start, self.end)]
            else:
                if node.event == "coal":
                    # union of ancestral of children
                    # get all child regions
                    regions = [region for child in node.children
                               for region in
                               self.get_ancestral(child, parent=node)]
                    regions.sort()

                    node.data["ancestral"] = [
                        (min(r[0] for r in group), max(r[1] for r in group))
                        for group in groupby_overlaps(regions)]

                elif node.event == "recomb":
                    # inherit all ancestral
                    assert len(node.children) == 1
                    node.data["ancestral"] = self.get_ancestral(
                        node.children[0], parent=node)

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
                    if reg[1] < node.pos:
                        # keep all regions fully left of recomb position
                        regions.append(reg)
                    elif reg[0] < node.pos:
                        # cut region
                        regions.append((reg[0], node.pos))
                elif side == 1:
                    if reg[0] > node.pos:
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

                            
                
    def get_marginal_tree(self, pos):
        """
        Returns the marginal tree of the ARG containing position 'pos'
        """

        # make new ARG to contain marginal tree
        tree = ARG()

        # populate tree with marginal nodes
        for node in self.postorder_marginal_tree(pos):
            tree.add(node.copy())
        
        # set parent and children
        for node2 in tree:
            node = self[node2.name]
            for parent in node.parents:
                if parent.name in tree.nodes:
                    parent2 = tree[parent.name]
                    node2.parents = [parent2]
                    parent2.children.append(node2)
                    break
            if len(node2.parents) == 0:
                tree.root = node2
        
        return tree


    def postorder_marginal_tree(self, pos):
        """
        Iterate postorder over the nodes in the marginal tree at position 'pos'
        """

        # initialize heap
        heap = [(node.age, node) for node in self.leaves()]
        seen = set([None])
        
        # add all ancestor of lineages
        while len(heap) > 0:
            age, node = heapq.heappop(heap)
            yield node

            # find correct marginal parent
            # add parent to lineages if it has not been seen before
            parent = self.get_local_parent(node, pos)
            if parent not in seen:
                heapq.heappush(heap, (parent.age, parent))
                seen.add(parent)



    def get_local_parent(self, node, pos):
        """Return the local parent of 'node' for position 'pos'"""
        
        if (node.event == "gene" or node.event == "coal"):
            if len(node.parents) > 0:
                return node.parents[0]
            else:
                return None
        elif node.event == "recomb":
            return node.parents[0 if pos < node.pos else 1]
        else:
            raise Exception("unknown event '%s'" % node.event)
        
        
    def get_tree(self, pos=None):
        """
        Returns a treelib.Tree() object representing the ARG if it is a tree

        if 'pos' is given, return a treelib.Tree() for the marginal tree at
        position 'pos'.
        """

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
            

    def prune(self, remove_single=True):
        """
        Prune ARG to only those nodes with ancestral sequence
        """

        prune = set(node for node in self if len(node.data["ancestral"]) == 0)

        # remove pruned nodes
        for node in list(self):
            if node in prune:
                self.remove(node)

        # remove pruned edges
        for node in self:
            for parent in list(node.parents):
                if len(self.get_ancestral(node, parent=parent)) == 0:
                    parent.children.remove(node)
                    node.parents.remove(parent)

        # remove single children
        if remove_single:
            remove_single_lineage(self)
            
        # set root
        # TODO: may need to actually use self.roots
        for node in self:
            if len(node.parents) == 0:
                self.root = node
                break
            

        

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

    arg = ARG()

    # make leaves
    lineages  = set((arg.add(CoalNode(arg.new_name())), 1)
                     for i in xrange(k))

    # process events
    for t, event in izip(times, events):
        if event == "coal":
            node = arg.add(CoalNode(arg.new_name(), age=t, event=event))
            a, b = random.sample(lineages, 2)
            lineages.remove(a)
            lineages.remove(b)
            node.children = [a[0], b[0]]
            a[0].parents.append(node)
            b[0].parents.append(node)
            lineages.add((node, 1))
            
        elif event == "recomb":
            node = arg.add(CoalNode(arg.new_name(), age=t, event=event))
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


def get_recomb_pos(arg):
    """
    Returns a sorted list of an ARG's recombination positions
    """
    rpos = [node.pos for node in
            arg if node.event == "recomb"]
    rpos.sort()
    return rpos


def iter_recomb_blocks(arg, start=None, end=None):
    """
    Iterates over the recombination blocks of an ARG
    """

    if start is None:
        start = arg.start
    if end is None:
        end = arg.end

    a = start
    b = start
    for pos in get_recomb_pos(arg):
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
    
    for a,b in iter_recomb_blocks(arg, start, end):
        yield arg.get_marginal_tree((a+b) / 2.0)
    

    
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


def remove_single_lineage(arg):
    """
    Remove unnecessary nodes with single parent and single child
    """
    for node in list(arg):
        if len(node.children) == 1 and len(node.parents) == 1:
            child = node.children[0]
            parent = node.parents[0]
            arg.remove(node)
            child.parents.append(parent)
            parent.children.append(child)

#=============================================================================
# mutations

def sample_mutations(arg, r):
    """

    r -- recombination rate (recomb/locus/gen)
    """

    mutations = []

    locsize = arg.end - arg.start

    for node in arg:
        for parent in node.parents:
            for region in arg.get_ancestral(node, parent=parent):
                frac = (region[1] - region[0]) / locsize
                dist = parent.age - node.age
                t = parent.age
                while True:
                    t -= random.expovariate(r * frac)
                    if t < node.age:
                        break
                    pos = random.uniform(region[0], region[1])
                    mutations.append((node, parent, pos, t))

    return mutations
    




#=============================================================================
# helper functions
        
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


