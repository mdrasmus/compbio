"""

   Code for the DLCoal model
   (duplications, losses, and coalescence)

"""

from __future__ import division

# python libs
import os
import random
from itertools import chain, izip
from math import *

# rasmus libs
from rasmus import util, treelib

# compbio libs
from . import birthdeath
from . import coal
from . import phylo

# spidir libs
import spidir
import spidir.topology_prior


#reload(phylo)


#=============================================================================
# reconciliation

def dlcoal_recon(tree, stree, gene2species,
                 n, duprate, lossrate,
                 pretime=None, premean=None,
                 nsearch=1000,
                 maxdoom=20, nsamples=100,
                 search=phylo.TreeSearchNni):
    """
    Perform reconciliation using the DLCoal model

    Returns (maxp, maxrecon) where 'maxp' is the probability of the
    MAP reconciliation 'maxrecon' which further defined as

    maxrecon = {'coal_recon': coal_recon,
                'locus_tree': locus_tree,
                'locus_recon': locus_recon,
                'locus_events': locus_events,
                'daughters': daughters}
    
    """

    reconer = DLCoalRecon(tree, stree, gene2species,
                          n, duprate, lossrate,
                          pretime=pretime, premean=premean,
                          maxdoom=maxdoom, nsamples=nsamples)
    reconer.set_proposer(DLCoalReconProposer(search=search))
    return reconer.recon(nsearch)



class DLCoalRecon (object):

    def __init__(self, tree, stree, gene2species,
                 n, duprate, lossrate,
                 pretime=None, premean=None,
                 maxdoom=20, nsamples=100,
                 name_internal="n"):

        # init coal tree
        self.coal_tree = tree
        self.stree = stree
        self.gene2species = gene2species
        self.n = n
        self.duprate = duprate
        self.lossrate = lossrate
        self.pretime = pretime
        self.premean = premean
        self.maxdoom = maxdoom
        self.nsamples = nsamples
        self.name_internal = name_internal

        self.proposer = DLCoalReconProposer()


    def set_proposer(self, proposer):
        """Set the proposal algorithm"""
        self.proposer = proposer
        

    def recon(self, nsearch=1000):
        """Perform reconciliation"""
        
        self.init_search()
        for i in xrange(nsearch):
            print "search", i
            proposal = self.proposer.next_proposal()
            p = self.eval_proposal(proposal)
            self.eval_search(p, proposal)
        
        # rename locus tree nodes
        rename_nodes(self.maxrecon["locus_tree"], self.name_internal)
        
        return self.maxp, self.maxrecon


    def init_search(self):
        """Initialize new search"""

        # init locus tree as congruent to coal tree
        # equivalent to assuming no ILS
        self.proposer.set_reconer(self)
        self.proposer.set_locus_tree(self.coal_tree.copy())

        self.maxp = - util.INF
        self.maxrecon = None


    def next_proposal(self):
        """Returns next proposal"""
        self.proposal.next_proposal()


    def eval_proposal(self, proposal):
        """Compute probability of proposal"""

        # compute recon probability
        phylo.add_implied_spec_nodes(proposal["locus_tree"], self.stree,
                                     proposal["locus_recon"],
                                     proposal["locus_events"])
        p = prob_dlcoal_recon_topology(self.coal_tree,
                                       proposal["coal_recon"],
                                       proposal["locus_tree"],
                                       proposal["locus_recon"],
                                       proposal["locus_events"],
                                       proposal["daughters"],
                                       self.stree, self.n,
                                       self.duprate, self.lossrate,
                                       self.pretime, self.premean,
                                       maxdoom=self.maxdoom,
                                       nsamples=self.nsamples,
                                       add_spec=False)
        treelib.remove_single_children(proposal["locus_tree"])
        phylo.subset_recon(proposal["locus_tree"], proposal["locus_recon"])

        return p


    def eval_search(self, p, proposal):
        """Evaluate a proposal for search"""

        if p > self.maxp:
            self.maxp = p
            self.maxrecon = proposal

            # search with a new copy
            self.proposer.accept()
        else:
            self.proposer.reject()



class DLCoalReconProposer (object):

    def __init__(self, search=phylo.TreeSearchNni):
        self.reconer = None
        self.locus_search = search(None)
        

    def set_reconer(self, reconer):
        self.reconer = reconer

    def set_locus_tree(self, locus_tree):
        self.locus_search.set_tree(locus_tree)

    def next_proposal(self):        
        self.locus_search.propose()
        
        # TODO: propose other reconciliations beside LCA
        locus_tree = self.locus_search.get_tree().copy()
        phylo.recon_root(locus_tree, self.reconer.stree,
                         self.reconer.gene2species,
                         newCopy=False)
        locus_recon = phylo.reconcile(locus_tree, self.reconer.stree,
                                      self.reconer.gene2species)
        locus_events = phylo.label_events(locus_tree, locus_recon)

        # propose daughters (TODO)
        daughters = set()

        # propose coal recon (TODO: propose others beside LCA)
        coal_recon = phylo.reconcile(self.reconer.coal_tree,
                                     locus_tree, lambda x: x)

        recon = {"coal_recon": coal_recon,
                 "locus_tree": locus_tree,
                 "locus_recon": locus_recon,
                 "locus_events": locus_events,
                 "daughters": daughters}
        return recon

    #def get_proposal(self):
    #    


    def accept(self):
        self.locus_search.set_tree(self.locus_search.get_tree().copy())
    

    def reject(self):
        self.locus_search.revert()



#=============================================================================
# probability functions for DLCoal model

def prob_dlcoal_recon_topology(coal_tree, coal_recon,
                               locus_tree, locus_recon, locus_events,
                               daughters,
                               stree, n, duprate, lossrate,
                               pretime=None, premean=None,
                               maxdoom=20, nsamples=100,
                               add_spec=True):
    """
    Probability of a reconcile gene tree in the DLCoal model.

    coal_tree    -- coalescent tree
    coal_recon   -- reconciliation of coalescent tree to locus tree
    locus_tree   -- locus tree (has dup-loss)
    locus_recon  -- reconciliation of locus tree to species tree
    locus_events -- events dict for locus tree
    stree        -- species tree
    n            -- population sizes in species tree
    duprate      -- duplication rate
    lossrate     -- loss rate

    You must also specify one of the following
    pretime      -- starting time before species tree
    premean      -- mean starting time before species tree

    Note: locus tree must have implied speciation nodes present
    """

    dups = phylo.count_dup(locus_tree, locus_events)

    # ensure implicit speciations are present
    if add_spec:
        phylo.add_implied_spec_nodes(locus_tree, stree,
                                     locus_recon, locus_events)
    
    # init popsizes for locus tree
    stree_popsizes = coal.init_popsizes(stree, n)
    popsizes = {}
    for node in locus_tree:
        popsizes[node.name] = stree_popsizes[locus_recon[node].name]


    # duploss probability

    #util.tic("top")
    dl_prob = spidir.calc_birth_death_prior(locus_tree, stree, locus_recon,
                                            duprate, lossrate,
                                            maxdoom=maxdoom)
    #util.toc()
    
    # daughters probability
    d_prob = dups * log(.5)


    # integrate over duplication times using sampling
    prob = 0.0
    #util.tic("int")
    for i in xrange(nsamples):
        # sample duplication times

        locus_times = spidir.topology_prior.sample_dup_times(
            locus_tree, stree, locus_recon, duprate, lossrate, pretime,
            premean,
            events=locus_events)
        assert len(locus_times) == len(locus_tree.nodes), (
            len(locus_times), len(locus_tree.nodes))
        birthdeath.set_dists_from_timestamps(locus_tree, locus_times)

        # coal topology probability
        coal_prob = prob_coal_recon_topology(coal_tree, coal_recon,
                                             locus_tree, popsizes, daughters)
        
        prob += exp(coal_prob)
        print coal_prob
    #util.toc()

    return dl_prob + d_prob + util.safelog(prob / nsamples, -util.INF)




def prob_coal_recon_topology(tree, recon, locus_tree, n, daughters):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a locus_tree 'locus_tree',
    population sizes 'n', and daughters set 'daughters'
    """

    # init population sizes
    popsizes = coal.init_popsizes(locus_tree, n)

    # log probability
    lnp = 0.0

    nodes = set(tree.postorder())

    # init reverse reconciliation
    rev_recon = {}
    for node, snode in recon.iteritems():
        if node not in nodes:
            raise Exception("node '%s' not in tree" % node.name)
        rev_recon.setdefault(snode, []).append(node)

    # init lineage counts
    lineages = {}
    for snode in locus_tree:
        if snode.is_leaf():
            lineages[snode] = len([x for x in rev_recon[snode]
                                   if x.is_leaf()])
        else:
            lineages[snode] = 0

    # iterate through species tree branches
    for snode in locus_tree.postorder():
        if snode.parent:
            # non root branch
            u = lineages[snode]

            # subtract number of coals in branch
            v = u - len([x for x in rev_recon.get(snode, [])
                         if not x.is_leaf()])            
            lineages[snode.parent] += v

            if snode not in daughters:
                try:
                    lnp += util.safelog(
                        coal.prob_coal_counts(u, v, snode.dist,
                                              popsizes[snode.name]),
                        -util.INF)
                except:
                    print u, v, snode.dist, popsizes[snode.name]
                    raise
            else:
                assert v == 1
            lnp -= log(coal.num_labeled_histories(u, v))
        else:
            # normal coalesent
            u = lineages[snode]
            lnp -= log(coal.num_labeled_histories(u, 1))

    
    # correct for topologies H(T)
    # find connected subtrees that are in the same species branch
    subtrees = []
    subtree_root = {}
    for node in tree.preorder():
        if node.parent and recon[node] == recon[node.parent]:
            subtree_root[node] = subtree_root[node.parent]
        else:
            subtrees.append(node)
            subtree_root[node] = node

    # find leaves through recursion
    def walk(node, subtree, leaves):
        if node.is_leaf():
            leaves.append(node)
        elif (subtree_root[node.children[0]] != subtree and
              subtree_root[node.children[1]] != subtree):
            leaves.append(node)
        else:
            for child in node.children:
                walk(child, subtree, leaves)

    # apply correction for each subtree
    for subtree in subtrees:
        leaves = []
        for child in subtree.children:
            walk(subtree, subtree, leaves)
        if len(leaves) > 2:
            lnp += log(birthdeath.num_topology_histories(subtree, leaves))

    return lnp



#=============================================================================
# sampling from the DLCoal model

def sample_dlcoal(stree, n, duprate, lossrate, namefunc=lambda x: x,
                  remove_single=True, name_internal="n",
                  minsize=0):
    """Sample a gene tree from the DLCoal model"""

    # generate the locus tree
    while True:
        locus_tree, locus_recon, locus_events = \
                    birthdeath.sample_birth_death_gene_tree(
            stree, duprate, lossrate)
        if len(locus_tree.leaves()) >= minsize:
            break

    if len(locus_tree.nodes) <= 1:
        # total extinction
        coal_tree = treelib.Tree()
        coal_tree.make_root()
        coal_recon = {coal_tree.root: locus_tree.root}
        daughters = set()
    else:
        # simulate coalescence
        
        # choose daughter duplications
        daughters = set()
        for node in locus_tree:
            if locus_events[node] == "dup":
                daughters.add(node.children[random.randint(0, 1)])

        coal_tree, coal_recon = sample_multicoal_tree(locus_tree, n,
                                                      daughters=daughters,
                                                      namefunc=namefunc)

        # clean up coal tree
        if remove_single:
            treelib.remove_single_children(coal_tree)
            phylo.subset_recon(coal_tree, coal_recon)

    if name_internal:
        rename_nodes(coal_tree, name_internal)
        rename_nodes(locus_tree, name_internal)


    # store extra information
    extra = {"locus_tree": locus_tree,
             "locus_recon": locus_recon,
             "locus_events": locus_events,
             "coal_tree": coal_tree,
             "coal_recon": coal_recon,
             "daughters": daughters}

    return coal_tree, extra


def rename_nodes(tree, prefix="n"):
    """Rename nodes that all names are strings"""
    for node in list(tree.postorder()):
        if isinstance(node.name, int):
            name2 = prefix + str(node.name)
            while name2 in tree.nodes:
                name2 = prefix + str(tree.new_name())
            tree.rename(node.name, name2)




def sample_multicoal_tree(stree, n, leaf_counts=None,
                          daughters=set(),
                          namefunc=None):
    """
    Returns a gene tree from a multi-species coalescence process
    n -- population size (int or dict)
         If n is a dict it must map from species name to population size
    """

    # initialize vector for how many genes per extant species
    if leaf_counts is None:
        leaf_counts = dict((l, 1) for l in stree.leaf_names())

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l, 1) for l in stree.leaf_names())
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    # initialize population sizes
    popsizes = coal.init_popsizes(stree, n)

    # init gene counts
    counts = dict((n.name, 0) for n in stree)
    counts.update(leaf_counts)

    # init reconciliation, events
    recon = {}

    subtrees = {}

    # loop through species tree
    for snode in stree.postorder():        
        # simulate population for one branch
        k = counts[snode.name]
        if snode in daughters:
            # daughter branch, use bounded coalescent
            subtree = coal.sample_coal_tree_bounded(
                k, popsizes[snode.name], snode.dist, capped=True)
            lineages = set(subtree.root)
        elif snode.parent:            
            # non basal branch
            subtree, lineages = coal.sample_coal_tree_fixed(
                k, popsizes[snode.name], snode.dist, capped=True)
        else:
            # basal branch
            subtree = coal.sample_coal_tree(k, popsizes[snode.name])
            lineages = set(subtree.root)
        subtrees[snode] = (subtree, lineages)
        if snode.parent:
            counts[snode.parent.name] += len(lineages)
        for node in subtree:
            recon[node] = snode


    # stitch subtrees together
    tree = treelib.Tree()

    # add all nodes to total tree
    for snode, (subtree, lineages) in subtrees.iteritems():
        tree.merge_names(subtree)
        if snode.parent:
            tree.remove(subtree.root)        
            del recon[subtree.root]
    
    for snode in stree:
        if not snode.is_leaf():
            subtree, lineages = subtrees[snode]

            # get lineages from child subtrees
            lineages2 = chain(*[subtrees[child][1]
                                for child in snode.children])

            # ensure leaves are randomly attached
            leaves = subtree.leaves()
            random.shuffle(leaves)

            # stitch leaves of the subtree to children subtree lineages
            for leaf, lineage in izip(leaves, lineages2):
                tree.add_child(leaf, lineage)


    # set root
    tree.root = subtrees[stree.root][0].root

    # name leaves
    for leaf in tree.leaves():
        tree.rename(leaf.name, namefunc(recon[leaf].name))
        
    return tree, recon



#=============================================================================
# Input/Output


def dlcoal_sims(outdir, nsims, stree, n, duprate, lossrate,
                start=0,
                **options):
    
    for i in xrange(start, nsims):
        outfile = phylo.phylofile(outdir, str(i), "")
        util.makedirs(os.path.dirname(outfile))
        print "simulating", outfile

        # sample a new tree from DLCoal model
        coal_tree, ex = sample_dlcoal(stree, n, duprate, lossrate, **options)
        write_dlcoal_recon(outfile, coal_tree, ex)










def write_dlcoal_recon(filename, coal_tree, extra,
                       exts={"coal_tree": ".coal.tree",
                             "coal_recon": ".coal.recon",
                             "locus_tree": ".locus.tree",
                             "locus_recon": ".locus.recon",
                             "daughters": ".daughters"
                             },
                       filenames={}):
    """Writes a reconciled gene tree to files"""

    # coal
    coal_tree.write(filenames.get("coal_tree", filename + exts["coal_tree"]),
                    rootData=True)
    phylo.write_recon_events(
        filenames.get("coal_recon", filename + exts["coal_recon"]),
        extra["coal_recon"], noevent="none")

    # locus
    extra["locus_tree"].write(
        filenames.get("locus_tree", filename + exts["locus_tree"]),
        rootData=True)
    phylo.write_recon_events(
        filenames.get("locus_recon", filename + exts["locus_recon"]),
        extra["locus_recon"], extra["locus_events"])

    util.write_list(
        filenames.get("daughters", filename + exts["daughters"]),
        [x.name for x in extra["daughters"]])



def read_dlcoal_recon(filename, stree,
                      exts={"coal_tree": ".coal.tree",
                            "coal_recon": ".coal.recon",
                            "locus_tree": ".locus.tree",
                            "locus_recon": ".locus.recon",
                            "daughters": ".daughters"
                            },
                      filenames={}):
    """Reads a reconciled gene tree from files"""

    extra = {}

    # trees
    coal_tree = treelib.read_tree(
        filenames.get("coal_tree", filename + exts["coal_tree"]))
    extra["locus_tree"] = treelib.read_tree(
        filenames.get("locus_tree", filename + exts["locus_tree"]))

    # recons
    extra["coal_recon"], junk = phylo.read_recon_events(
        filenames.get("coal_recon", filename + exts["coal_recon"]),
        coal_tree, extra["locus_tree"])
    extra["locus_recon"], extra["locus_events"] = phylo.read_recon(
        filenames.get("locus_recon", filename + exts["locus_recon"]),
        extra["locus_tree"], stree)


    extra["daughters"] = set(
        extra["locus_tree"].nodes[x] for x in util.read_strings(
        filenames.get("daughters", filename + exts["daughters"])))

    return coal_tree, extra
    


#=============================================================================
# old code

def write_dlcoal_recon2(filename, coal_tree, extra,
                       exts={"coal_tree": ".coal.tree",
                             "coal_recon": ".coal.recon",
                             "locus_tree": ".locus.tree",
                             "locus_recon": ".locus.recon",
                             "locus_events": ".locus.events",
                             "daughters": ".daughters"
                             },
                       filenames={}):
    """Writes a reconciled gene tree to files"""

    # coal
    coal_tree.write(filenames.get("coal_tree", filename + exts["coal_tree"]),
                    rootData=True)
    phylo.write_recon(
        filenames.get("coal_recon", filename + exts["coal_recon"]),
        extra["coal_recon"])

    # locus
    extra["locus_tree"].write(
        filenames.get("locus_tree", filename + exts["locus_tree"]),
        rootData=True)
    phylo.write_recon(
        filenames.get("locus_recon", filename + exts["locus_recon"]),
        extra["locus_recon"])
    phylo.write_events(
        filenames.get("locus_events", filename + exts["locus_events"]),
        extra["locus_events"])

    util.write_list(
        filenames.get("daughters", filename + exts["daughters"]),
        [x.name for x in extra["daughters"]])



def read_dlcoal_recon2(filename, stree,
                      exts={"coal_tree": ".coal.tree",
                            "coal_recon": ".coal.recon",
                            "locus_tree": ".locus.tree",
                            "locus_recon": ".locus.recon",
                            "locus_events": ".locus.events",
                            "daughters": ".daughters"
                            },
                      filenames={}):
    """Reads a reconciled gene tree from files"""

    extra = {}

    # trees
    coal_tree = treelib.read_tree(
        filenames.get("coal_tree", filename + exts["coal_tree"]))
    extra["locus_tree"] = treelib.read_tree(
        filenames.get("locus_tree", filename + exts["locus_tree"]))

    # recons
    extra["coal_recon"] = phylo.read_recon(
        filenames.get("coal_recon", filename + exts["coal_recon"]),
        coal_tree, extra["locus_tree"])
    extra["locus_recon"] = phylo.read_recon(
        filenames.get("locus_recon", filename + exts["locus_recon"]),
        extra["locus_tree"], stree)
    extra["locus_events"] = phylo.read_events(
        filenames.get("locus_events", filename + exts["locus_events"]),
        extra["locus_tree"])


    extra["daughters"] = set(
        extra["locus_tree"].nodes[x] for x in util.read_strings(
        filenames.get("daughters", filename + exts["daughters"])))

    return coal_tree, extra
    














#=============================================================================
# OLD


def dlcoal_recon_old(tree, stree, gene2species,
                 n, duprate, lossrate,
                 pretime=None, premean=None,
                 nsearch=1000,
                 maxdoom=20, nsamples=100,
                 search=phylo.TreeSearchNni):
    """
    Perform reconciliation using the DLCoal model

    Returns (maxp, maxrecon) where 'maxp' is the probability of the
    MAP reconciliation 'maxrecon' which further defined as

    maxrecon = {'coal_recon': coal_recon,
                'locus_tree': locus_tree,
                'locus_recon': locus_recon,
                'locus_events': locus_events,
                'daughters': daughters}
    
    """

    # init coal tree
    coal_tree = tree

    # init locus tree as congruent to coal tree
    # equivalent to assuming no ILS
    locus_tree = coal_tree.copy()

    maxp = - util.INF
    maxrecon = None

    # init search
    locus_search = search(locus_tree)

    for i in xrange(nsearch):       
        # TODO: propose other reconciliations beside LCA
        locus_tree2 = locus_tree.copy()
        phylo.recon_root(locus_tree2, stree, gene2species, newCopy=False)
        locus_recon = phylo.reconcile(locus_tree2, stree, gene2species)
        locus_events = phylo.label_events(locus_tree2, locus_recon)

        # propose daughters (TODO)
        daughters = set()

        # propose coal recon (TODO: propose others beside LCA)
        coal_recon = phylo.reconcile(coal_tree, locus_tree2, lambda x: x)

        # compute recon probability
        phylo.add_implied_spec_nodes(locus_tree2, stree,
                                     locus_recon, locus_events)
        p = prob_dlcoal_recon_topology(coal_tree, coal_recon,
                                       locus_tree2, locus_recon, locus_events,
                                       daughters,
                                       stree, n, duprate, lossrate,
                                       pretime, premean,
                                       maxdoom=maxdoom, nsamples=nsamples,
                                       add_spec=False)
        treelib.remove_single_children(locus_tree2)

        if p > maxp:
            maxp = p
            maxrecon = {"coal_recon": coal_recon,
                        "locus_tree": locus_tree2,
                        "locus_recon": locus_recon,
                        "locus_events": locus_events,
                        "daughters": daughters}
            locus_tree = locus_tree2.copy()
            locus_search.set_tree(locus_tree)
        else:
            locus_search.revert()

        # perform local rearrangement to locus tree
        locus_search.propose()




    return maxp, maxrecon




