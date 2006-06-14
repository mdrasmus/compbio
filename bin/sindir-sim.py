#!/usr/bin/env python

import sys, time, os, random, math
from rasmus import synteny, synphylweb, bionj, ensembl, phyloutil, env
from rasmus import fasta, util, genomeio, genomeutil, phylip
from rasmus import muscle, phylip, algorithms, clustalw, alignlib
from rasmus import sindirlib


options = [
    ["p:", "params=", "params", "<params>"],
    ["s:", "stree=", "stree", "<species tree>"],
    ["S:", "smap=", "smap", "<gene2species map>"],
    ["q:", "seq=", "seq", "<seq>"],    
    ["D:", "dupprob=", "dupprob", "<duplication probability>"],
    ["P:", "predup=", "predup", "<number of duplication before speciations>"],
    ["L:", "lossprob=", "lossprob", "<loss probability>"],
    ["A:", "alpha=", "alpha", "<kimura alpha>"],
    ["B:", "beta=", "beta", "<kimura beta>"],
    ["G:", "gc=", "gc", "<GC content>"],
    ["r:", "baserate=", "baserate", "<baserate>"],
    ["t:", "outtree=", "outtree", "<output tree filename extension>"],
    ["a:", "outalign=", "outalign", "<output align filename extension>"],
    ["d:", "outdist=", "outdist", "<output distance filename extension>"],
    ["n:", "ntrees=", "ntrees", "<number of trees to produce>"]
]


conf = util.parseOptions(sys.argv, options, quit=True)


BASE2INT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}

INT2BASE = ["A", "C", "G", "T"]

KIMURA_MATRIX = [
    ['r', 's', 'u', 's'],
    ['s', 'r', 's', 'u'],
    ['u', 's', 'r', 's'],
    ['s', 'u', 's', 'r']
]


def main(conf):
    env.addEnvPaths("DATAPATH")
    
    # parse options
    conf["dupnum"] = 0
    if "dupprob" not in conf:
        conf["dupprob"] = [0]
    else:
        conf["dupprob"][-1] = float(conf["dupprob"][-1])
        if conf["dupprob"][-1] >= 1:
            conf["dupnum"] = int(conf["dupprob"][-1])
        
    conf["lossnum"] = 0
    if "lossprob" not in conf:
        conf["lossprob"] = [0]
    else:
        conf["lossprob"][-1] = float(conf["lossprob"][-1])
        if conf["lossprob"][-1] >= 1:
            conf["lossnum"] = int(conf["lossprob"][-1])
    
    if "predup" not in conf:
        conf["predup"] = [0]
    else:
        conf["predup"][-1] = int(conf["predup"][-1])
    
    if "alpha" not in conf:
        conf["alpha"] = [.25]
    else:
        conf["alpha"][-1] = float(conf["alpha"][-1])
    
    if "beta" not in conf:
        conf["beta"] = [.25]
    else:
        conf["beta"][-1] = float(conf["beta"][-1])
    
    if "gc" not in conf:
        conf["gc"] = [.42]
    else:
        conf["gc"][-1] = float(conf["gc"][-1])
    
    if "baserate" not in conf:
        conf["baserate"] = [1.0]
    else:
        conf["baserate"][-1] = float(conf["baserate"][-1])
    
    if "ntrees" not in conf:
        conf["ntrees"] = [1]
    else:
        conf["ntrees"][-1] = int(conf["ntrees"][-1])
    
    if conf["seq"][-1].isdigit():
        seq = []
        gc = conf["gc"][-1]
        for i in range(int(conf["seq"][-1])):
            j = random.random() 
            if j < gc:
                if j < gc / 2:
                    seq.append("G")
                else:
                    seq.append("C")
            else:
                if j < gc + (1 - gc) / 2:
                    seq.append("A")
                else:
                    seq.append("T")
        rootseq = "".join(seq)
    else:
        rootseq = conf["seq"][-1]
    
    
    # adjustment for baserate (empirically determined)
    conf["baserate"][-1] = conf["baserate"][-1] / .9
    
    
    stree = algorithms.readTree(env.findFile(conf["stree"][-1]))
    params = sindirlib.readParams(conf["params"][-1])
    
    # read species2gene
    smap = util.readDelim(env.findFile(conf["smap"][-1]))
    species2gene = {}
    
    for prefix, genome in smap:
        species2gene[genome] = prefix
    
    
    # simulate
    for i in range(conf["ntrees"][-1]):
        util.log("simulating", i)
        
        while True:
            tree = simTree(conf, rootseq, stree, params, species2gene)
            algorithms.removeSingleChildren(tree)
            
            # make another tree if this one is too simple
            if len(tree.leaves()) >= 4:
                break
        
        # write output
        tree.write(str(i) + conf["outtree"][-1])

        seqs = fasta.FastaDict()
        out = file(str(i) + conf["outalign"][-1], "w")
        for leaf in tree.leaves():
            print >>out, ">%s" % leaf.name
            print >>out, leaf.data["seq"]
            seqs[leaf.name] = leaf.data["seq"]
        out.close()
        
        phylip.dnadist(seqs, str(i) + conf["outdist"][-1], verbose=False)



def simTree(conf, rootseq, stree, params, species2gene):
    tree = algorithms.Tree()
    tree.makeRoot()
    
    # hard-coded pre-dup lengths
    params[stree.root.name] = [.25, .125]
    
    
    tree.root.data["seq"] = rootseq
    
    # determine where losses and dups should be
    dups = util.Dict(1, 0)    
    for i in range(conf["dupnum"]):
        while True:
            x = int(random.random() * len(stree.nodes))
            if stree.nodes.values()[x] != stree.root: break
        dups[stree.nodes.values()[x]] += 1    
    
    
    
    def walk(snode, gparent, midpoint):       
        # decide if this branch is lost
        if random.random() < conf["lossprob"][-1] < 1.0:
            return    
    
        # decide if duplication occurs
        if random.random() < conf["dupprob"][-1] < 1.0 or \
           dups[snode] > 0:
            dup = True
            dups[snode] -= 1
            
            # choose dup point
            endpoint = midpoint + (1 - midpoint) * random.random()
        else:
            dup = False
            endpoint = 1.0
        
        # choose branch length
        branchlen = simBranchLen(conf, conf["baserate"][-1], 
                                 params[snode.name], midpoint, endpoint)
        
        
        # mutate sequence
        seq2 = simSeq(conf, gparent.data["seq"], branchlen)
        
        # setup new gene node
        if snode.isLeaf() and not dup:
            prefix = species2gene[snode.name]
            gchild = algorithms.TreeNode(prefix + "_sim_" + str(tree.newName()))
        else:
            gchild = algorithms.TreeNode(tree.newName())
        gchild.dist = branchlen
        gchild.data["seq"] = seq2
        tree.addChild(gparent, gchild)
        
        # recurse
        if dup:
            walk(snode, gchild, endpoint)
            walk(snode, gchild, endpoint)
        else:
            for schild in snode.children:
                walk(schild, gchild, 0)
        
        # if all children have been lost and we recon to a internal 
        # species node, then we should be removed as well
        if len(gchild.children) == 0 and not snode.isLeaf():
            tree.remove(gchild)
    
    if conf["predup"][-1] == 0:
        for schild in stree.root.children:
            walk(schild, tree.root, 0)
    else:
        dups[stree.root] = conf["predup"][-1]
        walk(stree.root, tree.root, 0)
    
    
    # apply numbered losses
    losses = util.Dict(1, 0)
    for i in range(conf["lossnum"]):
        while True:
            x = int(random.random() * len(tree.nodes))
            if tree.nodes.values()[x] != tree.root and \
               losses[tree.nodes.values()[x]] == 0: break
        losses[tree.nodes.values()[x]] = 1

    for key, val in losses.items():
        if val == 1:
            tree.removeTree(key)
    
    # find exposed ancestors
    while True:
        #print tree.leaveNames()
        exposed = filter(lambda x: type(x.name) == int, tree.leaves())
        if len(exposed) > 0:
            #print ">>", [x.name for x in exposed]
            for node in exposed:
                if node == tree.root:
                    # try to make another tree
                    return simTree(conf, rootseq, stree, params, species2gene)
                tree.removeTree(node)
        else:
            break
    
    # remove any redundant single children nodes
    algorithms.removeSingleChildren(tree)
    
    assert int not in [type(x.name) for x in tree.leaves()]
    
    return tree



def simBranchLen(conf, baserate, param, midpoint, endpoint):
    blen = 0
    
    frac = endpoint - midpoint
    
    while blen <= 0:
        blen = baserate * random.normalvariate(frac * param[0], frac * param[1])
    
    return blen


def simSeq(conf, seq, branchlen):
    alpha = conf["alpha"][-1]
    beta = conf["beta"][-1]

    seq2 = []
    for base in seq:
        seq2.append(kimura(base, alpha, beta, branchlen))
    
    return "".join(seq2)



def kimura(base, alpha, beta, time):
    probs = {
        's': .25 * (1 - math.e**(-4 * beta * time)),
        'u': .25 * (1 + math.e**(-4 * beta * time)
                      - 2*math.e**(-2*(alpha+beta)*time))
    }
    probs['r'] =  1 - 2*probs['s'] - probs['u']
    
    cdf = 0
    row = KIMURA_MATRIX[BASE2INT[base]]
    pick = random.random()
    
    for i in range(4):
        cdf += probs[row[i]]
        if cdf >= pick:
            return INT2BASE[i]
    
    assert False, "probabilities do not add to one"
        


main(conf)
