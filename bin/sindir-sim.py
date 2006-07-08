#!/usr/bin/env python

# python libs
import sys, time, os, random, math

# rasmus libs
from rasmus import env, treelib
from rasmus import fasta, util
from rasmus import sindirlib, stats

# scipy libs
from scipy import matrixmultiply as mm
from scipy import transpose, array, dot
from scipy.linalg import expm, logm, eig



options = [
"""\
    Simulates sequence evolution using SINDIR's model. 
""",

    ["s:", "stree=", "stree", "<species tree>",
        {"single": True}],
    ["S:", "smap=", "smap", "<gene2species map>"],
    ["p:", "params=", "params", "<params>",
        {"single": True}],
    ["r:", "baserate=", "baserate", "<baserate>",
        {"single": True,
         "default": 1.0,
         "parser": float}],        
    
    "Sequence generation",
    ["n:", "ntrees=", "ntrees", "<number of trees to produce>",
        {"single": True,
         "default": 1,
         "parser": int}],    
    ["l:", "seqlen=", "seqlen", "<length of starting sequence>",
        {"single": True,
         "default": 0,
         "parser": int}],
    ["q:", "seqfile=", "seqfile", "<FASTA file of starting sequences>",
        {"default": []}],
    ["G:", "gc=", "gc", "<GC content>",
        {"single": True,
         "default": .5,
         "parser": float}],
    
    "Duplication/Loss configuration",
    ["D:", "dupprob=", "dupprob", "<duplication probability>",
        {"single": True,
         "default": 0,
         "parser": float}],
    ["P:", "predup=", "predup", "<number of duplication before speciations>",
        {"single": True,
         "default": 0,
         "parser": int}],
    ["L:", "lossprob=", "lossprob", "<loss probability>",
        {"single": True,
         "default": 0,
         "parser": float}],
    
    "Kimura model parameters",
    ["A:", "alpha=", "alpha", "<kimura alpha>",
        {"single": True,
         "default": .25,
         "parser": float}],
    ["B:", "beta=", "beta", "<kimura beta>",
        {"single": True,
         "default": .25,
         "parser": float}],
    
    "Codon simulation",
    ["", "csm=", "csm", "<Codon Substitution Matrix>",
        {"single": True}],
    ["", "csmbranchlen=", "csmbranchlen", "<branch length for species in CSM>",
        {"single": True,
         "parser": float}],
    
    "Output extensions",
    ["t:", "outtree=", "outtree", "<output tree filename extension>",
        {"single": True,
         "default": ".tree"}],
    ["a:", "outalign=", "outalign", "<output align filename extension>",
        {"single": True,
         "default": ".align"}],
    ["d:", "outdist=", "outdist", "<output distance filename extension>",
        {"single": True,
         "default": ".dist"}],
    
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
    if conf["dupprob"] >= 1:
        conf["dupnum"] = int(conf["dupprob"])
    else:
        conf["dupnum"] = 0
        
    if conf["lossprob"] >= 1:
        conf["lossnum"] = int(conf["lossprob"])
    else:
        conf["lossnum"] = 0
    
    
    
    # adjustment for baserate (empirically determined)
    #conf["baserate"] = conf["baserate"] / .9
    
    
    stree = treelib.readTree(env.findFile(conf["stree"]))
    params = sindirlib.readParams(conf["params"])
    
    # read species2gene
    for f in conf["smap"]:
        smap = util.readDelim(env.findFile(f))
    species2gene = {}
    
    for prefix, genome in smap:
        species2gene[genome] = prefix.replace("*", "")
    
    # read starting sequences
    rootseqs = []
    for f in conf["seqfile"]:
        rootseqs.extend(fasta.readFasta(f).values())
    rootseqs = [s.replace("-", "") for s in rootseqs]
    
    
    # load CSM
    if "csm" in conf:
        labels, conf["csm"] = loadCsm(conf["csm"])
        
        # remove AA from labels
        labels = [x[:x.index("/")] for x in labels]
        conf["csmlabels"] = labels
        
        # make csm lookup
        conf["csmlookup"] = {}
        for i, label in enumerate(labels):
            conf["csmlookup"][label] = i
        
        # find rate matrix
        conf["csmrate"] = logm(conf["csm"]) / conf["csmbranchlen"]
        
    
    
    # simulate
    for i in range(conf["ntrees"]):
        util.log("simulating", i)
        
        # create initial sequence
        if conf["seqlen"] > 0:
            seq = []
            gc = conf["gc"]
            for k in range(conf["seqlen"]):
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
            rootseq = rootseqs[i % len(rootseqs)]
        
        
        while True:
            tree = simTree(conf, rootseq, stree, params, species2gene)
            treelib.removeSingleChildren(tree)
            
            # make another tree if this one is too simple
            if len(tree.leaves()) >= 4:
                break
        
        # write output
        tree.write(str(i) + conf["outtree"])

        seqs = fasta.FastaDict()
        out = file(str(i) + conf["outalign"], "w")
        for leaf in tree.leaves():
            print >>out, ">%s" % leaf.name
            print >>out, leaf.data["seq"]
            seqs[leaf.name] = leaf.data["seq"]
        out.close()
        
        #phylip.dnadist(seqs, str(i) + conf["outdist"], verbose=False)


def real(x):
    if isinstance(x, float):
        return x
    else:
        return x.real


def simTree(conf, rootseq, stree, params, species2gene):
    tree = treelib.Tree()
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
        if random.random() < conf["lossprob"] < 1.0:
            return    
    
        # decide if duplication occurs
        if random.random() < conf["dupprob"] < 1.0 or \
           dups[snode] > 0:
            dup = True
            dups[snode] -= 1
            
            # choose dup point
            endpoint = midpoint + (1 - midpoint) * random.random()
        else:
            dup = False
            endpoint = 1.0
        
        # choose branch length
        branchlen = simBranchLen(conf, conf["baserate"], 
                                 params[snode.name], midpoint, endpoint)
        
        
        # mutate sequence
        if "csm" in conf:
            seq2 = simSeqCsm(conf, gparent.data["seq"], branchlen)
        else:
            seq2 = simSeq(conf, gparent.data["seq"], branchlen)
        
        # setup new gene node
        if snode.isLeaf() and not dup:
            prefix = species2gene[snode.name]
            gchild = treelib.TreeNode(prefix + "_sim_" + str(tree.newName()))
        else:
            gchild = treelib.TreeNode(tree.newName())
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
    
    if conf["predup"] == 0:
        for schild in stree.root.children:
            walk(schild, tree.root, 0)
    else:
        dups[stree.root] = conf["predup"]
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
    treelib.removeSingleChildren(tree)
    
    assert int not in [type(x.name) for x in tree.leaves()]
    
    return tree



def simBranchLen(conf, baserate, param, midpoint, endpoint):
    blen = 0
    
    frac = endpoint - midpoint
    
    while blen <= 0:
        blen = baserate * random.normalvariate(frac * param[0], frac * param[1])
    
    return blen


def simSeq(conf, seq, branchlen):
    alpha = conf["alpha"]
    beta = conf["beta"]
    
    seq2 = []
    for base in seq:
        seq2.append(kimura(base, alpha, beta, branchlen))
    
    return "".join(seq2)


def simSeqCsm(conf, seq, branchlen):
    seq2 = []
    
    csm2 = util.map2(real, expm(branchlen * conf["csmrate"]))
    lookup = conf["csmlookup"]
    labels = conf["csmlabels"]
    
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        codon2 = labels[stats.sample(csm2[lookup[codon]])]
        seq2.append(codon2)
    
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
    

def readCsm(filename):
    data = util.readDelim(filename)

    labels = data[0]
    mat = []

    for i, row in enumerate(data[1:]):
        assert row[0] == labels[i]
        mat.append(map(lambda x: int(float(x)), row[1:]))

    return labels, mat


def normMatrix(mat):
    rowsum = map(sum, mat)
    
    return [[x/float(tot) for x in row]
            for row, tot in zip(mat, rowsum)]

def loadCsm(filename):
    labels, csm_counts = readCsm(filename)
    return labels, normMatrix(csm_counts)


main(conf)
