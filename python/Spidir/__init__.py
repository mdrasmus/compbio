#
# SPIDIR library
#
# note: SINDIR was a codename, it may still be present in the code.
#

# python libs
import math, StringIO, copy, random, sys, time


# rasmus libs
from rasmus import stats
from rasmus import tablelib
from rasmus import treelib
from rasmus import util

from rasmus.bio import bionj
from rasmus.bio import fasta
from rasmus.bio import phylo
from rasmus.bio import phylip

from rasmus.vis import treesvg




# SPIDIR libs
#from Spidir import Search
#from Spidir import Likelihood
#from Spidir.Debug import *

import Search
import Likelihood
from Debug import *

try:
    import pyspidir
except ImportError:
    pyspidir = None

# scipy libs
# (needed for numerical integration and least square error fitting)
try:
    import scipy
    import scipy.linalg
    import scipy.integrate
    import scipy.optimize
except:
    if pyspidir == None:
        raise Exception("Must install SciPy before using SPIDIR")



#-------------------------------------------------------------------------------
# SPIDIR input/output
#-------------------------------------------------------------------------------

def writeParams(filename, params):
    """Write SPIDIR model parameters to a file"""
    
    out = file(filename, "w")
    
    keys = util.sort(params.keys())
    
    for key in keys:
        values = params[key]
        print >>out, "%s\t%s" % (str(key), "\t".join(map(str,values)))


def readParams(filename):
    """Read SPIDIR model parameters to a file"""
    
    infile = file(filename)
    params = {}
    
    for line in infile:
        tokens = line.split("\t")
        key = tokens[0]
        values = tokens[1:]
        if key[0].isdigit():
            key = int(key)
        params[key] = map(float, values)
        
    return params


def readLabels(filename):
    """Read gene names from a file"""
    
    if filename.endswith(".fasta") or \
       filename.endswith(".fa") or \
       filename.endswith(".align"):
        labels = fasta.readFasta(filename).keys()
    else:
        labels = util.readStrings(filename)
    
    return labels


def writeTreeDistrib(out, lengths):
    out = util.openStream(out, "w")

    for node, lens in lengths.items():
        if len(lens) == 0 or max(lens) == min(lens):
            continue
        
        if isinstance(node, treelib.TreeNode):
            out.write(str(node.name))
        else:
            out.write(str(node))

        for length in lens:
            out.write("\t%f" % length)
        out.write("\n")


def treeDistrib2table(lengths, filenames=None):
    keys = []
    mat = []
    
    # build matrix of rates
    for node, lens in lengths.iteritems():
        if len(lens) == 0 or max(lens) == min(lens):
            continue
        
        if isinstance(node, treelib.TreeNode):
            keys.append(str(node.name))
        else:
            keys.append(str(node))
        
        mat.append(lens)

    # add filenames if they exist
    if filenames:
        mat.append(filenames)
        keys.append("filename")

    # transpose matrix
    mat = zip(* mat)
    
    rates = tablelib.Table(mat, headers=keys)
    
    return rates
    


def readTreeDistrib(filename):
    infile = util.openStream(filename)
    lengths = {}
    
    for line in infile:
        tokens = line.split("\t")
        name = tokens[0]
        
        if name.isdigit():
            name = int(name)
        
        lengths[name] = map(float, tokens[1:])
    
    return lengths


def outTreeFile(conf):
    return conf["out"] + ".tree"


def debugFile(conf):
    return conf["out"] + ".debug"


def drawParamTree(tree, params, *args, **kargs):
    tree = tree.copy()
    kargs.setdefault("legendScale", True)
    kargs.setdefault("xscale", 100)
    kargs.setdefault("yscale", 20)
    kargs.setdefault("minlen", 0)
    kargs.setdefault("maxlen", util.INF)
    
    if "labels" not in kargs:
        kargs["labels"] = {}
        
        for name in tree.nodes:
            if not tree.nodes[name].isLeaf():
                kargs["labels"][name] = str(name)
        
    
    # set branch lengths to means
    for name in tree.nodes:
        tree.nodes[name].dist = params[name][0]
    
    # draw basic tree
    tmargin = 10
    lmargin = 10
    canvas = treesvg.drawTree(tree, autoclose=False, 
                              tmargin=tmargin, lmargin=lmargin,
                              *args, **kargs)
    
    # draw variance
    coords = treelib.layoutTree(tree, kargs["xscale"], kargs["yscale"],
                                      kargs["minlen"], kargs["maxlen"])
    
    canvas.beginTransform(("translate", lmargin, tmargin))
    canvas.beginStyle("stroke-width: 3")
    for name, node in tree.nodes.iteritems():
        if node == tree.root:
            continue
        x, y = coords[node]
        
        if node.parent:
            parentx = coords[node.parent][0]
        else:
            parentx = 0
        
        varline = params[name][1] * (x - parentx) / params[name][0]
        
        canvas.line(x, y, max(parentx, x - varline), y, )
    canvas.endStyle()
    canvas.endTransform()
    
    canvas.endSvg()
        
    
    

#-------------------------------------------------------------------------------
# Branch length fitting
#-------------------------------------------------------------------------------


def setTreeDistances(conf, tree, distmat, genes):
    if isDebug(DEBUG_MED):
        util.tic("fit branch lengths")
    
    if pyspidir and "parsimony" in conf:
        # estimate branch lengths with parsimony
        parsimony_C(conf["aln"], tree)
        tree.data["error"] = sum(node.dist 
                                 for node in tree.nodes.itervalues())
    
    elif pyspidir and "mlhkydist" in conf:
        # estimate branch lengths with ML
        logl = mlhkydist_C(conf["aln"], tree, conf["bgfreq"], conf["tsvratio"], 
                           3*len(tree.nodes))
        tree.data["distlogl"] = logl
        tree.data["error"] = 0.0
    else:
        # perform LSE
        lse = phylo.leastSquareError(tree, distmat, genes)

        # catch unusual case that may occur in greedy search
        if sum(x.dist for x in tree.nodes.values()) == 0:
            for node in tree.nodes.values():
                node.dist = .01

        tree.data["error"] = math.sqrt(scipy.dot(lse.resids, lse.resids)) / \
                                       sum(x.dist for x in tree.nodes.values())

        setBranchError(conf, tree, lse.resids, lse.paths, lse.edges, lse.topmat)
        
    if isDebug(DEBUG_MED):
        util.toc()


def setBranchError(conf, tree, errors, dists, edges, topmat):
    """Assigns an error to each branch"""
    
    # init errors to zero
    for node in tree.nodes.values():
        node.data["error"] = 0
    
    npaths = util.Dict(default=0)
    
    # add up errors
    for i in xrange(len(topmat)):
        for j in xrange(len(edges)):
            if topmat[i][j] == 1:
                gene1, gene2 = edges[j]
                if tree.nodes[gene2].parent == tree.nodes[gene1]:
                    gene1, gene2 = gene2, gene1
                
                tree.nodes[gene1].data["error"] += \
                        util.safediv(abs(errors[i]), dists[i], 0)
                npaths[gene1] += 1
                
    
    # normalize by number of paths through edge
    for node in tree.nodes.itervalues():
        if node.name in npaths:
            node.data["error"] /= npaths[node.name]
        


def findSplits(network, leaves):
    """DEPRECATED: use phylo.findAllBranchSplits()"""
    
    return phylo.findAllBranchSplits(network, leaves)


def getSplit(tree):
    """DEPRECATED: use phylo.findBranchSplits()"""
    
    return phylo.findBranchSplits(tree)
              

def robinsonFouldsError(tree1, tree2):
    """DEPRECATED: use phylo.robinsonFouldsError()"""
    
    return phylo.robinsonFouldsError(tree1, tree2)



#=============================================================================
# Learning
#

def variance2(vals, u):
    return sum(map(lambda x: (x - u)**2, vals)) / float(len(vals)-1)

def sdev2(vals, u):
    return math.sqrt(variance2(vals, u))


def fitNormal2(lens):
    mu = stats.mean(lens)
    sigma = stats.sdev(lens)
    param, resid = stats.fitDistrib(stats.normalPdf, 
                                    [mu, sigma],
                                    lens,
                                    mu - 2*sigma,
                                    mu + 2*sigma,
                                    sigma / min(30, len(lens)/5))
    return param


def mleNormal(lens):
    mu = stats.mean(lens)
    sigma = stats.sdev(lens)
    return mu, sigma




def learnModel(trees, stree, gene2species, statsprefix="", filenames=None):
    util.tic("learn model")

    util.tic("find branch length distributions")
    lengths, used = phylo.findBranchDistrib(trees, stree, gene2species, False)
    debug("Total trees matching species topology: %d out of %d" % 
          (sum(used), len(trees)))
    util.toc()
    
    params = {}
    
    totlens = map(sum, zip(* lengths.values()))
    
    # print output stats
    if statsprefix != "":
        writeTreeDistrib(file(statsprefix + ".lens", "w"), lengths)
        rates = treeDistrib2table(lengths, filenames=filenames)
        rates.write(statsprefix + "_rates.tab")
    
    
    util.tic("fitting params")
    for node, lens in lengths.items():
        if len(lens) == 0 or max(lens) == min(lens):
            continue
        
        util.tic("fitting params for " + str(node.name))
        param = fitNormal2(util.vdiv(lens, totlens))
        
        params[node.name] = param
        util.toc()
    util.toc()
    
    # calc distribution of total tree length
    trees2 = util.mget(trees, util.findeq(True, used))
    lens = map(lambda x: sum(y.dist for y in x.nodes.values()), trees2)
    lens = filter(lambda x: x < 20, lens)
    mu = stats.mean(lens)
    lens = filter(lambda x: x < 2*mu, lens)
    mu = stats.mean(lens)
    sigma2 = stats.variance(lens)
    params["baserate"] = [mu*mu/sigma2, mu/sigma2]
    params[stree.root.name] = [0, 1]
    
    util.toc()
    
    return params




#=============================================================================
# new C-extension
#

def makePtree(tree):
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
        if a.isLeaf():
            if b.isLeaf():
                return 0
            else:
                return -1
        else:
            if b.isLeaf():
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


def treeLikelihood_C(conf, tree, recon, events, stree, params, generate, 
                     gene2species):
    """calculate likelihood of tree using C"""


    ptree, nodes, nodelookup = makePtree(tree)
    dists = [float(node.dist) for node in nodes]
    
    pstree, snodes, snodelookup = makePtree(stree)
    
    reconarray = [snodelookup[recon[node]] for node in nodes]
    eventslookup = {"gene": Likelihood.EVENT_GENE,
                    "spec": Likelihood.EVENT_SPEC,
                    "dup": Likelihood.EVENT_DUP}
    eventsarray = [eventslookup[events[node]] for node in nodes]
    
    gene2speciesarray = []
    for node in nodes:
        if node.isLeaf():
            gene2speciesarray.append(snodelookup[
                                     stree.nodes[gene2species(node.name)]])
        else:
            gene2speciesarray.append(-1)
    
    mu = [float(params[snode.name][0]) for snode in snodes]
    sigma = [float(params[snode.name][1]) for snode in snodes]
    
        
    ret = pyspidir.treelk(ptree, dists, pstree, 
                           gene2speciesarray,
                           mu, sigma, 
                           float(params["baserate"][0]), 
                           float(params["baserate"][1]),
                           generate,
                           float(conf.get("predupprob", 1.0)),
                           float(conf.get("dupprob", 1.0)))
    
    return ret


def parsimony(aln, usertree):
    return parsimony_C(aln, usertree)

def parsimony_C(aln, tree):    
    ptree, nodes, nodelookup = makePtree(tree)
    leaves = [x.name for x in nodes if isinstance(x.name, str)]
    seqs = util.mget(aln, leaves)
    
    dists = pyspidir.parsimony(ptree, seqs)
    
    for i in xrange(len(dists)):
        nodes[i].dist = dists[i]


def mlhkydist(aln, tree, bgfreq, ratio, maxiter):
    return mlhkydist_C(aln, tree, bgfreq, ratio, maxiter)


def mlhkydist_C(aln, tree, bgfreq, ratio, maxiter):
    ptree, nodes, nodelookup = makePtree(tree)
    leaves = [x.name for x in nodes if isinstance(x.name, str)]
    seqs = util.mget(aln, leaves)
    
    dists, logl = pyspidir.mlhkydist(ptree, seqs, bgfreq, ratio, maxiter)
    
    for i in xrange(len(dists)):
        nodes[i].dist = dists[i]
    
    return logl


def estGeneRate(tree, stree, params, gene2species):
    """calculate maximum a posteriori gene rate"""

    ptree, nodes, nodelookup = makePtree(tree)
    dists = [float(node.dist) for node in nodes]
    
    pstree, snodes, snodelookup = makePtree(stree)
    
    gene2speciesarray = []
    for node in nodes:
        if node.isLeaf():
            gene2speciesarray.append(snodelookup[
                                     stree.nodes[gene2species(node.name)]])
        else:
            gene2speciesarray.append(-1)
    
    mu = [float(params[snode.name][0]) for snode in snodes]
    sigma = [float(params[snode.name][1]) for snode in snodes]
    
        
    ret = pyspidir.est_generate(ptree, dists, pstree, 
                                gene2speciesarray,
                                mu, sigma, 
                                float(params["baserate"][0]), 
                                float(params["baserate"][1]))
    
    return ret
    

#=============================================================================
# branch length likelihood function

def treeLogLikelihood(conf, tree, stree, gene2species, params, baserate=None):
    conf.setdefault("bestlogl", -util.INF)
    
    if pyspidir == None or conf.get("python_only", False):
        return Likelihood.treeLogLikelihood_python(conf, tree, stree, 
                                                   gene2species, params, 
                                                   baserate=baserate, 
                                                   integration="fastsampling")

    # debug info
    if isDebug(DEBUG_MED):
        util.tic("find logl")
    

    # derive relative branch lengths
    #tree.clearData("logl", "extra", "fracs", "params", "unfold")
    recon = phylo.reconcile(tree, stree, gene2species)
    events = phylo.labelEvents(tree, recon)
    
    # determine if top branch unfolds
    if recon[tree.root] ==  stree.root and \
       events[tree.root] == "dup":
        for child in tree.root.children:
            if recon[child] != stree.root:
                child.data["unfold"] = True

    # top branch is "free"
    params[stree.root.name] = [0,0]
    this = util.Bundle(logl=0.0)
    
    if conf.get("generate_int", False):
        baserate = -99.0 # indicates in integration over gene rates is requested
    elif baserate == None:
        baserate = Likelihood.getBaserate(tree, stree, params, recon=recon)
        
    
    phylo.midrootRecon(tree, stree, recon, events, params, baserate)
    
    # calc likelihood in C
    this.logl = treeLikelihood_C(conf, tree, recon, events, stree, params, 
                                 baserate, gene2species)
    
    # calc probability of rare events
    tree.data["eventlogl"] = Likelihood.rareEventsLikelihood(conf, tree, stree, recon, events)
    
    # calc penality of error
    tree.data["errorlogl"] = tree.data.get("error", 0.0) * \
                             conf.get("errorcost", 0.0)
    this.logl += tree.data["errorlogl"]
    
    # add logl of sequence evolution
    this.logl += tree.data.get("distlogl", 0.0)
    
    if baserate == -99.0: # indicates in integration over gene rates is requested
        baserate = Likelihood.getBaserate(tree, stree, params, recon=recon) 
    
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    if isDebug(DEBUG_MED):
        util.toc()
        debug("\n\n")
        drawTreeLogl(tree, events=events)
        
    return this.logl



#-------------------------------------------------------------------------------
# Main SPIDIR algorithm function
#-------------------------------------------------------------------------------

def spidir(conf, distmat, labels, stree, gene2species, params):
    """Main function for the SPIDIR algorithm"""
    
    setDebug(conf["debug"])

    if isDebug(DEBUG_HIGH) and pyspidir:
        pyspidir.set_log(3, "")
        
    
    if "out" in conf:
        # create debug table
        conf["debugtab_file"] = file(conf["out"] + ".debug.tab", "w")
        
        debugtab = tablelib.Table(headers=["correct",
                                           "logl", "treelen", "baserate", 
                                           "error", "errorlogl", 
                                           "eventlogl", "tree",
                                           "topology", "species_hash"],
                                  types={"correct": bool,
                                         "logl": float, 
                                         "treelen": float, 
                                         "baserate": float, 
                                         "error": float, 
                                         "errorlogl": float,
                                         "eventlogl": float, 
                                         "tree": str,
                                         "topology": str,
                                         "species_hash": str})
        debugtab.writeHeader(conf["debugtab_file"])
        conf["debugtab"] = debugtab
    else:
        conf["debugfile"] = None
    
    
    trees = []
    logls = []
    tree = None
    visited = {}
    
    util.tic("SPIDIR")
    
    # do auto searches
    for search in conf["search"]:
        util.tic("Search by %s" % search)
        
        if search == "greedy":
            tree, logl = Search.searchGreedy(conf, distmat, labels, stree, 
                                      gene2species, params,
                                      visited=visited)
            
        elif search == "mcmc":
            tree, logl = Search.searchMCMC(conf, distmat, labels, stree, 
                                    gene2species, params, initTree=tree,
                                    visited=visited)
                                    
        elif search == "regraft":
            tree, logl = Search.searchRegraft(conf, distmat, labels, stree, 
                                    gene2species, params, initTree=tree,
                                    visited=visited, proposeFunc=Search.proposeTree3)
                                    
        elif search == "exhaustive":
            if tree == None:
                tree = phylo.neighborjoin(distmat, labels)
                tree = phylo.reconRoot(tree, stree, gene2species)
            
            tree, logl = Search.searchExhaustive(conf, distmat, labels, tree, stree, 
                                          gene2species, params, 
                                          depth=conf["depth"],
                                          visited=visited)
        elif search == "hillclimb":
            tree, logl = Search.searchHillClimb(conf, distmat, labels, stree, 
                                         gene2species, params, initTree=tree,
                                         visited=visited)
        
        elif search == "none":
            break
        else:
            raise SindirError("unknown search '%s'" % search)
        
        util.toc()
        
        Search.printMCMC(conf, "N/A", tree, stree, gene2species, visited)
        
        printVisitedTrees(visited)
        

    def evalUserTree(tree):        
        setTreeDistances(conf, tree, distmat, labels)
        logl = treeLogLikelihood(conf, tree, stree, gene2species, params)
        
        thash = phylo.hashTree(tree)
        if thash in visited:
            a, b, count = visited[thash]
        else:
            count = 0
        visited[thash] = [logl, tree.copy(), count+1]
        
        if isDebug(DEBUG_LOW):
            debug("\nuser given tree:")
            recon = phylo.reconcile(tree, stree, gene2species)
            events = phylo.labelEvents(tree, recon)
            drawTreeLogl(tree, events=events)        
    
    # eval the user given trees
    for treefile in conf["tree"]:
        tree = treelib.readTree(treefile)
        evalUserTree(tree)
    
    for topfile in conf["tops"]:
        infile = file(topfile)
        strees = []
        
        while True:
            try:
                strees.append(treelib.readTree(infile))
            except:
                break
        
        print len(strees)
        
        for top in strees:
            tree = phylo.stree2gtree(top, labels, gene2species)
            evalUserTree(tree)    
    
    if len(conf["tops"]) > 0:
        printVisitedTrees(visited)    
    
    
    
    # eval correcttree for debug only
    if "correcttree" in conf:
        tree = conf["correcttree"]
        setTreeDistances(conf, tree, distmat, labels)
        logl = treeLogLikelihood(conf, tree, stree, gene2species, params)
        
        if isDebug(DEBUG_LOW):
            debug("\ncorrect tree:")
            recon = phylo.reconcile(tree, stree, gene2species)
            events = phylo.labelEvents(tree, recon)
            drawTreeLogl(tree, events=events)
    
    
    util.toc()
    
    if len(visited) == 0:
        raise SindirError("No search or tree topologies given")
    
    
    if "correcthash" in conf:
        if conf["correcthash"] in visited:
            debug("SEARCH: visited correct tree")
        else:
            debug("SEARCH: NEVER saw correct tree")

    
    # return ML tree
    trees = [x[1] for x in visited.itervalues()]
    i = util.argmax([x.data["logl"] for x in trees])
    return trees[i], trees[i].data["logl"]
    



#
# testing
#
if __name__ == "__main__":
    import StringIO
    
    def floateq(a, b, accuracy=.0001):
        if b - accuracy <= a <= b + accuracy:
            print "pass"
            print a, "==", b
        else:
            print a, "!=", b
            raise Exception("not equal")
        
    
    
    def gene2species(name):
        return name[:1].upper()
    
    
    params = {"A": [4, 2],
              "B": [3, 1]}
              
    conf = {"debug": 0,
            "dupprob": .5,
            "lossprob": 1.0}
    
    
    
    stree = treelib.readTree(StringIO.StringIO("(A, B);"))
    
    
    # test 1
    print "\n\nTest 1"
    tree  = treelib.readTree(StringIO.StringIO("(a:3, b:2);"))
    logl = treeLogLikelihood(conf, tree, stree, gene2species, params, baserate=1)
    
    treelib.drawTreeLens(tree,scale=5)
    floateq(logl, log(stats.normalPdf(3, params["A"]) *
                      stats.normalPdf(2, params["B"])))
    
    
    # test 2
    print "\n\nTest 2"    
    tree  = treelib.readTree(StringIO.StringIO("((a1:2.5, a2:2):1, b:2);"))
    logl = treeLogLikelihood(conf, tree, stree, gene2species, params, baserate=1)
    
    treelib.drawTreeLens(tree,scale=5)
    floateq(logl, log(stats.normalPdf(2.5+1, params["A"])) +
                  log(stats.normalPdf(2+1, params["A"])) -
                  log(1.0 - stats.normalCdf(1, params["A"])) +
                  log(stats.normalPdf(2, params["B"])))


    print "\n\nTest 3"    
    tree  = treelib.readTree(StringIO.StringIO(
                             "(((a1:2.5, a2:2):1, a3:1.5):1.2, b:2);"))
    logl = treeLogLikelihood(conf, tree, stree, gene2species, params, baserate=1)
    
    treelib.drawTreeLens(tree,scale=5)
    floateq(logl, log(stats.normalPdf(2.5+1+1.2, params["A"])) +
                  log(stats.normalPdf(2+1+1.2, params["A"])) -
                  log(1.0 - stats.normalCdf(1+1.2, params["A"])) +
                  log(stats.normalPdf(1.5+1.2, params["A"])) -
                  log(1.0 - stats.normalCdf(1.2, params["A"])) +
                  log(stats.normalPdf(2, params["B"])))    




