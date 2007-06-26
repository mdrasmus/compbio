#
# SPIDIR library
#
# note: SPIDIR was a codename, it may still be present in the code.
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

from rasmus.vis import treevis

# scipy libs
# (needed for numerical integration and least square error fitting)
import scipy
import scipy.linalg
import scipy.integrate
import scipy.optimize


# SPIDIR libs
from Spidir import Search
from Spidir.Debug import *
from Spidir import pyspidir


# events
EVENT_GENE = 0
EVENT_SPEC = 1
EVENT_DUP = 2

# fractional branches
FRAC_NONE = 0
FRAC_DIFF = 1
FRAC_PARENT = 2
FRAC_NODE = 3



#-------------------------------------------------------------------------------
# SINDIR input/output
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
    canvas = treevis.drawTree(tree, autoclose=False, 
                              tmargin=tmargin, lmargin=lmargin,
                              *args, **kargs)
    
    # draw variance
    coords = treevis.layoutTree(tree, kargs["xscale"], kargs["yscale"],
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
    
    if "parsimony" in conf:
        parsimony_C(conf["aln"], tree)
        tree.data["error"] = 0.0
        return
    
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



#-------------------------------------------------------------------------------
# Learning
#-------------------------------------------------------------------------------

def variance2(vals, u):
    return sum(map(lambda x: (x - u)**2, vals)) / float(len(vals)-1)

def sdev2(vals, u):
    return math.sqrt(variance2(vals, u))

def fitNormal(lens):
    ndivs = int((max(lens) - min(lens)) / stats.mean(lens) *  40)
    x, y = util.hist(lens, ndivs)
    mu = x[util.argmax(y)]
    data = filter(util.withinfunc(0, mu*2), lens)

    if len(data) < 2:
        sigma = 1
    else:
        sigma = sdev2(data, mu)
    
    return mu, sigma

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



def fitParams(lengths, baserates, gene2species, fit=True):
    ntrees = len(lengths.values()[0])
    
    params = {}
    
    dist = util.distrib(baserates, width=.01)
    top = min(max(baserates, 10))
    param, resid = stats.fitCurve(dist[0], dist[1], stats.gammaPdf, [0,top])
    params["baserate"] = param
    
    
    util.tic("fitting params")
    for node, lens in lengths.items():
        if len(lens) == 0 or max(lens) == min(lens):
            params[node.name] = [0, 1]
            continue
        
        #util.tic("fitting params for " + str(node.name))
        
        lens = util.vdiv(lens, baserates)
        
        if fit:
            mu = stats.mean(lens)
            sigma = stats.sdev(lens)            
            param, resid = stats.fitDistrib(stats.normalPdf, 
                                            [mu, sigma],
                                            lens,
                                            mu - 3*sigma,
                                            mu + 3*sigma,
                                            sigma / 10)
            params[node.name] = param 
        else:
            ndivs = int((max(lens) - min(lens)) / stats.mean(lens) *  40)
            x, y = util.hist(lens, ndivs)
            mu = x[util.argmax(y)]
            data = filter(util.withinfunc(0, mu*2), lens)
            
            if len(data) < 2:
                sigma = 1
            else:
                sigma = sdev2(data, mu)
            
            #mu = stats.mean(lens)
            #sigma = stats.sdev(lens)
            
            params[node.name] = [mu, sigma]
        
        #util.toc()
    util.toc()
    
    return params


def dataLikelihood(lenmat, baserates, means, sdevs, baserateparam):
    logl = 0
    
    for i in range(len(lenmat)):
        for j in range(len(lenmat[i])):
            logl += log(stats.normalPdf(lenmat[i][j]/baserates[i], 
                                        [means[j], sdevs[j]]))
    
        # calc baserate logl
        logl += stats.gammaPdf(baserates[i], baserateparam)
    
    return logl


def dataLikelihoodG(lenmat, baserates, alphas, betas, baserateparam):
    logl = 0
    
    for i in range(len(lenmat)):
        for j in range(len(lenmat[i])):
            logl += log(stats.gammaPdf(lenmat[i][j]/baserates[i], 
                                        [alphas[j], betas[j]]))
    
    return logl


def mleBaserates(lengths, params, baserateparam):
    lenmat = zip(* lengths.values())
    keys = map(lambda x: x.name, lengths.keys())
    means, sdevs = zip(* util.mget(params, keys))
    baserates = []
    for i in xrange(len(lenmat)):
        baserates.append(mleBaserate(lenmat[i], means, sdevs, baserateparam))
    return baserates


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
    
    
def mleGamma(vals):
    mu = stats.mean(vals)
    sigma2 = stats.variance(vals)
    return [mu*mu/sigma2, mu/sigma2]


def learnModel2(lengths, gene2species, niters=10):
    lenmat = zip(* lengths.values())
    
    debug = util.Closure(
        baserates=[],
        params=[],
        logls=[])
    
    # init base rates
    baserates = map(sum, lenmat)
    debug.baserates.append(baserates)
    mu = stats.mean(baserates)
    sigma2 = stats.variance(baserates)
    baserateparam = [mu*mu/sigma2, mu/sigma2]
    
    
    # do EM
    for i in range(niters):
        print i
        params = {}
        for name, lens in lengths.items():
            if len(lens) == 0 or max(lens) == min(lens):
                continue
            #params[name] = fitNormal2(util.vdiv(lens, baserates))
            params[name] = mleNormal(util.vdiv(lens, baserates))
            #params[name] = mleGamma(util.vdiv(lens, baserates))
        debug.params.append(params)
        
        means, sdevs = zip(* params.values())
        
        baserates = []
        for j in xrange(len(lenmat)):
            baserates.append(mleBaserate3(lenmat[j], means, sdevs, 
                                          baserateparam))
        debug.baserates.append(baserates)
        
        # calc likelihood
        logl = dataLikelihood(lenmat, baserates, means, sdevs, 
                                baserateparam)
        util.log("log:", logl,
                 "mean baserate:", stats.mean(baserates))
        debug.logls.append(logl)

        
    
    return debug



#-------------------------------------------------------------------------------
# Likelihood calculation
#-------------------------------------------------------------------------------


def mleBaserate2(lens, means, sdevs, baserateparam):
    ind = range(len(means))
    ind.sort(lambda a, b: cmp(means[b], means[a]))
    ind = ind[:len(ind) / 2 + 1]
    means = util.mget(means, ind)
    sdevs = util.mget(sdevs, ind)
    lens = util.mget(lens, ind)
    
    vars = util.vmul(sdevs, sdevs)
    denom = sum(util.vdiv(util.vmul(means, lens), vars))
    
    if denom != 0:
        return sum(util.vdiv(util.vmul(lens, lens), vars)) / denom
    else:
        return 1.0
           


def mleBaserate(lens, means, sdevs, baserateparam):
    [alpha, beta] = baserateparam
    
    # use only best means and sdevs (highest means)
    
    """
    ind = range(len(means))
    ind.sort(lambda a, b: cmp(means[b], means[a]))
    ind = ind[:max(4, len(ind) / 2 + 1)]
    #ind = ind[:max(4, len(means)-2)]
    means = util.mget(means, ind)
    sdevs = util.mget(sdevs, ind)
    lens = util.mget(lens, ind)
    """
    
    # protect against zero
    ind = util.findgt(.0001, sdevs)
    lens = util.mget(lens, ind)
    means = util.mget(means, ind)
    sdevs = util.mget(sdevs, ind)
    
    a = (1 - alpha) / beta
    b = sum(means[i] * lens[i] / sdevs[i]**2
            for i in range(len(lens))) / beta
    c = - sum(lens[i] ** 2 / sdevs[i] ** 2
              for i in range(len(lens))) / beta
    
    #p = util.plotfunc(lambda x: x*x*x + a+x*x + b*x + c, -50, 50, .1)
    #p.plot([-50, 50], [0, 0], style="lines")
    #time.sleep(10)
    
    #print filter(lambda x: x>0, stats.solveCubic(a, b, c))
    return max(stats.solveCubic(a, b, c))


def mleBaserate3(lens, means, sdevs, baserateparam):
    [alpha, beta] = baserateparam
    
    # use only best means and sdevs (highest means)
    
    ind = range(len(means))
    ind.sort(lambda a, b: cmp(means[b], means[a]))
    ind = ind[:len(ind) / 4 + 1]
    means = util.mget(means, ind)
    sdevs = util.mget(sdevs, ind)
    lens = util.mget(lens, ind)
    
    
    # protect against zero
    ind = util.findgt(.0001, sdevs)
    lens = util.mget(lens, ind)
    means = util.mget(means, ind)
    sdevs = util.mget(sdevs, ind)
    
    a = (1 - alpha) / beta
    b = sum(means[i] * lens[i] / sdevs[i]**2
            for i in range(len(lens))) / beta
    c = - sum(lens[i] ** 2 / sdevs[i] ** 2
              for i in range(len(lens))) / beta
    
    #print filter(lambda x: x>0, stats.solveCubic(a, b, c))
    return max(stats.solveCubic(a, b, c))


def mleBaserateG(lens, alphas, betas, baserateparam):
    [alpha, beta] = baserateparam
    
    # protect against zero
    #ind = util.findgt(.0001, sdevs)
    #lens = util.mget(lens, ind)
    #alphas = util.mget(means, ind)
    #betas = util.mget(sdevs, ind)
    
    nom = 0
    denom = 0
    
    for i in xrange(len(lens)):
        nom += lens[i] * betas[i]
        denom += alphas[i] - 1
    
    b = nom / denom
    return b



def getExtraBranches(root, recon, events, stree):
    extraBranches = {}

    # determine if any extra branches exist
    def markExtras(node):
        if recon[node] == stree.root and \
           events[node] == "dup":
            for child in node.children:
                if recon[child] != stree.root:
                    extraBranches[child] = 1
                    child.data["extra"] = 1
        node.recurse(markExtras)
    markExtras(root)
     
    return extraBranches


def getBaserate(tree, stree, params, recon=None, gene2species=None):
    if recon == None:
        assert gene2species != None
        recon = phylo.reconcile(tree, stree, gene2species)
    events = phylo.labelEvents(tree, recon)
    
    extraBranches = getExtraBranches(tree.root, recon, events, stree)
    
    lens = []
    means = []
    sdevs = []
    
    # process each child of subtree root
    def walk(node, depths, sroot, extra):
        # save depth of node
        if recon[node] != recon[tree.root]:  #stree.root:
            depths[node] = node.dist + depths[node.parent]
        else:
            # ignore branch length of free branches
            depths[node] = depths[node.parent]
        
        
        # record presence of extra in path
        extra = extra or ("extra" in node.data)
        
        
        if events[node] == "dup":
            # recurse within dup-only subtree
            #   therefore pass depths and sroot unaltered
            node.recurse(walk, depths, sroot, extra)
        else:
            # we are at subtree leaf
            
            # figure out species branches that we cross
            # get total mean and variance of this path            
            mu = 0
            sigma2 = 0            
            snode = recon[node]
            
            # branch is also free if we do not cross any more species
            # don't estimate baserates from extra branches
            if snode != sroot and not extra:
                
                while snode != sroot and snode != stree.root:
                    mu += params[snode.name][0]
                    sigma2 += params[snode.name][1]**2
                    snode = snode.parent
                assert abs(sigma2) > .00000001, "sigma too small"
                sigma = math.sqrt(sigma2)
                
                # save dist and params
                lens.append(depths[node])
                means.append(mu)
                sdevs.append(sigma)
            
            # continue recursion, but with new depths and sroot
            for child in node.children:
                walk(child, depths={node: 0}, sroot=recon[node], extra=False)
    
    
    for child in tree.root.children:
        walk(child, depths={tree.root: 0}, sroot=recon[tree.root], extra=False)
    
    
    baserate = mleBaserate(lens, means, sdevs, params["baserate"])
    
    #util.printcols(zip(lens, means, sdevs, util.vdiv(lens, means)))
    
    #baserate = mleBaserate2(lens, means, sdevs, params["baserate"])
        
    return baserate



def rareEventsLikelihood(conf, tree, stree, recon, events):
    logl = 0.0
    
    for node, event in events.items():
        if recon[node] == stree.root and \
           event == "dup":
            logl += log(conf["predupprob"])
        
        if event == "dup":
            logl += log(conf["dupprob"])
        
        #if event == "dup":
        #    logl += log(conf["dupprob"])
        
    #nloss = len(phylo.findLoss(tree, stree, recon))
    #logl += nloss * log(conf["lossprob"])
    
    return logl




### TESTING INTEGRATION

def branchLikelihood(dist, fracs, params):
    totmean = 0.0
    totvar = 0.0
    
    for frac, (mu, sigma) in zip(fracs, params):
        #assert frac > 0.0, frac
        frac = max(0.0, frac)  # ensure frac greater than zero
        totmean += frac * mu
        totvar  += frac * sigma * sigma
    
    # don't let variance get too low
    if totvar < .0000000001:
        return -util.INF
    
    return log(stats.normalPdf(dist, [totmean, math.sqrt(totvar)]))



def setMidpointsRandom(node, events, recon, midpoints, wholeTree=True):
    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup" and \
           node.parent != None:
            if recon[node] == recon[node.parent]:
                # if im the same species branch as my parent 
                # then he is my last midpoint
                lastpoint = midpoints[node.parent]
            else:
                # im the first on this branch so the last midpoint is zero
                lastpoint = 0.0
            
            # pick a midpoint uniformly after the last one
            midpoints[node] = lastpoint + \
                        (random.random() * (1 - lastpoint))
        else:
            # genes or speciations reconcile exactly to the end of the branch
            # gene tree roots also reconcile exactly to the end of the branch
            midpoints[node] = 1.0
    
        # recurse within dup-only subtree
        if events[node] == "dup" or wholeTree:
            node.recurse(walk)
    
    walk(node)


def countMidpointParameters(node, events):
    this = util.Closure(nvars = 0)

    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup":
            this.nvars += 1
        
        # recurse within dup-only subtree
        if events[node] == "dup":
            node.recurse(walk)
    
    walk(node)
    
    return this.nvars



def reconBranch(node, recon, events, stree, params, 
                midpoints, extraBranches):
    fracs = []
    snodes = []

    if recon[node] == recon[node.parent]:
        # we begin and end on same branch
        fracs.append(midpoints[node] - midpoints[node.parent])
        snodes.append(recon[node].name)
    else:
        # we begin and end on different branches

        # ending (most recent)
        if recon[node] != stree.root:
            fracs.append(midpoints[node])
            snodes.append(recon[node].name)

        # walk up until starting species branch
        snode = recon[node].parent
        snode2 = recon[node.parent]
        while snode != snode2:
            fracs.append(1) # full branch
            snodes.append(snode.name)
            snode = snode.parent

        # beginning branch
        if midpoints[node.parent] != 1.0 and \
           snode2 != stree.root:
            fracs.append(1.0 - midpoints[node.parent])
            snodes.append(snode2.name)


    # dealing with free mutations (extra branches)
    if "unfold" in node.data:
        dist = 2 * node.dist
    else:
        dist = node.dist
    
    return dist, fracs, snodes



def calcSubtreeLikelihood(root, recon, events, stree, params, 
                          midpoints, extraBranches, baserate):
    this = util.Bundle(
        logl=0.0
    )
    
    def walk(node):
        # branch that recon to the species root is "free"
        if recon[node] != stree.root:
            dist, fracs, snodes = reconBranch(node, recon, events, stree, 
                                     params, midpoints, extraBranches)
            params2 = [params[x] for x in snodes]

            dist /= baserate
            
            if "extra" in node.data and dist > 0:
                mean = util.vdot(fracs, [params[x][0] for x in snodes])
                dist = min(dist, mean)
            
            # add likelihood of branch to subtree total
            this.logl += branchLikelihood(dist, fracs, params2)
            
            
            # debug saving
            node.data["logl"] = "x"
            node.data["params"] = params2
            if "fracs" not in node.data:
                node.data["fracs"] = []
            node.data["fracs"].append(fracs)           
            
            
        # recurse within dup-only subtree
        if events[node] == "dup":
            node.recurse(walk)
    walk(root)
    
    return this.logl


def setMidpoints(node, events, recon, midpoints, kvars):
    this = util.Bundle(i = 0)
    
    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup":
            #if recon[node] == recon[node.parent]:
            #    # if im the same species branch as my parent 
            #    # then he is my last midpoint
            #    lastpoint = midpoints[node.parent]
            #else:
            #    # im the first on this branch so the last midpoint is zero
            #    lastpoint = 0.0
            
            # pick a midpoint uniformly after the last one
            #midpoints[node] = lastpoint + \
            #            (kvars[this.i] * (1 - lastpoint))
            midpoints[node] = kvars[this.i]
            this.i += 1
        else:
            # genes or speciations reconcile exactly to the end of the branch
            midpoints[node] = 1.0
        
        # recurse within dup-only subtree
        if events[node] == "dup":
            node.recurse(walk)
    
    walk(node)


def makeKDepend(node, events, recon):
    kdepend = {}
    korder = {}
    korderrev = {}
    
    this = util.Closure(i = 0)

    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup":
            if recon[node] == recon[node.parent]:
                # if im the same species branch as my parent 
                # then he is my last midpoint
                kdepend[node] = node.parent
            else:
                # im the first on this branch so the last midpoint is zero
                kdepend[node] = None
                
            korder[node] = this.i
            korderrev[this.i] = node
            this.i += 1
    
            # recurse within dup-only subtree
            node.recurse(walk)
    
    walk(node)
    
    return kdepend, korder, korderrev


def subtreeLikelihood(conf, root, recon, events, stree, params, baserate,
                      integration="quad"):
    midpoints = {}
    extraBranches = getExtraBranches(root, recon, events, stree)
    
    this = util.Closure(ncalls=0, depth=0, printing=0)    
    
    
    #quad = scipy.integrate.quad
    
    def quad(func, start, end, n):
        def func2(xs):
            return map(func, xs)
        val, err = scipy.integrate.fixed_quad(func2, start, end, n=n)
        err = 0
        return (val, err)
        
    """
    def quad2(func, start, end, epsrel=None):
        tot = 0
        n = 0
        if end <= start:
            return [0, 666]
        step = (end - start) / 20.
        for i in util.frange(start, end, step):
            tot += func(i) * step
        
        return [tot, 666]
    """
    
    print "INT", integration
    
    def integrate(node, nvars):
        kvars = [0] * nvars
        
        kdepend, korder, korderrev = makeKDepend(node, events, recon)
        def nodename(node):
            if node:
                return node.name
            else:
                return None
        
        ## debugging integration
        #print util.mapdict(kdepend, 
        #                   keyfunc=nodename,
        #                   valfunc=nodename)
        #print util.mapdict(korder, 
        #                   keyfunc=nodename)
        #print util.mapdict(korderrev, 
        #                   valfunc=nodename)
        
        
        def func(k):
            kvars[this.depth] = k
            this.depth += 1
            
            
            if this.depth == nvars:
                this.ncalls += 1
                setMidpoints(child, events, recon, midpoints, kvars)
                this.depth -= 1
                
                #if this.printing < 0:
                #    debug("int:", kvars)
                #    this.printing = 100
                #else:
                #    this.printing -= 1
                
                return math.e** \
                 calcSubtreeLikelihood(node, recon, events, stree, params, 
                                       midpoints, extraBranches, baserate)
            else:
                depnode = kdepend[korderrev[this.depth]]                
                if depnode == None:
                    lowk = 0
                else:
                    lowk = kvars[korder[depnode]] + .0001
                
                n = max(1, 10 / (2**this.depth))
                #n = 20
                
                ret = quad(func, lowk, .9999, n=n)[0] / \
                      (1.0 - lowk)
                this.depth -= 1
                return ret
        return quad(func, .0001, .9999, n=20)
    
    
    if integration == "quad":
        # process each child of subtree root
        logl = 0.0
        midpoints[root] = 1.0
        for child in root.children:
            # integration is only needed if child is dup
            if events[child] != "dup":
                setMidpoints(child, events, recon, midpoints, [])
                clogl = calcSubtreeLikelihood(child, recon, events, stree, params, 
                          midpoints, extraBranches, baserate)
                child.data["logl"] = clogl
                logl += clogl
            else:
                # integrate over midpoints
                nvars = countMidpointParameters(child, events)
                val, err = integrate(child, nvars)
                #if conf["debug"]:
                #    #debug("integration error:", (log(err + val) - log(val)))
                #    debug("int nvars:", nvars)
                clogl = log(val)
                child.data["logl"] = clogl
                logl += clogl
    
    # do random sampling integration
    if integration == "sampling":        
        logl = 0.0

        midpoints[root] = 1.0
        for child in root.children:
            # integration is only needed if child is dup
            if events[child] != "dup":
                setMidpoints(child, events, recon, midpoints, [])
                clogl = calcSubtreeLikelihood(child, recon, events, stree, params, 
                          midpoints, extraBranches, baserate)
                child.data["logl"] = clogl
                logl += clogl
            else:
                for samples in [100]: #[50, 100, 250, 500]:        
                    val = 0.0            

                    for i in xrange(samples):
                        setMidpointsRandom(child, events, recon, midpoints, False)
                        val += math.exp(calcSubtreeLikelihood(child, recon, events, 
                                                         stree, params, 
                                                         midpoints, extraBranches, 
                                                         baserate))

                    clogl = log(val / float(samples))

                child.data["logl"] = clogl
                logl += clogl
    
    
    if integration == "fastsampling":
        # do fast integration
        logl3 = 0.0
        midpoints[root] = 1.0
        for child in root.children:
            # integration is only needed if child is dup
            if events[child] != "dup":
                #setMidpoints(child, events, recon, midpoints, [])
                #clogl = calcSubtreeLikelihood(child, recon, events, stree, params, 
                #          midpoints, extraBranches, baserate)
                
                node = child
                if recon[node] != stree.root:
                    startparams, startfrac, midparams, \
                        endparams, endfrac, kdepend = \
                        reconBranch2(node, recon, events, params)
                
                    setMidpoints(child, events, recon, midpoints, [])            
                    clogl = branchLikelihood2(node.dist / baserate, 
                                              node, midpoints, 
                                              startparams, startfrac,
                                              midparams, endparams, 
                                              endfrac)
                else:
                    clogl = 0.0
            else:
                startparams = {}
                startfrac = {}
                midparams = {}
                endparams = {}
                endfrac = {}
                kdepend = {}

                # recon subtree
                nodes = []
                def walk(node):
                    nodes.append(node)
                    startparams[node], startfrac[node], midparams[node], \
                        endparams[node], endfrac[node], kdepend[node] = \
                        reconBranch2(node, recon, events, params)

                    if events[node] == "dup":
                        for child in node.children:
                            walk(child)
                walk(child)
                

                for samples in [100]: #[50, 100, 250, 500]:        
                    val = 0.0            

                    for i in xrange(samples):
                        setMidpointsRandom2(child, events, recon, midpoints)                

                        val2 = 0.0
                        for node in nodes:
                            if recon[node] != stree.root:
                                v = branchLikelihood2(node.dist / baserate, 
                                              node, midpoints, 
                                              startparams[node], startfrac[node],
                                              midparams[node], endparams[node], 
                                              endfrac[node])
                                val2 += v
                                #print 'v', v

                        val += math.exp(val2)
                    clogl = log(val / float(samples))
                    
            child.data["logl"] = clogl
            logl3 += clogl 
        logl = logl3
    
    
    #logl=logl2
    #print "logl", logl, logl2, logl3
    
    return logl


def reconBranch2(node, recon, events, params):

    # set fractional branches
    if recon[node] == recon[node.parent]:
        # start reconciles to a subportion of species branch
        if events[node] == "dup":
            # only case k's are dependent
            startfrac = FRAC_DIFF # k[node] - k[node.parent]
            kdepend = node.parent
        else:
            startfrac = FRAC_PARENT # 1.0 - k[node.parent]
            kdepend = None
        startparams = params[recon[node].name]

        # there is only one frac
        endfrac = FRAC_NONE
        endparams = None
    else:
        kdepend = None

        if events[node.parent] == "dup":            
            # start reconciles to last part of species branch
            startfrac = FRAC_PARENT # 1.0 - k[node.parent]
            startparams = params[recon[node.parent].name]
        else:
            startfrac = FRAC_NONE
            startparams = None

        if events[node] == "dup":
            # end reconciles to first part of species branch
            endfrac = FRAC_NODE # k[node]
            endparams = params[recon[node].name]
        else:    
            # end reconcile to at least one whole species branch
            endfrac = FRAC_NONE
            endparams = None
    
    # set midparams
    if recon[node] == recon[node.parent]:
        # we begin and end on same branch
        # there are no midparams
        midparams = None
    else:
        # we begin and end on different branches
        totmean = 0.0
        totvar = 0.0

        # determine most recent species branch which we fully recon to
        if events[node] == "dup":
            snode = recon[node].parent
        else:
            snode = recon[node]

        # walk up species spath until starting species branch
        # starting species branch is either fractional or NULL
        parent_snode = recon[node.parent]
        while snode != parent_snode:
            totmean += params[snode.name][0]
            totvar += params[snode.name][1] ** 2
            snode = snode.parent

        midparams = [totmean, math.sqrt(totvar)]

    return startparams, startfrac, midparams, endparams, endfrac, kdepend


def branchLikelihood2(dist, node, k, startparams, startfrac,
                      midparams, endparams, endfrac):
    totmean = 0.0
    totvar  = 0.0
    
    #print k[node], startfrac, midparams, endfrac, endparams
    
    if startfrac == FRAC_DIFF:
        totmean += (k[node] - k[node.parent]) * startparams[0]
        totvar  += (k[node] - k[node.parent]) * startparams[1] ** 2
    elif startfrac == FRAC_PARENT:
        totmean += (1.0 - k[node.parent]) * startparams[0]
        totvar  += (1.0 - k[node.parent]) * startparams[1] ** 2
    #else startfrac == FRAC_NONE:
    #    pass
    
    if midparams != None:
        totmean += midparams[0]
        totvar  += midparams[1] ** 2
    
    if endfrac == FRAC_PARENT:
        totmean += (1.0 - k[node.parent]) * endparams[0]
        totvar  += (1.0 - k[node.parent]) * endparams[1] ** 2
    elif endfrac == FRAC_NODE:
        totmean += k[node] * endparams[0]
        totvar  += k[node] * endparams[1] ** 2
    #else endfrac == FRAC_NONE:
    #    pass
    
    if totvar <= 0.0:
        print "!!!!"
        print k[node], k[node.parent]
        print startfrac, startparams, midparams, endfrac, endparams
    
    
    # handle partially-free branches and unfold
    if "unfold" in node.data:
        dist *= 2;
    
    # augment a branch if it is partially free
    if "extra" in node.data:
        if dist > totmean:
            dist = totmean
    
    try:
        return log(stats.normalPdf(dist, [totmean, math.sqrt(totvar)]))
    except:
        print >>sys.stderr, dist, node.name, \
                  k, startparams, startfrac, midparams, endparams, endfrac
        raise



def setMidpointsRandom2(node, events, recon, midpoints, wholeTree=False):
    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup" and \
           node.parent != None:
            if recon[node] == recon[node.parent]:
                # if im the same species branch as my parent 
                # then he is my last midpoint
                lastpoint = midpoints[node.parent]
            else:
                # im the first on this branch so the last midpoint is zero
                lastpoint = 0.0
            
            # pick a midpoint uniformly after the last one
            midpoints[node] = lastpoint + \
                        (random.random() * (1 - lastpoint))
        else:
            # genes or speciations reconcile exactly to the end of the branch
            # gene tree roots also reconcile exactly to the end of the branch
            midpoints[node] = 1.0
    
        # recurse within dup-only subtree
        if events[node] == "dup" or wholeTree:
            node.recurse(walk)
    
    walk(node)


def treeLogLikelihood_old(conf, tree, stree, gene2species, params, baserate=None):
    #conf["accuracy"] = .01
    conf.setdefault("bestlogl", -util.INF)

    #estlogl = treeLogLikelihood_est(conf, tree, stree, gene2species, params, baserate=None)
    #if "integrate" not in conf:
    #    #return estlogl
    #    pass
    #else:
    
    # skip trees that are estimated to be very bad
    #if estlogl < conf["bestlogl"] - 50:
    #    debug("SKIP integration, est logl=%f" % estlogl)
    #    return estlogl
    #conf["bestlogl"] = max(conf["bestlogl"], estlogl)
            
    #Search.printMCMC(conf, "est", tree, stree, gene2species, {})


    # debug info
    if isDebug(DEBUG_MED):
        util.tic("find logl")
    

    # derive relative branch lengths
    tree.clearData("logl", "extra", "fracs", "params", "unfold")
    recon = phylo.reconcile(tree, stree, gene2species)
    events = phylo.labelEvents(tree, recon)

    # determine if top branch unfolds
    if recon[tree.root] ==  stree.root and \
       events[tree.root] == "dup":
        for child in tree.root.children:
            if recon[child] != stree.root:
                child.data["unfold"] = True
    
    if baserate == None:
        baserate = getBaserate(tree, stree, params, recon=recon)
    
    # top branch is "free"
    params[stree.root.name] = [0,0]
    this = util.Bundle(logl=0.0)
    
    # calc likelihood in C
    this.logl = treelikelihood_C(tree, recon, events, stree, params, baserate)
    
    
    # recurse through indep sub-trees
    #def walk(node):
    #    if events[node] == "spec" or \
    #       node == tree.root:
    #        this.logl += subtreeLikelihood(conf, node, recon, events, 
    #                                       stree, params, baserate)
    #    node.recurse(walk)
    #walk(tree.root)
    
    #print "cmp", this.logl, this.logl2

    
    # calc probability of rare events
    tree.data["eventlogl"] = rareEventsLikelihood(conf, tree, stree, recon, events)
    this.logl += tree.data["eventlogl"]
    
    # calc penality of error
    tree.data["errorlogl"] = tree.data["error"] * conf["errorcost"]
    this.logl += tree.data["errorlogl"]

    # family rate likelihood
    if conf["famprob"]:
        this.logl += log(stats.gammaPdf(baserate, params["baserate"]))
    
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    #print >> conf["intcmp"], "%f\t%f" % (this.logl, estlogl)
    #conf["intcmp"].flush()
    
    
    if isDebug(DEBUG_MED):
        util.toc()
        debug("\n\n")
        drawTreeLogl(tree, events=events)
    
    return this.logl





#=============================================================================
# new C-extension

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
    eventslookup = {"gene": EVENT_GENE,
                    "spec": EVENT_SPEC,
                    "dup": EVENT_DUP}
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
                           float(tree.data["error"]),
                           float(conf["predupprob"]),
                           float(conf["dupprob"]),
                           float(conf["errorcost"]))
    
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
    
    
    #treelib.drawTreeLens(tree)
    
    

def treeLogLikelihood(conf, tree, stree, gene2species, params, baserate=None):
    conf.setdefault("bestlogl", -util.INF)
    
    if "python_only" in conf and conf["python_only"]:
        return treeLogLikelihood_python(conf, tree, stree, gene2species, params, 
                                        baserate=baserate)

    if "lkcmp" in conf and conf["lkcmp"]:
        pysampling = treeLogLikelihood_python(conf, tree, stree, gene2species, params, 
                                     baserate=baserate, integration="sampling")

        pyfastsampling = treeLogLikelihood_python(conf, tree, stree, gene2species, params,
                                  baserate=baserate, integration="fastsampling")


        pyquad = treeLogLikelihood_python(conf, tree, stree, gene2species, params, 
                                          baserate=baserate, integration="quad")
    # debug info
    if isDebug(DEBUG_MED):
        util.tic("find logl")
    

    # derive relative branch lengths
    #tree.clearData("logl", "extra", "fracs", "params", "unfold")
    recon = phylo.reconcile(tree, stree, gene2species)
    events = phylo.labelEvents(tree, recon)
    
    phylo.midrootRecon(tree, stree, recon, events, params)
    
    # determine if top branch unfolds
    if recon[tree.root] ==  stree.root and \
       events[tree.root] == "dup":
        for child in tree.root.children:
            if recon[child] != stree.root:
                child.data["unfold"] = True
    
    if baserate == None:
        baserate = getBaserate(tree, stree, params, recon=recon)
    
    # top branch is "free"
    params[stree.root.name] = [0,0]
    this = util.Bundle(logl=0.0)
    
    # test out C parsimony
    if "aln" in conf and "parsimony" in conf:
        print "HERE"
        parsimony_C(conf["aln"], tree)
    
    
    # calc likelihood in C
    this.logl = treeLikelihood_C(conf, tree, recon, events, stree, params, 
                                 baserate, gene2species)
    
    # calc probability of rare events
    tree.data["eventlogl"] = rareEventsLikelihood(conf, tree, stree, recon, events)
    #this.logl += tree.data["eventlogl"]
    
    # calc penality of error
    tree.data["errorlogl"] = tree.data["error"] * conf["errorcost"]
    #this.logl += tree.data["errorlogl"]
    
    # generate
    #if conf["famprob"]:
    #    this.logl += log(stats.gammaPdf(baserate, params["baserate"]))    
    
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    if isDebug(DEBUG_MED):
        util.toc()
        debug("\n\n")
        drawTreeLogl(tree, events=events)
    
    if "lkcmp" in conf:
        tree.write(sys.stdout, oneline=True)
        print
        print "LKCMP\t%f\t%f\t%f\t%f" % \
                (pyquad, pysampling, pyfastsampling, this.logl)
    
    return this.logl


def treeLogLikelihood_python(conf, tree, stree, gene2species, params, 
                             baserate=None, integration="quad"):

    # debug info
    if isDebug(DEBUG_MED):
        util.tic("find logl")
    
    # derive relative branch lengths
    tree.clearData("logl", "extra", "fracs", "params", "unfold")
    recon = phylo.reconcile(tree, stree, gene2species)
    events = phylo.labelEvents(tree, recon)

    # determine if top branch unfolds
    if recon[tree.root] ==  stree.root and \
       events[tree.root] == "dup":
        for child in tree.root.children:
            if recon[child] != stree.root:
                child.data["unfold"] = True
    
    if baserate == None:
        baserate = getBaserate(tree, stree, params, recon=recon)
    
    # top branch is "free"
    params[stree.root.name] = [0,0]
    this = util.Closure(logl=0.0)
    
    # recurse through indep sub-trees
    def walk(node):
        if events[node] == "spec" or \
           node == tree.root:
            this.logl += subtreeLikelihood(conf, node, recon, events, 
                                           stree, params, baserate, 
                                           integration=integration)
        node.recurse(walk)
    walk(tree.root)
    
    
    # calc probability of rare events
    tree.data["eventlogl"] = rareEventsLikelihood(conf, tree, stree, recon, events)
    this.logl += tree.data["eventlogl"]
    
    # calc penality of error
    tree.data["errorlogl"] = tree.data["error"] * conf["errorcost"]
    this.logl += tree.data["errorlogl"]

    # family rate likelihood
    if conf["famprob"]:
        this.logl += log(stats.gammaPdf(baserate, params["baserate"]))
    
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    #print >> conf["intcmp"], "%f\t%f" % (this.logl, estlogl)
    #conf["intcmp"].flush()
    
    
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

    
    #if "integrate" in conf and "out" in conf:
    #    conf["intcmp"] = file(conf["out"] + ".int", "w")
        
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
        if True: #sum(node.dist for node in tree.nodes.values()) == 0.0: # or True:
            debug("fitting distances")     
            setTreeDistances(conf, tree, distmat, labels)
        else:
            debug("use distances from file")
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




