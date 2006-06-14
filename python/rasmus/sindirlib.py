#
# SINDIR library
#
#


# rasmus libs
#from rasmus import algorithms
from rasmus import fasta
from rasmus import matrix
from rasmus import phyloutil
from rasmus import stats
from rasmus import treelib
from rasmus import util

# python libs
import math, StringIO, copy, random, sys


# scipy libs
# (needed for numerical integration and least square error fitting)
import scipy
import scipy.linalg
import scipy.integrate
import scipy.optimize



#-------------------------------------------------------------------------------
# debugging variables and functions
#-------------------------------------------------------------------------------
DEBUG = sys.stdout


def debug(* text):
    print >>DEBUG, " ".join(map(str, text))


def setDebugStream(stream):
    globals()["DEBUG"] = stream


def drawTreeLogl(tree, out=DEBUG, events={}, baserate=1.0):
    labels = {}
    
    if "baserate" in tree.data:
        baserate = tree.data["baserate"]
    
    for node in tree.nodes.values():
        notes = ""
        if "extra" in node.data:
            notes += "E"
        if "unfold" in node.data:
            notes += "U"
        
        if "logl" in node.data:
            if isinstance(node.data["logl"], float):
                labels[node.name] = "[%s]\n%.3f (%.3f) %s" % \
                    (node.name, node.dist, node.data["logl"], notes)
                #logl += node.data["logl"]
            else:
                labels[node.name] = "[%s]\n%.3f (%s) %s" % \
                    (node.name, node.dist, str(node.data["logl"]), notes)

        else:
            labels[node.name] = "[%s]\n%.3f (*) %s" % \
                (node.name, node.dist, notes)
        
        if "params" in node.data:
            fracs = map(stats.mean, zip(* node.data["fracs"]))
            mean = sum(util.vmul(util.cget(node.data["params"], 0), fracs))
            sdev = sum(util.vmul(util.cget(node.data["params"], 1), fracs))
            
            mean *= baserate
            sdev *= baserate
            
            labels[node.name] += "\n%.3f %.3f" % (mean, sdev)
        
        if node in events:
            labels[node.name] += " %s" % events[node]
        
    if "logl" in tree.data:
        debug("logl: %f" % tree.data["logl"])
    treelib.drawTree(tree, minlen=20, labels=labels, spacing=4, 
                        labelOffset=-3, out=out)



class SindirError (Exception):
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    def __str__(self):
        return str(self.msg)



#-------------------------------------------------------------------------------
# SINDIR input/output
#-------------------------------------------------------------------------------

def writeParams(filename, params):
    """Write SINDIR model parameters to a file"""
    
    out = file(filename, "w")
    
    keys = util.sort(params.keys())
    
    for key in keys:
        values = params[key]
        print >>out, "%s\t%s" % (str(key), "\t".join(map(str,values)))


def readParams(filename):
    """Read SINDIR model parameters to a file"""
    
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


#-------------------------------------------------------------------------------
# Branch length fitting
#-------------------------------------------------------------------------------

def neighborjoin(distmat, genes):
    tree = treelib.Tree()
    leaves = {}
    dists = util.Dict(2, None)
    restdists = {}
    
    
    # initialize distances
    for i in range(len(genes)):
        r = 0
        for j in range(len(genes)):
            dists[genes[i]][genes[j]] = distmat[i][j]
            r += distmat[i][j]
        restdists[genes[i]] = r / (len(genes) - 2)
        
    # initialize leaves
    for gene in genes:
        tree.add(treelib.TreeNode(gene))
        leaves[gene] = 1
    
    # join loop
    while len(leaves) > 2:       
        # search for closest genes
        low = 1e1000
        lowpair = (None, None)
        leaveslst = leaves.keys()

        for i in range(len(leaves)):
            for j in range(i+1, len(leaves)):
                gene1, gene2 = leaveslst[i], leaveslst[j]
                dist = dists[gene1][gene2] - restdists[gene1] - restdists[gene2]
                
                if dist < low:
                    low = dist
                    lowpair = (gene1, gene2)
        
        # join gene1 and gene2
        gene1, gene2 = lowpair
        parent = treelib.TreeNode(tree.newName())
        tree.addChild(parent, tree.nodes[gene1])
        tree.addChild(parent, tree.nodes[gene2])
        
        # set distances
        tree.nodes[gene1].dist = (dists[gene1][gene2] + restdists[gene1] - 
                                  restdists[gene2]) / 2.0
        tree.nodes[gene2].dist = dists[gene1][gene2] - tree.nodes[gene1].dist
        
        # gene1 and gene2 are no longer leaves
        del leaves[gene1]
        del leaves[gene2]
        
        gene3 = parent.name
        r = 0
        for gene in leaves:
            dists[gene3][gene] = (dists[gene1][gene] + dists[gene2][gene] -
                                  dists[gene1][gene2]) / 2.0
            dists[gene][gene3] = dists[gene3][gene]
            r += distmat[i][j]
        leaves[gene3] = 1
        
        if len(leaves) > 2:
            restdists[gene3] = r / (len(leaves) - 2)
    
    # join the last two genes into a tribranch
    gene1, gene2 = leaves.keys()
    if type(gene1) == str:
        gene1, gene2 = gene2, gene1
    tree.addChild(tree.nodes[gene1], tree.nodes[gene2])
    tree.nodes[gene2].dist = dists[gene1][gene2]
    tree.root = tree.nodes[gene1]
    
    return tree



def findSplits(network, leaves):
    # find vertice and edge visit history
    start = network.keys()[0]

    openset = [start]
    closedset = {}
    
    vhistory = []
    ehistory = []
    elookup = util.Dict(1, [])
    
    
    while len(openset) > 0:
        vertex = openset.pop()
        
        vhistory.append(vertex)
        
        if len(vhistory) > 1:
            edge = tuple(util.sort(vhistory[-2:]))        
            ehistory.append(edge)
            elookup[edge].append(len(ehistory) - 1)
        
        # skip closed vertices
        if vertex in closedset:
            continue
        
        for v in network[vertex].keys():
            if v not in closedset:
                openset.append(vertex)            
                openset.append(v)
        

        # close new vertex
        closedset[vertex] = 1
    
    
    # use histories to define half each split
    splits = {}
    for edge in elookup:
        set1 = {}
        
        start, end = elookup[edge]
        for i in range(start+1, end+1):
            if vhistory[i] in leaves:
                set1[vhistory[i]] = 1
        
        # fill in other half of splits using complement
        set2 = {}
        for v in leaves:
            if v not in set1:
                set2[v] = 1
        
        if edge[0] == vhistory[start]:
            splits[edge] = [set2, set1]
        else:
            splits[edge] = [set1, set2]
        
    
    return splits



def setTreeDistances(conf, tree, distmat, genes):
    if not treelib.isRooted(tree):
        tree.addChild(treelib.TreeNode(tree.newName()), tree.root)
        tree.root = tree.root.parent
    network = treelib.tree2graph(tree)
        
    # create pairwise dist array
    dists = []
    for i in xrange(len(genes)):
        for j in xrange(i+1, len(genes)):
            dists.append(distmat[i][j])
    
    # find how edges split vertices
    splits = findSplits(network, util.makeset(genes))
    edges = splits.keys()
    
    # create topology matrix
    topmat = matrix.makeMatrix(len(dists), len(edges))
    
    vlookup = util.list2lookup(genes)
    n = len(genes)
    for e in xrange(len(edges)):
        set1, set2 = splits[edges[e]]
        for gene1 in set1:
            for gene2 in set2:
                i, j = util.sort([vlookup[gene1], vlookup[gene2]])
                index = i*n-i*(i+1)/2+j-i-1
                topmat[index][e] = 1
    
        
    A = scipy.array(topmat)
    d = scipy.array(dists)
    b,resids,rank,singlars = scipy.linalg.lstsq(A, d)
    
    if conf["debug"]:
        error = scipy.matrixmultiply(A, b) - scipy.transpose([d])
        debug("distance error:", max(max(abs(error))))
    
    for i in xrange(len(edges)):
        gene1, gene2 = edges[i]
        if tree.nodes[gene2].parent == tree.nodes[gene1]:
            gene1, gene2 = gene2, gene1
        if len(b.shape) == 2:
            tree.nodes[gene1].dist = float(b[i][0])
        else:
            tree.nodes[gene1].dist = float(b[i])
    

    
    if len(tree.root.children) == 1:
        tree.root = tree.root.children[0]
        tree.remove(tree.root.parent)
        tree.root.parent = None




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


def fitParams(lengths, baserates, gene2species, fit=False):
    ntrees = len(lengths.values()[0])
    
    params = {}
    
    dist = util.distrib(baserates, size=.01)
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
            dist = util.distrib(lens, size=.01)
            param, resid = stats.fitCurve(dist[0], dist[1], stats.normalPdf, [0,1])
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


def mleBaserates(lengths, params, baserateparam):
    lenmat = zip(* lengths.values())
    keys = map(lambda x: x.name, lengths.keys())
    means, sdevs = zip(* util.sublist(params, keys))
    baserates = []
    for i in xrange(len(lenmat)):
        baserates.append(mleBaserate(lenmat[i], means, sdevs, baserateparam))
    return baserates


def learnModel(trees, stree, gene2species, statsprefix=""):
    util.tic("learn model")

    util.tic("find branch length distributions")
    trees2, lengths = phyloutil.findBranchDistrib(trees, stree, gene2species,
                                                  False)
    debug("Total trees matching species topology: %d out f %d" % 
          (len(trees2), len(trees)))
    util.toc()
    
    params = {}
    
    totlens = map(sum, zip(* lengths.values()))
    
    if statsprefix != "":
        out = file(statsprefix + ".lens", "w")
        
        for node, lens in lengths.items():
            if len(lens) == 0 or max(lens) == min(lens):
                continue
        
            out.write(str(node.name))
            
            for length in lens:
                out.write("\t%f" % length)
            out.write("\n")
        
        out.close()
    
    
    
    util.tic("fitting params")
    for (node, lens), totlen in zip(lengths.items(), totlens):
        if len(lens) == 0 or max(lens) == min(lens):
            continue
        
        util.tic("fitting params for " + str(node.name))
        
        #ndivs = int(max(lens) / .001)
        #dist = util.distrib(lens, size=.001)
        #param, resid = stats.fitCurve(dist[0], dist[1], stats.normalPdf, [1,1])
        #param, resid = stats.fitCurve(dist[0], dist[1], stats.gammaDistrib, [1,1])
        
        param = fitNormal(util.vdivs(lens, totlen))
        
        params[node.name] = param
        util.toc()
    util.toc()

    # calc distribution of total tree length
    lens = map(lambda x: sum(y.dist for y in x.nodes.values()), trees2)
    lens = filter(lambda x: x < 20, lens)
    mu = stats.mean(lens)
    lens = filter(lambda x: x < 2*mu, lens)
    mu = stats.mean(lens)
    sigma2 = stats.variance(lens)
    params["baserate"] = [mu*mu/sigma2, mu/sigma2]
    params[stree.root.name] = [0, 1]
    #util.writeVector("treelens", lens)
    
    util.toc()
    
    return params
    

def learnModel2(lengths, gene2species, niters=10, fit=True):
    lenmat = zip(* lengths.values())
    keys = map(lambda x: x.name, lengths.keys())

    # init base rates
    baserates = map(sum, lenmat)
    
    baseratesList = [baserates]
    paramsList = []
    
    # fit baserate distribution
    dist = util.distrib(baserates, size=.2)
    baserateparam, resid = stats.fitCurve(dist[0], dist[1], stats.gammaPdf, [1,1])
    
    
    # do EM
    for i in range(niters):
        params = fitParams(lengths, baserates, gene2species, fit=fit)
        means, sdevs = zip(* util.sublist(params, keys))
        
        paramsList.append(params)
        
        baserates = []
        for i in xrange(len(lenmat)):
            baserates.append(mleBaserate(lenmat[i], means, sdevs, baserateparam))
        
        
        #factor = stats.mean(util.vdiv(baseratesList[0], baserates))
        #baserates = [x*factor for x in baserates]
        
        baseratesList.append(baserates)
        
        # calc likelihood
        util.log(dataLikelihood(lenmat, baserates, means, sdevs, baserateparam))
        
    
    return paramsList, baseratesList



#-------------------------------------------------------------------------------
# Likelihood calculation
#-------------------------------------------------------------------------------


def mleBaserate2(lens, means, sdevs, baserateparam):
    vars = util.vmul(sdevs, sdevs)
    return sum(util.vdiv(util.vmul(lens, lens), vars)) / \
           sum(util.vdiv(util.vmul(means, lens), vars))


def mleBaserate(lens, means, sdevs, baserateparam):
    [alpha, beta] = baserateparam
    
    # protect against zero
    ind = util.findgt(.0001, sdevs)
    lens = util.sublist(lens, ind)
    means = util.sublist(means, ind)
    sdevs = util.sublist(sdevs, ind)
    
    a = (1 - alpha) / beta
    b = sum(means[i] * lens[i] / sdevs[i]**2
            for i in range(len(lens))) / beta
    c = - sum(lens[i] ** 2 / sdevs[i] ** 2
              for i in range(len(lens))) / beta
    
    #print filter(lambda x: x>0, stats.solveCubic(a, b, c))
    return max(stats.solveCubic(a, b, c))


def log(x):
    """Safe logarithm function"""
    
    if x <= 0:
        return -1e100
    else:
        return math.log(x)



def branchLikelihood(dist, fracs, params):
    totmean = 0.0
    totvar = 0.0
    
    for frac, (mean, sigma) in zip(fracs, params):
        totmean += frac * mean
        totvar  += frac * sigma * sigma
    
    # don't let variance get too low
    if abs(totvar) < .000001:
        totvar = .000001
    
    return log(stats.normalPdf(dist, [totmean, math.sqrt(totvar)]))


def setMidpoints(node, events, recon, midpoints, kvars):
    this = util.Closure(i = 0)

    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup":
            if recon[node] == recon[node.parent]:
                # if im the same species branch as my parent 
                # then he is my last midpoint
                lastpoint = midpoints[node.parent]
            else:
                # im the first on this branch so the last midpoint is zero
                lastpoint = 0.0
            
            # pick a midpoint uniformly after the last one
            midpoints[node] = lastpoint + \
                        (kvars[this.i] * (1 - lastpoint))
            this.i += 1
        else:
            # genes or speciations reconcile exactly to the end of the branch
            midpoints[node] = 1.0
    
        # recurse within dup-only subtree
        if events[node] == "dup":
            node.recurse(walk)
    
    walk(node)


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
    this = util.Closure(
        logl=0.0
    )
    
    def walk(node):
        # branch that recon to the species root are "free"
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


def subtreeLikelihood(conf, root, recon, events, stree, params, baserate):
    midpoints = {}
    extraBranches = getExtraBranches(root, recon, events, stree)
    
    this = util.Closure(ncalls=0)    
    
    def integrate(node, nvars):
        kvars = []
        
        def func(k):
            kvars.append(k)
            
            if len(kvars) == nvars:
                this.ncalls += 1
                setMidpoints(child, events, recon, midpoints, kvars)
                kvars.pop()
                return math.e** \
                 calcSubtreeLikelihood(node, recon, events, stree, params, 
                                       midpoints, extraBranches, baserate)
            else:
                return scipy.integrate.quad(func, .05, .95, epsrel=conf["accuracy"])[0]
        return scipy.integrate.quad(func, .05, .95, epsrel=conf["accuracy"])
    
    
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
            if conf["debug"]:
                debug("integration error:", (log(err + val) - log(val)))
                debug("nvars:", nvars)
            clogl = log(val)
            child.data["logl"] = clogl
            logl += clogl
    
    if conf["debug"]:
        debug("logl calls:", this.ncalls)
    
    return logl


def subtreeLikelihood2(conf, root, recon, events, stree, params, baserate):
    midpoints = {}
    extraBranches = getExtraBranches(root, recon, events, stree)

    
    def integrate(node, nvars):
        kvars = []
        
        def func(k):
            kvars.append(k)
            
            if len(kvars) == nvars:
                setMidpoints(child, events, recon, midpoints, kvars)
                kvars.pop()
                return math.e** \
                 calcSubtreeLikelihood(node, recon, events, stree, params, 
                                       midpoints, extraBranches, baserate)
            else:
                return scipy.integrate.romberg(func, .01, .99, conf["accuracy"])
        return scipy.integrate.romberg(func, .01, .99, conf["accuracy"])
    
    
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
            clogl = log(integrate(child, nvars))
            child.data["logl"] = clogl
            logl += clogl
    
    return logl






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


def subtreeLikelihoodRandom(conf, root, recon, events, stree, params):
    midpoints = {}
    extraBranches = getExtraBranches(root, recon, events, stree)
    
    nsamples = 10
    
    # process each child of subtree root
    logl = 0.0
    midpoints[root] = 1.0
    for child in root.children:
        prob = 0.0
        for sample in range(nsamples):
            setMidpointsRandom(child, events, recon, midpoints)
            prob += math.e** \
             calcSubtreeLikelihood(child, recon, events, stree, params, 
                          midpoints, extraBranches)
            
            # integration is only needed if child is dup
            if events[child] != "dup":
                break
        clogl = log(prob / float(sample+1))
        child.data["logl"] = clogl
        logl += clogl
    
    return logl


def getBaserate(tree, stree, params, recon=None, gene2species=None):
    if recon == None:
        assert gene2species != None
        recon = phyloutil.reconcile(tree, stree, gene2species)
    events = phyloutil.labelEvents(tree, recon)
    
    extraBranches = getExtraBranches(tree.root, recon, events, stree)
    baserateparam = params["baserate"]
    
    
    
    def walk(node):
        dist, fracs, snodes = reconBranch(node, recon, events, stree, 
                params, midpoints, extraBranches)
        params2 = [params[x] for x in snodes]
        
        #print dist, fracs, snodes
        
        # do not try to use extra branches for estimating base rate
        if "extra" not in node.data and \
           recon[node] != stree.root:
            lens.append(dist)
            means.append(util.vdot(fracs, util.cget(params2, 0)))
            s = util.cget(params2, 1)
            v = util.vmul(s, s)
            sdevs.append(math.sqrt(util.vdot(fracs, v)))
        
        node.recurse(walk)
    
    totallens = []
    for i in range(10):
        lens = []
        means = []
        sdevs = []
        
        midpoints = {}
        setMidpointsRandom(tree.root, events, recon, midpoints, wholeTree=True)
        for child in tree.root.children:
            walk(child)
        
        totallens.append(mleBaserate(lens, means, sdevs, baserateparam))
    
    return stats.mean(totallens)


def getBaserate2(tree, stree, params, recon=None, gene2species=None):
    if recon == None:
        assert gene2species != None
        recon = phyloutil.reconcile(tree, stree, gene2species)
    
    lens = []
    means = []
    sdevs = []
    
    # TODO: this is probably wrong!  Need to do 
    for node in tree.nodes.values():
        if recon[node] == stree.root:
            continue
        lens.append(node.dist)
        means.append(params[recon[node].name][0])
        sdevs.append(params[recon[node].name][1])
    
    baserateparam = params["baserate"]
    
    totallen = mleBaserate(lens, means, sdevs, baserateparam)
    
    return totallen


def subtreeLikelihood3(conf, root, recon, events, stree, params, baserate):
    extraBranches = getExtraBranches(root, recon, events, stree)

    this = util.Closure(
        logl=0.0
    )
    
    depths = {root.parent: 0}
    marks = {root.parent: 1}
    
    sroot = recon[root.parent]
    
    # process each child of subtree root
    def walk(node, depth):        
        # save depth of node
        depths[node] = node.dist + depth
        
        
        if events[node] == "dup":
            # recurse within dup-only subtree
            node.recurse(walk, depth + node.dist)
        else:
            # we are at subtree leaf
            # compute likelihood of path from leaf to root

            # branch that recon to the species root are "free"
            if recon[node] == stree.root:
                return
            
            # figure out species branches that we cross
            # get total mean and variance of this path            
            mu = 0
            sigma = 0            
            snode = recon[node]
            while snode != sroot and snode != stree.root:
                mu += params[snode.name][0]
                sigma += params[snode.name][1]
                snode = snode.parent
            
            # find out how much of our path is not conditioned
            ptr = node
            while ptr not in marks:
                marks[ptr] = 1
                ptr = ptr.parent
            assert node != ptr
            
            assert abs(sigma) > .00000001
            
            if depths[ptr] == 0.0:
                denom = 1.0
            else:
                denom = 1 - stats.normalCdf(depths[ptr], [mu, sigma])
            
            print node.name, mu, sigma, depths[node], depths[ptr]
            
            if denom < .0000001:
                this.logl = -util.INF
            else:
                this.logl += log(stats.normalPdf(depths[node], [mu, sigma]) /
                                 denom)
            
            node.data["params"] = [[mu, sigma]]
            node.data["fracs"] = [[1]]
            
    walk(root, 0)
    
    # debug saving
    root.data["logl"] = this.logl
    
    return this.logl



def treeLogLikelihood(conf, tree, stree, gene2species, params):
    # derive relative branch lengths
    tree.clearData("logl", "extra", "fracs", "params", "unfold")
    recon = phyloutil.reconcile(tree, stree, gene2species)
    events = phyloutil.labelEvents(tree, recon)

    # determine if top branch unfolds
    if recon[tree.root] ==  stree.root and \
       events[tree.root] == "dup":
        for child in tree.root.children:
            if recon[child] != stree.root:
                child.data["unfold"] = True
    
    if conf["debug"]:
        print >>DEBUG, "\ntreelen: ", sum(x.dist for x in tree.nodes.values())
        util.tic("find baserate")
    baserate = getBaserate(tree, stree, params, recon=recon)
    if conf["debug"]:
        print >>DEBUG, "baserate:", baserate, log(stats.gammaPdf(baserate, [9.8, 8.4]))
        util.toc()    

    if conf["debug"]:
        util.tic("find logl")
    
    # top branch is "free"
    params[stree.root.name] = [0,0]
    this = util.Closure(logl=0.0)
    
    # recurse through indep sub-trees
    def walk(node):
        if events[node] == "spec" or \
           node == tree.root:
            #this.logl += subtreeLikelihood(conf, node, recon, events, 
            #                               stree, params, baserate)
            
            for child in node.children:
                this.logl += subtreeLikelihood3(conf, child, recon, events, 
                                                stree, params, baserate)
            
        node.recurse(walk)
    walk(tree.root)
    
    if conf["debug"]:
        util.toc()
    
    # calc probability of rare events
    for node, event in events.items():
        if event == "dup":
            this.logl += log(conf["dupprob"])
        elif event == "spec":
            this.logl += log(conf["specprob"])
    
    nloss = len(phyloutil.findLoss(tree, stree, recon))
    this.logl += nloss * log(conf["lossprob"])
    #this.logl += (len(tree.nodes) - 1) * log(1 - conf["lossprob"])
    
    
    # debugging information
    if conf["debug"]:
        debug("\n\n")
        drawTreeLogl(tree, events=events)
    
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    return this.logl








#-------------------------------------------------------------------------------
# Tree search
#-------------------------------------------------------------------------------

def proposeNni(tree, node1, node2, change=0):
    """Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      uncle    node1
               /  \
         child[0]  child[1]
    
    """
    
    # ensure node2 is parent of node1
    if node1.parent != node2:
        node1, node2 = node2, node1
    assert node1.parent == node2
    
    # find uncle
    uncle = 0 
    if node2.children[uncle] == node1:
        uncle = 1
    
    # swap parent pointers
    node1.children[change].parent = node2
    node2.children[uncle].parent = node1
    
    # swap child pointers
    node2.children[uncle], node1.children[change] = \
        node1.children[change], node2.children[uncle]


def proposeTree(tree):
    tree2 = tree.copy()
    
    # find edges for NNI
    nodes = tree2.nodes.values()
    nodes = filter(lambda x: not x.isLeaf() and 
                             x != tree2.root, nodes)
    edges = [(node, node.parent) for node in nodes]
    edge = edges[int(random.random() * len(edges))]
    
    proposeNni(tree2, edge[0], edge[1], int(round(random.random())))
    return tree2


def searchMCMC(conf, distmat, labels, stree, gene2species, params,
               initTree=None):
    # init with NJ
    
    if initTree != None:
        tree = initTree
    else:
        tree = neighborjoin(distmat, labels)
        tree = phyloutil.reconRoot(tree, stree, gene2species)
        setTreeDistances(conf, tree, distmat, labels)
    
    
    # init likelihood score
    top = treeLogLikelihood(conf, tree, stree, gene2species, params)
    toptree = tree.copy()
    thash = phyloutil.hashTree(tree, lambda x: x)
    
    visited = {thash: (top, tree.copy())}
    
    
    debug("\n=======================================")
    debug("i:", 0, "logl:", top, "top:", top)
    drawTreeLogl(tree)
    debug()
    debug()


    # tree search
    lastl = top
    for i in xrange(1, conf["iters"]):
        tree2 = proposeTree(tree)
        tree2 = proposeTree(tree2)
        
        # just for debug
        recon = phyloutil.reconcile(tree2, stree, gene2species)
        events = phyloutil.labelEvents(tree2, recon)
        
        #tree2 = proposeTree(tree2)
        #tree2 = phyloutil.reconRoot(tree2, stree, gene2species)
        
        thash = phyloutil.hashTree(tree2, lambda x: x)
        if thash in visited:
            logl, tree2 = visited[thash]
        else:
            setTreeDistances(conf, tree2, distmat, labels)      
            logl = treeLogLikelihood(conf, tree2, stree, gene2species, params)
            
        
        # store logl in visited
        if thash not in visited or logl > visited[thash][0]:
            visited[thash] = (logl, tree2.copy())
        
        
        # best yet tree
        if logl > top:
            debug("\n=======================================")
            debug("i:", i, "logl:", logl, "top:", top)
            drawTreeLogl(tree2, events=events)
            debug()
            debug()
        
            top = logl
            toptree = tree2.copy()


        # accept/reject
        if logl > lastl:
            # accept new tree
            tree = tree2
            lastl = logl
            if logl - top > log(.90):
                print >>DEBUG, "#",
            elif logl - top > log(.50):
                print >>DEBUG, "+",
            else:
                print >>DEBUG, "^",
            DEBUG.flush()
        else:
            # accept with a chance
            if logl - lastl > log(random.random()):
                print >>DEBUG, "v",
                tree = tree2
                lastl = logl
            else:
                print >>DEBUG, "_",
            DEBUG.flush()
        
    
    debug("\n\nmost visited trees out of %d: " % len(visited))
    visited = util.mapdict(visited, valfunc=lambda x: "%.4f" % x[0])
    util.printDictByValues(visited, num=30, spacing=4, 
                           compare=lambda a,b: cmp(float(b),float(a)), 
                           out=DEBUG)
    debug()

    return toptree, top





def searchGreedy(conf, distmat, labels, stree, gene2species, params):
    totalgenes = len(labels)
    ngenes = 2
    
    # create initial 2 gene tree (labels[0], labels[1])
    tree = treelib.Tree()
    tree.makeRoot()
    tree.addChild(tree.root, treelib.TreeNode(labels[0]))
    tree.addChild(tree.root, treelib.TreeNode(labels[1]))
        
    
    for ngenes in xrange(2, totalgenes):
        debug("adding", labels[ngenes])
        
        toplogl = -1e100
        toptree = None
        
        distmat2 = matrix.submatrix(distmat, range(ngenes+1), range(ngenes+1))
        labels2  = labels[:ngenes+1]
        
        
        # place new gene on every branch
        for name in tree.nodes:
            tree2 = tree.copy()
            node = tree2.nodes[name]

            if node == tree2.root:
                newnode = treelib.TreeNode(tree2.newName())
                tree2.addChild(newnode, tree2.root)
                tree2.root = newnode
                tree2.addChild(newnode, treelib.TreeNode(labels[ngenes]))
            else:
                parent = node.parent
                tree2.remove(node)
                newnode = treelib.TreeNode(tree2.newName())
                tree2.addChild(parent, newnode)
                tree2.addChild(newnode, node)
                tree2.addChild(newnode, treelib.TreeNode(labels[ngenes]))
            
            #tree2 = phyloutil.reconRoot(tree2, stree, gene2species)
            setTreeDistances(conf, tree2, distmat2, labels2)
            logl = treeLogLikelihood(conf, tree2, stree, gene2species, params)

            if logl >= toplogl:
                toplogl = logl
                toptree = tree2
        tree = toptree
        
        visited = {}
        tree, logl = searchExhaustive(conf, distmat2, labels2, 
                                      tree, stree, gene2species, params,
                                      visited=visited)
        if logl >= toplogl:
            toplogl = logl
            toptree = tree
        tree = toptree
        
        
        debug()
    
    return tree, toplogl



def searchExhaustive(conf, distmat, labels, tree, stree, gene2species, params,
                     depth=2, visited=None):
    if visited == None:
        visited = {}
    
    # find initial logl
    thash = phyloutil.hashTree(tree)
    setTreeDistances(conf, tree, distmat, labels)
    toplogl = treeLogLikelihood(conf, tree, stree, 
                             gene2species, params)
    visited[thash] = toplogl
    toptree = tree
    
    
    # try all NNI
    # find edges for NNI
    nodes = tree.nodes.values()
    nodes = filter(lambda x: not x.isLeaf() and 
                             x != tree.root, nodes)
    edges = [(node, node.parent) for node in nodes]

    for edge in edges:
        for change in (0,1):
            proposeNni(tree, edge[0], edge[1], change)
            
            # calc logl
            if depth > 1:
                tree, logl = searchExhaustive(conf, distmat, labels, 
                                              tree, stree, gene2species, params,
                                              depth=depth-1, visited=visited)
            else:
                thash = phyloutil.hashTree(tree)
                if thash in visited:
                    logl = visited[thash]
                else:
                    setTreeDistances(conf, tree, distmat, labels)
                    logl = treeLogLikelihood(conf, tree, stree, 
                                             gene2species, params)
                    visited[thash] = logl
            
            # save max logl
            if logl > toplogl:
                toplogl = logl
                toptree = tree.copy()

            # switch branch back
            proposeNni(tree, edge[0], edge[1], change)
    
    # debug
    debug("\n\nmost visited trees out of %d: " % len(visited))
    visited = util.mapdict(visited, valfunc=lambda x: "%.4f" % x)
    util.printDictByValues(visited, num=40, spacing=4, 
                           compare=lambda a,b: cmp(float(b),float(a)), 
                           out=DEBUG)
    
    return toptree, toplogl


#-------------------------------------------------------------------------------
# Main SINDIR algorithm function
#-------------------------------------------------------------------------------

def sindir(conf, distmat, labels, stree, gene2species, params):
    """Main function for the SINDIR algorithm"""

    trees = []
    logls = []
    tree = None
    
    util.tic("SINDIR")
    
    for search in conf["search"]:

        if search == "greedy":
            tree, logl = searchGreedy(conf, distmat, labels, stree, 
                                   gene2species, params)        
            
        elif search == "mcmc":
            tree, logl = searchMCMC(conf, distmat, labels, stree, 
                                    gene2species, params, initTree=tree)
        elif search == "exhaustive":
            if tree == None:
                tree = neighborjoin(distmat, labels)
                tree = phyloutil.reconRoot(tree, stree, gene2species)
                setTreeDistances(conf, tree, distmat, labels)
                
                drawTreeLogl(tree)
        
            tree, logl = searchExhaustive(conf, distmat, labels, tree, stree, 
                                    gene2species, params)
        elif search == "none":
            break
        else:
            raise SindirError("unknown search '%s'" % search)
        
        trees.append(tree)
        logls.append(logl)
        
    
    # eval the user given trees
    for treefile in conf["tree"]:
        tree2 = treelib.Tree()
        tree2.readNewick(treefile)        
        setTreeDistances(conf, tree2, distmat, labels)
        logl2 = treeLogLikelihood(conf, tree2, stree, gene2species, params)
    
        trees.append(tree2)
        logls.append(logl2)
    
        debug("\nuser given tree:")
        debug("logl:", logl2)
        drawTreeLogl(tree2)
    
    util.toc()
    
    if len(trees) == 0:
        raise SindirError("No search or tree topologies given")
    
    return trees[util.argmax(logls)], max(logls)
