#
# SPIDIR library
#
# note: SPIDIR was a codename, it may still be present in the code.
#

# python libs
import math, StringIO, copy, random, sys


# rasmus libs
from rasmus import bionj
from rasmus import fasta
from rasmus import matrix
from rasmus import phyloutil
from rasmus import stats
from rasmus import tablelib
from rasmus import treelib
from rasmus import util
from rasmus import phylip
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


def treeDistrib2table(lengths):
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
        low = util.INF
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

def makeVector(array):
    if len(array.shape) == 2:
        if array.shape[0] == 1:
            return array[0]
        else:
            return scipy.transpose(array)[0]
    else:
        return array

def setTreeDistances(conf, tree, distmat, genes):
    if isDebug(DEBUG_MED):
        util.tic("fit branch lengths")
    
    if treelib.isRooted(tree):
        rootedge = sorted([x.name for x in tree.root.children])
        treelib.unroot(tree, newCopy=False)
    else:
        rootedge = None
    
            
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
    
    # force non-negative branch lengths
    b = [max(float(x), 0) for x in makeVector(b)]
    #b = [float(x) for x in makeVector(b)]
    
    d2 = makeVector(scipy.dot(A, b))
    resids = (d2 - d).tolist()
    d = d.tolist()
    
    
    # recreate rooting branches
    if rootedge != None:
        # restore original rooting
        if tree.nodes[rootedge[0]].parent == tree.nodes[rootedge[1]]:
            treelib.reroot(tree, rootedge[0], newCopy=False)
        else:
            treelib.reroot(tree, rootedge[1], newCopy=False)
    
        # find root edge in edges
        for i in xrange(len(edges)):
            if sorted(edges[i]) == rootedge:
                break
                
        edges[i] = [rootedge[0], tree.root.name]
        edges.append([rootedge[1], tree.root.name])
        b[i] /= 2.0
        b.append(b[i])
        resids[i] /= 2.0
        resids.append(resids[i])
        d[i] /= 2.0
        d.append(d[i])
        
        for row in topmat:
            row.append(row[i])

    
    for i in xrange(len(edges)):
        gene1, gene2 = edges[i]
        if tree.nodes[gene2].parent == tree.nodes[gene1]:
            gene1, gene2 = gene2, gene1
        tree.nodes[gene1].dist = b[i]
    
    #for node in tree.nodes.values():
    #    assert node.dist >= 0
    
    
    # catch unusual case that may occur in greedy search
    if sum(x.dist for x in tree.nodes.values()) == 0:
        for node in tree.nodes.values():
            node.dist = .01
    
    tree.data["error"] = math.sqrt(scipy.dot(resids, resids)) / \
                                   sum(x.dist for x in tree.nodes.values())
    
    #util.writeVector("resids", resids)
    #util.writeVector("dists", d)
    #util.writeVector("dists2", d2)
    
    setBranchError(conf, tree, resids, d, edges, topmat)
        
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
        
    
def setBranchError2(conf, tree, pathErrors, edges, topmat):
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
                
                tree.nodes[gene1].data["error"] += pathErrors[i]**2
                npaths[gene1] += 1
                
    
    # normalize by number of paths through edge
    for node in tree.nodes.itervalues():
        if node.name in npaths:
            node.data["error"] /= npaths[node.name]
        


def getSplit(tree):
    splits = findSplits(treelib.tree2graph(tree), tree.leafNames())
    splits2 = {}
    
    for edge, sets in splits.iteritems():
        # skip external edges
        if len(sets[0]) == 1 or len(sets[1]) == 1:
            continue
        
        s = tuple(sorted([tuple(sorted(i.keys())) for i in sets]))
        splits2[edge] = s
    
    # if tree is rooted, remove duplicate edge
    if treelib.isRooted(tree):
        edge1 = tuple(sorted([tree.root.name, tree.root.children[0].name]))
        edge2 = tuple(sorted([tree.root.name, tree.root.children[1].name]))
        if edge1 > edge2:
            edge1, edge2 = edge2, edge1
        if edge1 in splits2 and edge2 in splits2:
            del splits2[edge1]
    
    return splits2
              

def robinsonFouldsError(tree1, tree2):
    splits1 = getSplit(tree1)
    splits2 = getSplit(tree2)

    overlap = set(splits1.values()) & set(splits2.values())
    
    #assert len(splits1) == len(splits2)

    return 1 - (len(overlap) / float(max(len(splits1), len(splits2))))


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
    debug("Total trees matching species topology: %d out of %d" % 
          (len(trees2), len(trees)))
    util.toc()
    
    params = {}
    
    totlens = map(sum, zip(* lengths.values()))
    
    # print output stats
    if statsprefix != "":
        writeTreeDistrib(file(statsprefix + ".lens", "w"), lengths)
        rates = treeDistrib2table(lengths)
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
    
    ind = range(len(means))
    ind.sort(lambda a, b: cmp(means[b], means[a]))
    #ind = ind[:len(ind) / 4 + 1]
    ind = ind[:max(4, len(means)-2)]
    means = util.mget(means, ind)
    sdevs = util.mget(sdevs, ind)
    lens = util.mget(lens, ind)
    
    
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


def mleBaserateG(lens, alphas, betas, baserateparam):
    [alpha, beta] = baserateparam
    
    # protect against zero
    #ind = util.findgt(.0001, sdevs)
    #lens = util.sublist(lens, ind)
    #alphas = util.sublist(means, ind)
    #betas = util.sublist(sdevs, ind)
    
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
        recon = phyloutil.reconcile(tree, stree, gene2species)
    events = phyloutil.labelEvents(tree, recon)
    
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
        
        if event == "dup":
            logl += log(conf["dupprob"])
        
    #nloss = len(phyloutil.findLoss(tree, stree, recon))
    #logl += nloss * log(conf["lossprob"])
    
    return logl


def subtreeLikelihood_est(conf, tree, root, recon, events, stree, params, baserate):
    this = util.Closure(
        logl=0.0,
        isExtra = False
    )
    
    extraBranches = getExtraBranches(root, recon, events, stree)
    depths = {root.parent: 0}
    marks = {root.parent: 1}
    sroot = recon[root.parent]
    
    
    # process each child of subtree root
    def walk(node, extra):        
        # save depth of node
        if recon[node] != recon[tree.root]:  #stree.root:
            depths[node] = node.dist + depths[node.parent]
        else:
            # ignore branch length of free branches
            depths[node] = depths[node.parent]
        
        # remember if extra node is in path
        if "extra" in node.data:
            extra = node
        
        
        if events[node] == "dup":
            # recurse within dup-only subtree
            node.recurse(walk, extra)
        else:
            # we are at subtree leaf
            # compute likelihood of path from leaf to root
            
            # figure out species branches that we cross
            # get total mean and variance of this path
            mu = 0
            sigma2 = 0
            snode = recon[node]
            
            # branch is free if we do not cross any more species
            if snode == sroot:
                return
            
            # sum means and variances along path
            while snode != sroot and snode != stree.root:
                mu += params[snode.name][0]
                sigma2 += params[snode.name][1]**2
                snode = snode.parent
            assert abs(sigma2) > .00000001, "sigma too small"
            sigma = math.sqrt(sigma2)
            
            
            # find out how much of our path is conditioned
            ptr = node
            while ptr not in marks:
                marks[ptr] = 1
                ptr = ptr.parent
            assert node != ptr
            condDist = depths[ptr]
            
            if condDist == 0.0:
                # if no distance to condition on denominator is 1.0
                logdenom = log(1.0)
            else:
                logdenom = log(1 - stats.normalCdf(condDist/baserate, [mu, sigma]))
            
            # determine dist of total path
            dist = max(depths[node], condDist)
            #dist = depths[node]
            
            
            # handle extra branches
            if extra != None:
                if "unfold" in extra.data:
                    dist += extra.dist
                
                # determine desired shrink
                target = min(mu, max(dist/baserate,0)) * baserate
                shrink = dist - target
                
                # determine how much shrink is allowed
                if "unfold" in extra.data:
                    extradist = max(2 * extra.dist, 0)
                else:
                    extradist = max(extra.dist, 0)
                shrink = min(shrink, extradist)
                
                if condDist == 0.0:
                    dist -= shrink
                else:
                    condDist -= shrink
            
            
            lognom = log(stats.normalPdf(dist/baserate, [mu, sigma]))
            
            
            if logdenom == -util.INF or \
               lognom   == util.INF:
                logl = -util.INF
                this.logl = -util.INF
            else:
                logl = lognom - logdenom
                this.logl += logl
            
                        
            """
            print "\t".join(["%10s" % str(node.name), 
                             "%.3f" % dist, 
                             "%.3f |" % condDist,
                             "%.3f" % (dist / baserate), 
                             "%.3f |" % (condDist / baserate), 
                             "%.3f" % mu, 
                             "%.3f |" % sigma,
                             "%.3f" % logdenom])
            """
            
            
            # debug saving
            node.data["logl"] = logl
            
            if this.logl > 1e10:
                debug(dist, condDist, baserate, mu, sigma, logl, logdenom)
                raise Exception("logl too high")
            
            node.data["params"] = [[mu, sigma]]
            node.data["fracs"] = [[1]]
            
    walk(root, None)    
    
    
    return this.logl




def branchLikelihoods(conf, tree, recon, events, stree, params, baserate):
    this = util.Closure(logl=0.0)

    # determine if top branch unfolds
    if recon[tree.root] ==  stree.root and \
       events[tree.root] == "dup":
        for child in tree.root.children:
            if recon[child] != stree.root:
                child.data["unfold"] = True    
    
    # determine if top branch needs to slide root
    if recon[tree.root] ==  stree.root and \
       len(tree.root.children) == 2 and \
       events[tree.root] == "spec" and \
       events[tree.root.children[0]] == "spec" and \
       events[tree.root.children[1]] == "spec":
        
        spath1 = 0
        snode1 = recon[tree.root.children[0]]
        while snode1 != stree.root:
            spath1 += params[snode1.name][0]
            snode1 = snode1.parent
        
        spath2 = 0
        snode2 = recon[tree.root.children[1]]
        while snode2 != stree.root:
            spath2 += params[snode2.name][0]
            snode2 = snode2.parent
        
        ratio = spath1 / float(spath1 + spath2)
        tot = tree.root.children[0].dist + tree.root.children[1].dist
        
        tree.root.children[0].dist = tot * ratio
        tree.root.children[1].dist = tot * (1 - ratio)
    
    
    # recurse through indep sub-trees
    def walk(node):
        if events[node] == "spec" or \
           node == tree.root:
            for child in node.children:
                this.logl += subtreeLikelihood_est(conf, tree, child, recon, events, 
                                                stree, params, baserate)
        node.recurse(walk)
    walk(tree.root)
    
    return this.logl


def treeLogLikelihood_est(conf, tree, stree, gene2species, params, baserate=None):
    # reconcile the gene tree
    # determine all events
    tree.clearData("logl", "extra", "fracs", "params", "unfold")
    recon = phyloutil.reconcile(tree, stree, gene2species)
    events = phyloutil.labelEvents(tree, recon)
    
    # determine baserate
    if baserate == None:
        baserate = getBaserate(tree, stree, params, recon=recon)
    
    

    
    # calc branch length likelihoods
    this = util.Closure(logl=0.0)
    this.logl = branchLikelihoods(conf, tree, recon, events, 
                                  stree, params, baserate)
    
    # calc baserate likelihood
    this.logl += log(stats.gammaPdf(baserate, params["baserate"]))
    
    # calc probability of rare events
    tree.data["eventlogl"] = rareEventsLikelihood(conf, tree, stree, recon, events)
    this.logl += tree.data["eventlogl"]
    
    # calc penality of error
    tree.data["errorlogl"] = tree.data["error"] * conf["errorcost"]
    this.logl += tree.data["errorlogl"]
    
    
    # debugging information
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    
    
    
    return this.logl

"""

def treeLogLikelihoodAllRoots(conf, tree, stree, gene2species, params, 
                              baserate=None):
    # find reconciliation that minimizes loss
    toplogl = -util.INF
    toproot = None
          
    # make an unrooted copy of gene tree
    tree = treelib.unroot(tree)
    
    # determine graph and possible roots
    mat = treelib.tree2graph(tree)
    newroots = util.sort(tree.nodes.keys())
    newroots.remove(tree.root.name)
    
    # try rooting on everything
    for root in newroots:
        
        tree2 = treelib.reroot(tree, root, mat)
        
        logl = treeLogLikelihood(conf, tree2, stree, gene2species, params, 
                                 baserate=baserate)
        
        # keep track of min loss
        if logl > toplogl:
            toplogl = logl
            toproot = root
    
    print toproot, toplogl
            
    # root tree by minroot
    toptree = treelib.reroot(tree, toproot)
    toplogl = treeLogLikelihood(conf, toptree, stree, gene2species, params, 
                                baserate=baserate)
    
    return toptree, toplogl

"""



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


def setMidpoints(node, events, recon, midpoints, kvars):
    this = util.Closure(i = 0)
    
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


def subtreeLikelihood(conf, root, recon, events, stree, params, baserate):
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
                
                if this.printing < 0:
                    debug("int:", kvars)
                    this.printing = 100
                else:
                    this.printing -= 1
                
                return math.e** \
                 calcSubtreeLikelihood(node, recon, events, stree, params, 
                                       midpoints, extraBranches, baserate)
            else:
                depnode = kdepend[korderrev[this.depth]]                
                if depnode == None:
                    lowk = 0
                else:
                    lowk = kvars[korder[depnode]] + .0001
                
                n = max(3, 10 / (2**this.depth))
                #n = 20
                
                ret = quad(func, lowk, .9999, n=n)[0] / \
                      (1.0 - lowk)
                this.depth -= 1
                return ret
        return quad(func, .0001, .9999, n=20)
    
    
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
                #debug("integration error:", (log(err + val) - log(val)))
                debug("int nvars:", nvars)
            clogl = log(val)
            child.data["logl"] = clogl
            logl += clogl
    
    if conf["debug"]:
        debug("int logl calls:", this.ncalls)
    
    return logl



def treeLogLikelihood(conf, tree, stree, gene2species, params, baserate=None):
    #conf["accuracy"] = .01
    conf.setdefault("bestlogl", -util.INF)

    estlogl = treeLogLikelihood_est(conf, tree, stree, gene2species, params, baserate=None)
    if "integrate" not in conf:
        return estlogl
    else:
        # skip trees that are estimated to be very bad
        if estlogl < conf["bestlogl"] - 50:
            debug("SKIP integration, est logl=%f" % estlogl)
            return estlogl
        conf["bestlogl"] = max(conf["bestlogl"], estlogl)
            
    #Search.printMCMC(conf, "est", tree, stree, gene2species, {})


    # debug info
    if isDebug(DEBUG_MED):
        util.tic("find logl")
    

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
                                           stree, params, baserate)
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
    
    print >> conf["intcmp"], "%f\t%f" % (this.logl, estlogl)
    conf["intcmp"].flush()
    
    
    if isDebug(DEBUG_MED):
        util.toc()
        debug("\n\n")
        drawTreeLogl(tree, events=events)
    
    return this.logl


### END TESTING INTEGRATION







#-------------------------------------------------------------------------------
# Main SPIDIR algorithm function
#-------------------------------------------------------------------------------

def spidir(conf, distmat, labels, stree, gene2species, params):
    """Main function for the SPIDIR algorithm"""
    
    setDebug(conf["debug"])

    
    if "integrate" in conf and "out" in conf:
        conf["intcmp"] = file(conf["out"] + ".int", "w")
        
    if "out" in conf:
        # create debug table
        conf["debugtab_file"] = file(conf["out"] + ".debug.tab", "w")
        
        debugtab = tablelib.Table(headers=["logl", "treelen", "baserate", 
            "error", "eventlogl", "tree"])
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
                tree = neighborjoin(distmat, labels)
                tree = phyloutil.reconRoot(tree, stree, gene2species)
            
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
        
        thash = phyloutil.hashTree(tree)
        if thash in visited:
            a, b, count = visited[thash]
        else:
            count = 0
        visited[thash] = [logl, tree.copy(), count+1]
        
        if isDebug(DEBUG_LOW):
            debug("\nuser given tree:")
            recon = phyloutil.reconcile(tree, stree, gene2species)
            events = phyloutil.labelEvents(tree, recon)
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
            tree = phyloutil.stree2gtree(top, labels, gene2species)
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
            recon = phyloutil.reconcile(tree, stree, gene2species)
            events = phyloutil.labelEvents(tree, recon)
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
    if True:
        trees = [x[1] for x in visited.itervalues()]
        i = util.argmax([x.data["logl"] for x in trees])
        return trees[i], trees[i].data["logl"]
    
    """
    # find best tree
    if False:
        errorfactor = 1.3
        minerror = min(x[1].data["error"] for x in visited.itervalues())

        # find all tree with error near minerror
        goodtrees = filter(lambda x: x[1].data["error"] < errorfactor * minerror,
                           visited.itervalues())

        # find best tree as max logl in good trees
        i = util.argmax([x[1].data["logl"] for x in goodtrees])

    if False:
        errorcutoff = conf["maxerror"]
        trees = [x[1] for x in visited.values()]
        errors = [tree.data["error"] for tree in trees]
        

        # find all tree with acceptable error
        goodind = util.find(lambda err: err < errorcutoff, errors)
        if len(goodind) > 0:
            goodtrees = util.mget(trees, goodind)
        else:
            # default to all trees if all errors are high
            debug("WARNING: high error rate in all trees found")
            goodtrees = trees

        # find best tree as max logl in good trees
        i = util.argmax([x.data["logl"] for x in goodtrees])
        return goodtrees[i], goodtrees[i].data["logl"]
    
    if False:
        # find best consensus tree
        items = visited.values()
        items.sort(key=lambda x: x[1].data["logl"], reverse=True)
        mat = [[x[1], x[2]] for x in items[:conf["toptrees"]]]
        trees, counts = zip(* mat)

        phylip.writeBootTrees(conf["out"] + ".trees", trees, counts=counts)
        tree = phylip.consense(trees, counts=counts, verbose=False, args="r\ny")
        return tree, 0
    """


def consensusTree(trees, counts):
    splits = util.Dict(default=0)
    
    genes = util.sort(trees[0].leaveNames())
    
    # count up splits
    for tree, count in zip(trees, counts):
        network = treelib.tree2graph(treelib.unroot(tree))
        splits2 = findSplits(network, util.makeset(tree.leaveNames()))
        
        print len(splits2)
        
        for key, (set1, set2) in splits2.iteritems():
            if len(set1) > len(set2):
                set1, set2 = set2, set1
            splitkey = tuple([int(gene in set1) for gene in genes])
            splits[splitkey] += count
    
    splits = splits.items()
    splits.sort(key=lambda x: x[1], reverse=True)
    
    half = len(trees) / 2.0
    if util.count(lambda x: x[1] >= half, splits):
        debug("consensus exists")
    
    # print splits
    if isDebug(DEBUG_LOW):
        mat = [genes + ["COUNT"]]
        for key, val in splits:
            mat.append(list(key))
            mat[-1].append(val)
        util.printcols(mat, out=DEBUG)

    
    return tree, tree.data["logl"]





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




