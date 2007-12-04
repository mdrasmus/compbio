
from rasmus.bio import phylo

import Spidir
from Spidir import Search
from Spidir.Debug import *


# events
EVENT_GENE = 0
EVENT_SPEC = 1
EVENT_DUP = 2

# fractional branches
FRAC_NONE = 0
FRAC_DIFF = 1
FRAC_PARENT = 2
FRAC_NODE = 3




#=============================================================================
# gene rate estimation

def mleBaserate(lens, means, sdevs, baserateparam):
    [alpha, beta] = baserateparam
    
    # use only best means and sdevs (highest means)
    ind = range(len(means))
    ind.sort(lambda a, b: cmp(means[b], means[a]))
    ind = ind[:max(4, len(ind) / 2 + 1)]
    #ind = ind[:max(4, len(means)-2)]
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
    
    roots = stats.solveCubic(a, b, c)
    
    def like(generate):
        if generate < 0:
            return -util.INF
        prod = 0
        for l, u, s in zip(lens, means, sdevs):
            prod += log(stats.normalPdf(l / generate, [u, s*s]))
        return log(stats.gammaPdf(generate, baserateparam)) + prod
    return roots[util.argmax(roots, like)]



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
    return baserate


#=============================================================================
# branch length likelihood


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


def rareEventsLikelihood(conf, tree, stree, recon, events):
    logl = 0.0
    
    for node, event in events.items():
        if recon[node] == stree.root and \
           event == "dup":
            logl += log(conf["predupprob"])
        
        if event == "dup":
            logl += log(conf["dupprob"])
        
    #nloss = len(phylo.findLoss(tree, stree, recon))
    #logl += nloss * log(conf["lossprob"])
    
    return logl


#-------------------------------------------------------------------------------
# Likelihood calculation
#-------------------------------------------------------------------------------



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





def subtreeLikelihood(conf, root, recon, events, stree, params, baserate,
                      integration="fastsampling"):
    midpoints = {}
    extraBranches = getExtraBranches(root, recon, events, stree)
    
    this = util.Closure(ncalls=0, depth=0, printing=0)    
    
    
    if integration == "fastsampling":
        # do fast integration
        logl3 = 0.0
        midpoints[root] = 1.0
        for child in root.children:
            # integration is only needed if child is dup
            if events[child] != "dup":
                                
                node = child
                if recon[node] != stree.root:
                    startparams, startfrac, midparams, \
                        endparams, endfrac, kdepend = \
                        reconBranch(node, recon, events, params)
                
                    setMidpoints(child, events, recon, midpoints, [])            
                    clogl = branchLikelihood(node.dist / baserate, 
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
                        reconBranch(node, recon, events, params)

                    if events[node] == "dup":
                        for child in node.children:
                            walk(child)
                walk(child)
                

                for samples in [100]:
                    val = 0.0            

                    for i in xrange(samples):
                        setMidpointsRandom(child, events, recon, midpoints)                

                        val2 = 0.0
                        for node in nodes:
                            if recon[node] != stree.root:
                                v = branchLikelihood(node.dist / baserate, 
                                              node, midpoints, 
                                              startparams[node], startfrac[node],
                                              midparams[node], endparams[node], 
                                              endfrac[node])
                                val2 += v

                        val += math.exp(val2)
                    clogl = log(val / float(samples))
                    
            child.data["logl"] = clogl
            logl3 += clogl 
        logl = logl3
    
    return logl


def reconBranch(node, recon, events, params):

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


def branchLikelihood(dist, node, k, startparams, startfrac,
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



def setMidpointsRandom(node, events, recon, midpoints, wholeTree=False):
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

def treeLogLikelihood_python(conf, tree, stree, gene2species, params, 
                             baserate=None, integration="fastsampling"):

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

    phylo.midrootRecon(tree, stree, recon, events, params, baserate)
    
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
    tree.data["errorlogl"] = tree.data.get("error", 0.0) * \
                             conf.get("errorcost", 0.0)
    this.logl += tree.data["errorlogl"]

    # family rate likelihood
    if conf["famprob"]:
        this.logl += log(stats.gammaPdf(baserate, params["baserate"]))
    
    tree.data["baserate"] = baserate
    tree.data["logl"] = this.logl
    
    
    if isDebug(DEBUG_MED):
        util.toc()
        debug("\n\n")
        drawTreeLogl(tree, events=events)
    
    return this.logl


