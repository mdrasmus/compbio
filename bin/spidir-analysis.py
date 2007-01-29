#!/usr/bin/env python

import sys, os
from rasmus import util, env, treelib, stats, tablelib
import Spidir

import rpy


options = [
  ["s:", "stree=", "stree", "<species tree>",
    {"single": True}],
  ["p:", "params=", "params", "<sindir params file>",
    {"single": True}],
  ["l:", "lens=", "lens", "<sindir branch lengths file>",
    {"single": True}],
  ["o:", "outdir=", "outdir", "<output directory>",
    {"single": True,
     "default": "."}],
  ["b:", "nbins=", "nbins", "<number of bins>",
    {"single": True,
     "parser": int,
     "default": 20}]
]


# parse options
conf = util.parseOptions(sys.argv, options, quit=True)


# read data
env.addEnvPaths("DATAPATH")
stree = treelib.readTree(env.findFile(conf["stree"]))
params = Spidir.readParams(env.findFile(conf["params"]))
lens = Spidir.readTreeDistrib(env.findFile(conf["lens"]))
totals = map(sum, zip(* lens.values()))
rlens = util.mapdict(lens, valfunc=lambda x: util.vidiv(x, totals))

fitting = tablelib.Table(headers=["name", "type", "param1", "param2", "pval"])


def plotAbsLens(name, lens, low, high, step):
    p = util.Gnuplot()
    p.enableOutput(False)
    
    lens = filter(util.withinfunc(low, high), lens)
    
    lens2 = filter(util.gtfunc(.001), lens)
    mu = stats.mean(lens2)
    sigma2 = stats.variance(lens2)
    prms = [mu*mu/sigma2, mu/sigma2]
    
    try:
        prms, resid = stats.fitDistrib(stats.gammaPdf, prms, 
                                       lens2, low, high, step)
        prms = map(float, prms)
        resid = float(resid)
    except:
        print "FIT EXCEPTION"
        pass
    
    util.plotdistrib(lens, low=low, width=step, plot=p)
    p.plotfunc(lambda x: stats.gammaPdf(x, prms), low, high, step/4.)
    
    
    bins, obs = util.hist(lens, low=low, width=step)
    tot = len(lens)
    ext = map(lambda x: tot * (stats.gammaCdf(x+step, prms) - 
                               stats.gammaCdf(x, prms)), bins)
    
    ind = util.find(lambda x: x>5, ext)
    ind = ind[1:]
    obs = util.mget(obs, ind)
    bins = util.mget(bins, ind)
    ext = util.mget(ext, ind)
    
    try:
        chisq = rpy.r.chisq_test(obs, p=util.oneNorm(ext))
        pval = float(chisq["p.value"])
    except rpy.RException, e:
        print e
        pval = 0.0
    
    
    util.logger("%s\t%.3e\t%s" % (name, pval, str(pval > .05)))
    
    p.set(xmin=low, xmax=high, 
          xlab="sub/site",
          main="%s absolute branch lengths (a=%f, b=%f, Pval=%e)" % 
          (str(name), prms[0], prms[1], pval))
    
    return p, util.Bundle(pval=pval, prms=prms)


def plotRelLens(name, params, lens, low, high, step):
    p = util.Gnuplot()
    p.enableOutput(False)
    
    lens = filter(util.withinfunc(low, high), lens)
    lens2 = filter(util.withinfunc(.00001, high), lens)
    
    prms, resid = stats.fitDistrib(stats.normalPdf, params, 
                                   lens2, low, high, step)
    prms = map(float, prms)
    resid = float(resid)  
    util.plotdistrib(lens, low=low, width=step, plot=p)
    p.plotfunc(lambda x: stats.normalPdf(x, prms), low, high, step/4.)
    
    bins, obs = util.hist(lens, low=low, width=step)
    tot = len(lens)
    ext = map(lambda x: tot * (stats.normalCdf(x+step, prms) - 
                               stats.normalCdf(x, prms)), bins)
    
    ind = util.find(lambda x: x>5, ext)
    ind = ind[1:]
    obs = util.mget(obs, ind)
    bins = util.mget(bins, ind)
    ext = util.mget(ext, ind)
    
    
    try:
        chisq = rpy.r.chisq_test(obs, p=util.oneNorm(ext))
        pval = float(chisq["p.value"])
    except rpy.RException, e:
        print e
        pval = 0.0
    
    #pval = rpy.r.shapiro_test(lens2)["p.value"]
    
    util.logger("%s\t%.3e\t%s" % (name, pval, str(pval > .05)))
    
    #p.data[1].options["plab"] = "best fit"
    p.set(xmin=low, xmax=high, 
          xlab="rel sub/site",
          main="%s relative branch lengths (mean=%f, sdev=%f, pval=%e)" % 
          (str(name), params[0], params[1], pval))
    
    #q=util.plot(bins, obs, style="lines")
    #q.plot(bins, ext, style="lines")
    #q.save("tmp_" + str(name) + ".ps")
    
    return p, util.Bundle(pval=pval, prms=prms)


def getCorrMatrix(lens, stree, leaves, useroot=True, get=None):
    corrmat = []

    ntrees = len(lens.values()[0])

    for leaf1 in leaves:
        corrmat.append([])
    
        for leaf2 in leaves:
            node1 = stree.nodes[leaf1]
            node2 = stree.nodes[leaf2]
            ancester = treelib.lca([node1, node2])
        
            dists1 = [0] * ntrees
            dists2 = [0] * ntrees

            
            stop = [ancester]
            if not useroot:
                stop.extend(stree.root.children)
            

            while node1 not in stop:
                for i in range(ntrees):
                    dists1[i] += lens[node1.name][i]
                node1 = node1.parent

            while node2 not in stop:
                for i in range(ntrees):
                    dists2[i] += lens[node2.name][i]
                node2 = node2.parent

            if get != None and \
               get == (leaf1, leaf2):
                return dists1, dists2
            
            c = stats.corr(dists1, dists2)

            if c == util.INF:
                c = 1.0
            corrmat[-1].append(c)

    return corrmat



    

# output params
os.system("mkdir -p %s" % conf["outdir"])
os.system("viewparam.py -p %s -s %s -l 500 > %s/params.txt" %
          (conf["params"], conf["stree"], conf["outdir"]))
Spidir.drawParamTree(stree, params, xscale=2000, 
                     filename="%s/params-tree.svg" % conf["outdir"])

# make table version of parameters
tab = tablelib.Table(headers=["branch", "mean", "sdev", "tight"])
for name in params:
    if name in (1, "baserate"):
        continue
    tab.add(branch=str(name),
            mean=params[name][0],
            sdev=params[name][1],
            tight=params[name][1] / params[name][0])
tab.sort(col="branch")
tab.write("%s/params.tab" % conf["outdir"])

    

# make directories
os.system("mkdir -p %s/abs" % conf["outdir"])
os.system("mkdir -p %s/rel" % conf["outdir"])
os.system("mkdir -p %s/corr" % conf["outdir"])
os.system("mkdir -p %s/corr/abs" % conf["outdir"])
os.system("mkdir -p %s/corr/rel" % conf["outdir"])

if 1:
    # plot all abs branch distributions
    util.tic("plot absolute branch lengths")
    for name in stree.nodes:
        if name not in lens:
            util.log("skipping abs '%s'" % str(name))
            continue
        
        low = 0
        high = 3 * stats.mean(lens[name])
        step = (high - low) / conf["nbins"]
        p, fit = plotAbsLens(name, lens[name], low, high, step)
        p.enableOutput()
        p.save(os.path.join(conf["outdir"], "abs/%s.ps" % str(name)))
        
        fitting.append({"name": str(name),
                        "type": "abs",
                        "param1": fit.prms[0],
                        "param2": fit.prms[1],
                        "pval": fit.pval})
        
    util.toc()

if 1:

    # plot all rel branch distributions
    util.tic("plot relative branch lengths")
    for name in stree.nodes:
        if name not in lens:
            util.log("skipping rel '%s'" % str(name))
            continue
        
        low = 0
        high = 3 * stats.mean(rlens[name])
        step = (high - low) / conf["nbins"]
        p, fit = plotRelLens(name, params[name], rlens[name], low, high, step)
        p.enableOutput()
        p.save(os.path.join(conf["outdir"], "rel/%s.ps" % str(name)))
        
        fitting.append({"name": str(name),
                        "type": "rel",
                        "param1": fit.prms[0],
                        "param2": fit.prms[1],
                        "pval": fit.pval})
        
    util.toc()


fitting.write(os.path.join(conf["outdir"], "fitting.tab"))


# total branch length
util.log("plot total tree lengths")
low = 0
high = 3 * stats.mean(totals)
step = (high - low) / conf["nbins"]
p, fit = plotAbsLens("family rate", totals, low, high, step)
p.enableOutput()
p.save(os.path.join(conf["outdir"], "family-rates.ps"))


#####################################################
# correlation matrix

# get in order traversal of nodes
keys = []
def walk(node):
    if node.isLeaf():
        keys.append(node.name)
    else:
        walk(node.children[0])
        if node.name in lens:
            keys.append(node.name)
        walk(node.children[1])
walk(stree.root)

#
# abs all branches
#
mat = util.mget(lens, keys)
rmat = util.mget(rlens, keys)

colormap = util.ColorMap([[-1, util.red],
                          [-.5, util.orange],
                          [-.3, util.yellow],
                          [-.2, util.green],
                          [-.1, util.blue],
                          [0, (1, 1, 1, 1)],
                          [.1, util.blue],
                          [.2, util.green],
                          [.3, util.yellow],
                          [.5, util.orange],
                          [1, util.red]])




util.tic("abs correlation")
corrmat = stats.corrmatrix(mat)
util.heatmap(corrmat, width=20, height=20,
             rlabels = keys, clabels=keys,
             xmargin=100, ymargin=100,
             labelSpacing=8,
             labelPadding=4,
             colormap=colormap,
             showVals=True,
             display=False,
             filename=os.path.join(conf["outdir"], "corr/abs.svg"))
util.toc()

# make png
os.system("convert %s/corr/abs.svg %s/corr/abs.png" % 
          (conf["outdir"], conf["outdir"]))

# write table
tab = tablelib.Table(corrmat, headers=map(str, keys))
tab.write(os.path.join(conf["outdir"], "corr/abs.tab"))


#
# rel all branches
#
util.tic("rel correlation")
rcorrmat = stats.corrmatrix(rmat)
util.heatmap(rcorrmat, width=20, height=20,
             rlabels = keys, clabels=keys,
             xmargin=100, ymargin=100,
             labelSpacing=8,
             labelPadding=4,
             colormap=colormap,
             showVals=True,
             display=False,
             filename=os.path.join(conf["outdir"], "corr/rel.svg"))
util.toc()

# make png
os.system("convert %s/corr/rel.svg %s/corr/rel.png" % 
          (conf["outdir"], conf["outdir"]))

# write table
tab = tablelib.Table(rcorrmat, headers=map(str, keys))
tab.write(os.path.join(conf["outdir"], "corr/rel.tab"))



#
# abs paths
#

util.tic("abs paths correlation")
keys = stree.leafNames()
corrmat = getCorrMatrix(lens, stree, keys, False)
util.heatmap(corrmat, width=20, height=20,
             rlabels = keys, clabels=keys,
             xmargin=100, ymargin=100,
             labelSpacing=8,
             labelPadding=4,
             colormap=colormap,
             showVals=True,
             display=False,
             filename=os.path.join(conf["outdir"], "corr/abs_paths.svg"))
util.toc()

# make png
os.system("convert %s/corr/abs_paths.svg %s/corr/abs_paths.png" % 
          (conf["outdir"], conf["outdir"]))

# write table
tab = tablelib.Table(corrmat, headers=map(str, keys))
tab.write(os.path.join(conf["outdir"], "corr/abs_paths.tab"))


#
# rel paths
#

util.tic("abs rel correlation")
keys = stree.leafNames()
rcorrmat = getCorrMatrix(rlens, stree, keys, False)
util.heatmap(rcorrmat, width=20, height=20,
             rlabels = keys, clabels=keys,
             xmargin=100, ymargin=100,
             labelSpacing=8,
             labelPadding=4,
             colormap=colormap,
             showVals=True,
             display=False,
             filename=os.path.join(conf["outdir"], "corr/rel_paths.svg"))
util.toc()

# make png
os.system("convert %s/corr/rel_paths.svg %s/corr/rel_paths.png" % 
          (conf["outdir"], conf["outdir"]))

# write table
tab = tablelib.Table(corrmat, headers=map(str, keys))
tab.write(os.path.join(conf["outdir"], "corr/rel_paths.tab"))



sys.exit(0)

################################################################
util.tic("abs scatter plots")
for i in range(len(keys)):
    for j in range(i+1, len(keys)):
        util.log("plot %s vs %s" % (str(keys[i]), str(keys[j])))
        
        p = util.Gnuplot()
        p.enableOutput(False)
        p.plot(lens[keys[i]], lens[keys[j]],
               xmax = 3 * stats.mean(lens[keys[i]]),
               ymax = 3 * stats.mean(lens[keys[j]]),
               xmin = 0,
               ymin = 0,
               main = "abs branch length correlation (corr=%f)" %
                      stats.corr(lens[keys[i]], lens[keys[j]]),
               xlab = str(keys[i]),
               ylab = str(keys[j]))
        p.enableOutput()
        p.save(os.path.join(conf["outdir"], "corr/abs/%s_%s.ps" % 
                            (str(keys[i]), str(keys[j]))))
util.toc()


util.tic("rel scatter plots")
for i in range(len(keys)):
    for j in range(i+1, len(keys)):
        util.log("plot %s vs %s" % (str(keys[i]), str(keys[j])))
        
        p = util.Gnuplot()
        p.enableOutput(False)
        p.plot(rlens[keys[i]], rlens[keys[j]],
               xmax = 3 * stats.mean(rlens[keys[i]]),
               ymax = 3 * stats.mean(rlens[keys[j]]),
               xmin = 0,
               ymin = 0,
               main = "rel branch length correlation (corr=%f)" %
                      stats.corr(rlens[keys[i]], rlens[keys[j]]),
               xlab = str(keys[i]),
               ylab = str(keys[j]))
        p.enableOutput()
        p.save(os.path.join(conf["outdir"], "corr/rel/%s_%s.ps" % 
                            (str(keys[i]), str(keys[j]))))
util.toc()
