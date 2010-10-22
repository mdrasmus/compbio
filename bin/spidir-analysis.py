#!/usr/bin/env python

import sys, os, random
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

fitting = tablelib.Table(headers=["name", "type", "param1", "param2", "pval",
                                  "fit"])


def chisqFit(name, data, func, low, high, step, perc):
    bins, obs = util.hist(data, low=low, width=step)
    tot = len(data)
    ext = map(lambda x: perc * tot * (func(x+step) - func(x)),
              bins)
    
    ind = util.find(lambda x: x>5, ext)
    obs = util.mget(obs, ind)
    bins = util.mget(bins, ind)
    ext = util.mget(ext, ind)
    
    try:
        chisq = rpy.r.chisq_test(obs, p=util.oneNorm(ext))
        pval = float(chisq["p.value"])
    except rpy.RException, e:
        print e
        pval = 0.0
    
    qqplot = util.Gnuplot()
    qqplot.enableOutput(False)
    qqplot.plot(util.sort(obs), util.sort(ext),
               xlab="data", ylab="samples",
               main="%s abs qq-plot" % name)
    qqplot.plotdiag()
    
    return pval, qqplot


def ksFitting(name, data, func, low, high):
    obs = util.unique(data)
    samples = []    
    while len(samples) < len(obs):
        s = func()
        if low < s < high:
            samples.append(s)
    
    try:
        ret = rpy.r.ks_test(obs, samples)
        pval = float(ret["p.value"])
    except rpy.RException, e:
        print "KS:", e
        pval = 0.0
        
    
    qqplot = util.Gnuplot()
    qqplot.enableOutput(False)
    qqplot.plot(util.sort(obs), util.sort(samples),
               xlab="data", ylab="samples",
               main="%s abs qq-plot" % name)
    qqplot.plotdiag()

    
    return pval, qqplot



def plotAbsLens(name, lens, low, high, high2, step):
    p = util.Gnuplot()
    p.enableOutput(False)
    
    lens2 = filter(util.withinfunc(low, high), lens)
    lens3 = filter(util.withinfunc(low, high2), lens)
    perc = len(lens2) / float(len(lens))
    
    mu = stats.mean(lens2)
    sigma2 = stats.variance(lens2)
    prms = [mu*mu/sigma2, mu/sigma2]
    
    try:
        pass
        prms, resid = stats.fitDistrib(stats.gammaPdf, prms, 
                                       lens2, low, high, step, perc)
        prms = map(float, prms)
        resid = float(resid)        
    except:
        print "FIT EXCEPTION"
        pass
    
    util.plotdistrib(lens, low=0, width=step, plot=p)
    p.plotfunc(lambda x: stats.gammaPdf(x, prms), 0, high, step/4.)
    
    # goodnes of fit
    pval, qqplot = ksFitting(name, lens3, 
                    lambda: random.gammavariate(prms[0], 1./prms[1]), low, high2)
    
    util.logger("%s\t%.3e\t%s" % (name, pval, str(pval > .05)))
    
    p.set(xmin=0, xmax=high, 
          xlab="sub/site",
          main="%s absolute branch lengths (a=%.4f, b=%.4f, Pval=%.3e)" % 
          (str(name), prms[0], prms[1], pval))
    
    return p, util.Bundle(pval=pval, prms=prms, qqplot=qqplot)


def plotRelLens(name, params, lens, low, high, high2, step):
    p = util.Gnuplot()
    p.enableOutput(False)
    
    lens2 = filter(util.withinfunc(low, high), lens)
    lens3 = filter(util.withinfunc(low, high2), lens)
    perc = len(lens2) / float(len(lens))
    
    #prms, resid = stats.fitDistrib(stats.normalPdf, params, 
    #                               lens2, low, high, step)
    prms = params
    prms = map(float, prms)
    #resid = float(resid)  
    util.plotdistrib(lens, low=0, width=step, plot=p)
    p.plotfunc(lambda x: stats.normalPdf(x, prms), 0, high, step/4.)
    
    # goodness of fit
    pval, qqplot = ksFitting(name, lens2, 
                    lambda: random.normalvariate(prms[0], prms[1]), low, high)
    pval2, junk = ksFitting(name, lens3, 
                    lambda: random.normalvariate(prms[0], prms[1]), low, high2)
    #pval, qqplot = chisqFit(name, lens, 
    #                        lambda x: stats.normalCdf(x, prms), 
    #                        low, high, step, perc)
    
    util.logger("%s\t%.3e\t%.3e\t%s" % (name, pval, pval2, str(pval2 > .05)))
    
    p.set(xmin=0, xmax=high, 
          xlab="rel sub/site",
          main="%s relative branch lengths (mean=%.4f, sdev=%.4f, pval=%.3e)" % 
          (str(name), params[0], params[1], pval2))
    
    return p, util.Bundle(pval=pval, pval2=pval2, prms=prms, qqplot=qqplot)


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
os.system("convert %s/params-tree.svg %s/params-tree.png" %
          (conf["outdir"], conf["outdir"]))


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
    names = set(stree.nodes) & set(lens)

    for name in names:        
        #high = 3 * stats.median(lens[name])
        vals = util.sort(lens[name])
        low = vals[int(len(vals) * .015)]
        high2 = vals[int(len(vals) * .90)]
        high = 3 * stats.median(lens[name])
        
        if high == 0:
            high = high2
        
        step = (high - low) / conf["nbins"]
        p, fit = plotAbsLens(name, lens[name], low, high, high2, step)
        p.enableOutput()
        p.save(os.path.join(conf["outdir"], "abs/%s.ps" % str(name)))
        p.save(os.path.join(conf["outdir"], "abs/%s.png" % str(name)))
        
        fit.qqplot.enableOutput()
        fit.qqplot.save(os.path.join(conf["outdir"], 
                                     "abs/%s-qqplot.ps" % str(name)))
        
        fitting.append({"name": str(name),
                        "type": "abs",
                        "param1": fit.prms[0],
                        "param2": fit.prms[1],
                        "pval": fit.pval,
                        "fit": fit.pval >= .05})
    util.toc()


if 1:
    # plot all rel branch distributions
    util.tic("plot relative branch lengths")
    for name in stree.nodes:
        if name not in lens:
            util.log("skipping rel '%s'" % str(name))
            continue
        
        #high = 3 * stats.mean(rlens[name])
        vals = util.sort(lens[name])
        high =  params[name][0] + 4 * params[name][1]
        low = params[name][0] - 2 * params[name][1]
        high2 = params[name][0] + 2 * params[name][1]
        
        if high == 0: 
            high = high2
        
        step = (high - low) / conf["nbins"]
        p, fit = plotRelLens(name, params[name], rlens[name], low, high, high2, step)
        
        p.enableOutput()
        p.save(os.path.join(conf["outdir"], "rel/%s.ps" % str(name)))
        p.save(os.path.join(conf["outdir"], "rel/%s.png" % str(name)))
        
        fit.qqplot.enableOutput()
        fit.qqplot.save(os.path.join(conf["outdir"], 
                                     "rel/%s-qqplot.ps" % str(name)))
        
        fitting.append({"name": str(name),
                        "type": "rel",
                        "param1": fit.prms[0],
                        "param2": fit.prms[1],
                        "pval": fit.pval2,
                        "fit": fit.pval2 >= .05})
        
    util.toc()


fitting.write(os.path.join(conf["outdir"], "fitting.tab"))
print >>file(os.path.join(conf["outdir"], "fitting.tab.txt"), "w"), \
      repr(fitting)



# total branch length
if 1:
    util.log("plot total tree lengths")
    low = 0
    high = 3 * stats.mean(totals)
    step = (high - low) / conf["nbins"]
    p, fit = plotAbsLens("family rate", totals, low, high, max(totals), step)
    p.enableOutput()
    p.save(os.path.join(conf["outdir"], "family-rates.ps"))
    p.save(os.path.join(conf["outdir"], "family-rates.png"))


#####################################################
# correlation matrix

# get in order traversal of nodes
keys = []
def walk(node):
    if node.is_leaf():
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

if 1:
    mus = map(stats.mean, mat)
    sigmas = map(stats.sdev, mat)
    sel = [sum(int(x > mu + 6 * sigma) for x, mu, sigma in zip(col, mus, sigmas)) == 0
       for col in zip(* mat)]
    ind = util.findeq(True, sel)

    print "%f of trees kept" % (len(ind) / float(len(sel)))

    mat = [util.mget(row, ind) for row in mat]
    rmat = [util.mget(row, ind) for row in rmat]



colormap = util.ColorMap([[-1, util.red],
                          [-.5, util.orange],
                          [-.3, util.yellow],
                          [-.2, util.green],
                          [-.1, util.blue],
                          [0, util.blue],
                          [.1, util.blue],
                          [.2, util.green],
                          [.3, util.yellow],
                          [.5, util.orange],
                          [1, util.red]])

# write correlation legend
util.makeColorLegend(os.path.join(conf["outdir"], "corr/colorscale.svg"), 
                     colormap, -1, 1, .1, 
                     width=100, height=10)


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
keys = stree.leaf_names()
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

util.tic("rel paths correlation")
keys = stree.leaf_names()
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
