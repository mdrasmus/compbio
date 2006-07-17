#!/usr/bin/env python

import sys, os
from rasmus import util, env, treelib, sindirlib, stats


options = [
  ["s:", "stree=", "stree", "<species tree>",
    {"single": True}],
  ["p:", "params=", "params", "<sindir params file>",
    {"single": True}],
  ["l:", "lens=", "lens", "<sindir branch lengths file>",
    {"single": True}],
  ["o:", "outdir=", "outdir", "<output directory>",
    {"single": True,
     "default": "."}]
]


# parse options
conf = util.parseOptions(sys.argv, options, quit=True)


# read data
env.addEnvPaths("DATAPATH")
stree = treelib.readTree(env.findFile(conf["stree"]))
params = sindirlib.readParams(env.findFile(conf["params"]))
lens = sindirlib.readTreeDistrib(env.findFile(conf["lens"]))
totals = map(sum, zip(* lens.values()))
rlens = util.mapdict(lens, valfunc=lambda x: util.vidiv(x, totals))



def plotAbsLens(name, lens, low, high, step):
    p = util.Gnuplot()
    p.enableOutput(False)
    p, prms, resid = stats.plotdistribFit(stats.gammaPdf, [1,1], 
                                          lens, low, high, step, plot=p)
    p.set(xmin=low, xmax=high, 
          xlab="sub/site",
          main="%s absolute branch lengths (a=%f, b=%f, resid=%f)" % 
          (str(name), prms[0], prms[1], resid))
    
    return p


def plotRelLens(name, params, lens, low, high, step):
    p = util.Gnuplot()
    p.enableOutput(False)
    p, prms, resid = stats.plotdistribFit(stats.normalPdf, params, 
                                          lens, low, high, step, plot=p)
    p.plotfunc(lambda x: stats.normalPdf(x, params), low, high, step,
               plab = "sindir param")
    p.data[1].options["plab"] = "best fit"
    p.set(xmin=low, xmax=high, 
          xlab="rel sub/site",
          main="%s relative branch lengths (mean=%f, sdev=%f, resid=%f)" % 
          (str(name), params[0], params[1], resid))
    
    return p
    

# output params
os.system("viewparam.py -p %s -s %s -l 500 > %s/params.txt" %
          (conf["params"], conf["stree"], conf["outdir"]))


# make directories
os.system("mkdir -p %s/abs" % conf["outdir"])
os.system("mkdir -p %s/rel" % conf["outdir"])
os.system("mkdir -p %s/corr" % conf["outdir"])
os.system("mkdir -p %s/corr/abs" % conf["outdir"])
os.system("mkdir -p %s/corr/rel" % conf["outdir"])


# plot all abs branch distributions
util.tic("plot absolute branch lengths")
for name in stree.nodes:
    if name not in lens:
        util.log("skipping abs '%s'" % str(name))
        continue

    util.log("plot abs '%s'" % str(name))

    low = 0
    high = 3 * stats.mean(lens[name])
    step = (high - low) / 50.
    p = plotAbsLens(name, lens[name], low, high, step)
    p.enableOutput()
    p.save(os.path.join(conf["outdir"], "abs/%s.ps" % str(name)))
util.toc()


# plot all rel branch distributions
util.tic("plot relative branch lengths")
for name in stree.nodes:
    if name not in lens:
        util.log("skipping abs '%s'" % str(name))
        continue

    util.log("plot rel '%s'" % str(name))

    low = 0
    high = 3 * stats.mean(rlens[name])
    step = (high - low) / 50.
    p = plotRelLens(name, params[name], rlens[name], low, high, step)
    p.enableOutput()
    p.save(os.path.join(conf["outdir"], "rel/%s.ps" % str(name)))
util.toc()



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

mat = util.mget(lens, keys)
rmat = util.mget(rlens, keys)

colormap = util.ColorMap([[-.1, util.blue],
                          [-.5, util.green],
                          [-.75, util.yellow],
                          [-1, util.red],
                          
                          [0, (0, 0, .7, 1)],
                          [.1, util.blue],
                          [.5, util.green],
                          [.75, util.yellow],
                          [1, util.red]])



util.tic("abs correlation")
util.heatmap(stats.corrmatrix(mat), width=20, height=20,
             rlabels = keys, clabels=keys,
             xmargin=100, ymargin=100,
             labelSpacing=8,
             labelPadding=4,
             colormap=colormap,
             showVals=True,
             display=False,
             filename=os.path.join(conf["outdir"], "corr/abs.svg"))
util.toc()

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


util.tic("rel correlation")
util.heatmap(stats.corrmatrix(rmat), width=20, height=20,
             rlabels = keys, clabels=keys,
             xmargin=100, ymargin=100,
             labelSpacing=8,
             labelPadding=4,
             colormap=colormap,
             showVals=True,
             display=False,
             filename=os.path.join(conf["outdir"], "corr/rel.svg"))
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
