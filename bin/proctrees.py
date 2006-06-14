#!/usr/bin/env python

from rasmus import util, algorithms, fasta, phyloutil, ensembl

import sys


options = [
    ["b:", "boot=", "boot", "AUTO<iters>"],
    ["t:", "stree=", "stree", "AUTO<species tree>"],
    ["d:", "dir=", "dir", "AUTO<trees directory>"],
    ["P", "nopart", "nopart", "AUTO"],
    ["s:", "minsize=", "minsize", "AUTO<minimum size trees to process"]
    ]
    

try:
    param, rest = util.parseArgs(sys.argv, options, "bootstrapped trees")
except:
    sys.exit(1)

if "boot" in param:
    iters = int(param["boot"][-1])


# read species tree
stree = algorithms.Tree()
stree.readNewick(param["stree"][-1])
genomes = stree.leaveNames()


if "dir" in param:
    files = util.listFiles(param["dir"][-1], ".boot.tree")
else:
    files = rest
    
if "minsize" in param:
    minsize = int(param["minsize"][-1])
else:
    minsize = 0

#util.globalTimer().maxdepth = 0

homology = phyloutil.Homology()

# process each gene tree
for treefile in files:
    util.tic("process %s" % treefile)
    
    infile = file(treefile)
    out = file(treefile.replace(".tree", ".part.tree"), "w")
    
    # process bootstrapped trees
    for i in xrange(iters):
        util.tic("iter %d" % i)
        
        # read tree
        tree = algorithms.Tree()
        tree.readNewick(infile)
        #print "full tree"
        #tree.write()
        
        # quick processing
        if len(tree.leaves()) < minsize:
            genes = tree.leaveNames()
            parts = util.Dict(1, [])
            for gene in genes:
                parts[ensembl.id2genome(gene)].append(gene)
            
            # count orthologs
            for a in xrange(len(parts)):
                for b in xrange(a+1, len(parts)):
                    for gene1 in parts.values()[a]:
                        for gene2 in parts.values()[b]:
                            homology.incOrtholog(gene1, gene2, 1)
            
            # count paralogs
            for part in parts.values():
                for a in xrange(len(part)):
                    for b in xrange(a+1, len(part)):
                        homology.incParalog(part[a], part[b], 1)
            
            util.toc()
            break
        
        # root tree and partition it
        tree = phyloutil.reconRoot(tree, stree)   
        
        if "nopart" in param:
            trees = [tree]
        else:     
            trees = phyloutil.partitionTree(tree, stree)
        
        util.log("parts %d" % len(trees))
        
        for tree2 in trees:
            if "nopart" in param:
                tree2 = phyloutil.reconRoot(tree2, stree)
            recon = phyloutil.reconcile(tree2, stree)
            
            tree2.writeNewick(out)
            
            #print "part tree"
            #tree2.write()
            #print "loss", len(phyloutil.findLoss(tree2, stree, recon))
            
            for gene in tree2.leaveNames():
                # count orthologs
                for genome in genomes:
                    if genome == ensembl.id2genome(gene):
                        continue
                    
                    ogenes = phyloutil.orthologs(gene, genome, tree2, 
                                                 stree, recon)
                    for gene2 in ogenes:
                        homology.incOrtholog(gene, gene2, .5 / iters)
                
                # count paralogs
                pgenes = phyloutil.paralogs(gene, tree2, recon)
                for gene2 in pgenes:
                    homology.incParalog(gene, gene2, .5 / iters)
        util.toc()   
    util.toc()        

# output
homology.write()

