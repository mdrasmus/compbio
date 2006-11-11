#!/usr/bin/env python

import os, sys

from rasmus import phyloutil, tablelib, sindirlib, spidirlib, util
from rasmus import treelib, genomeutil, env


options = [
    ["s:", "stree=", "stree", "<species tree>",
     {"single": True}],
    ["S:", "smap=", "smap", "<gene to species map>",
     {"single": True}],
    ["T:", "treeext=", "treeext", "<Tree extension>",
     {"single": True}],
    ["i", "stdin", "stdin", "",
     {"single": True}],
    ["p:", "params=", "params", "<spidir parameters file>",
     {"single": True}],
    ["P:", "parts=", "parts", "<partition file>",
     {"single": True}],
    ["o:", "output=", "output", "<output prefix>",
     {"single": True}],
    ]


conf = util.parseOptions(sys.argv, options, quit=True)



def drawEventTree(stree, out=sys.stdout):
    labels = {}
    for name, node in stree.nodes.iteritems():
        labels[name] = "[%s]\nD=%d,L=%d;\nG=%d;" % \
                       (str(name),
                        node.data['dup'], node.data['loss'],
                        node.data['genes'])

    treelib.drawTree(stree, labels=labels, minlen=15, spacing=4, labelOffset=-3,
                     out=out)



def countDupLoss(conf, treefiles, stree, gene2species, params):
    totalTree = stree.copy()
    singleTree = stree.copy()
    
    # setup table and headers
    headers = ["partid", "genes", "treelen"]
    if params:
        headers.append("famrate")
    headers += ["appear", "dup", "loss"]
    
    for name in stree.nodes:
        headers.extend(["%s-appear" % str(name),
                        "%s-dup" % str(name),
                        "%s-loss" % str(name),
                        "%s-genes" % str(name)])
    
    tab = tablelib.Table(headers=headers)
    
    
    
    # initalize counts to zero
    phyloutil.initDupLossTree(totalTree)
    
    
    util.tic("read trees")
    # count dup loss
    j = 0
    for i, f in enumerate(treefiles):
        if not os.path.exists(f):
            print "skipping", f
            continue
        print i, j, f
        j += 1
        
        tree = treelib.readTree(f)
        phyloutil.reconRoot(tree, stree, gene2species, newCopy=False)

        phyloutil.initDupLossTree(singleTree)
        dup, loss, appear = phyloutil.countDupLossTree(tree, singleTree, gene2species)
        dup, loss, appear = phyloutil.countDupLossTree(tree, totalTree, gene2species)

        phyloutil.countAncestralGenes(singleTree)
        
        row, junk = tree2row(singleTree)
        
        """
        row = {"partid": i,
               "genes": len(tree.leaves()),
               "appear": appear,               
               "dup": dup,
               "loss": loss
               }

        for name, node in singleTree.nodes.iteritems():
            row["%s-appear" % str(name)] = node.data['appear']
            row["%s-dup" % str(name)] = node.data['dup']
            row["%s-loss" % str(name)] = node.data['loss']
            row["%s-genes" % str(name)] = node.data['genes']
        """
        
        row["partid"] = i
        
        # add treelen
        row['treelen'] = sum(x.dist for x in tree.nodes.values())        
        
        # add family rate
        if params:
            row['famrate'] = spidirlib.getBaserate(tree, stree, params, \
                                                   gene2species=gene2species)
        
        tab.append(row)
    util.toc()

    # count ancestral genes
    phyloutil.countAncestralGenes(totalTree)
    
        
    return tab, totalTree




def readDupLossTree(stree, row):
    stree2 = stree.copy()

    for name, node in stree2.nodes.items():
        node.data['dup'] = row["%s-dup" % str(name)]
        node.data['loss'] = row["%s-loss" % str(name)]
        node.data['genes'] = row["%s-genes" % str(name)]
    return stree2


def tree2row(stree):
    row = {"genes": 0,
           "appear": 0,               
           "dup": 0,
           "loss": 0
           }
    
    headers = ["genes", "appear", "dup", "loss"]

    for name, node in stree.nodes.iteritems():
        row["%s-appear" % str(name)] = node.data['appear']
        row["%s-dup" % str(name)] = node.data['dup']
        row["%s-loss" % str(name)] = node.data['loss']
        row["%s-genes" % str(name)] = node.data['genes']
        
        headers.extend(["%s-appear" % str(name),
                        "%s-dup" % str(name),
                        "%s-loss" % str(name),
                        "%s-genes" % str(name)])
        
        row['appear'] += node.data['appear']
        row['dup'] += node.data['dup']
        row['loss'] += node.data['loss']
        row['genes'] += node.data['genes']
    
    return row, headers
    

def addDupLossTree(stree, row):
    for name, node in stree.nodes.items():
        if 'dup' not in node.data:
            node.data['dup'] = 0
            node.data['loss'] = 0
            node.data['genes'] = 0
            
        node.data['dup'] += row["%s-dup" % str(name)]
        node.data['loss'] += row["%s-loss" % str(name)]
        node.data['genes'] += row["%s-genes" % str(name)]
    



if 0:
    tottree = stree.copy()
    for row in tab:
        if 17 < row['size'] < 40:
            addDupLossTree(tottree, row)

# count duplications after speciations
if 0:
    tab2 = copy.copy(tab)
    tab2.headers.append("post-dup")
    tab2.headers.append("dup/gene")

    for row in tab2:
        row["post-dup"] = row["dup"] - row["1-dup"]
        row["dup/gene"] = row["dup"] / float(row["size"])



def main(conf):
    print "here"

    # read data
    env.addEnvPaths("DATAPATH")
    gene2species = genomeutil.readGene2species(env.findFile(conf["smap"]))
    stree = treelib.readTree(env.findFile(conf["stree"]))

    if "params" in conf:
        params = sindirlib.readParams(env.findFile(conf["params"]))
    else:
        params = None
    
    #terms = readGo(env.findFile("orf_geneontology.tab"))
    #parts = util.readDelim(env.findFile(conf["parts"]))
    
    
    # determine treefiles
    treefiles = []
    treefiles.extend(conf["REST"])

    if conf["stdin"]:  
        for line in sys.stdin:
            treefiles.append(line.rstrip())


    #names = [util.replaceExt(os.path.basename(x), conf["treeext"], "")
    #         for x in treefiles]
    
    tab, totalTree = countDupLoss(conf, treefiles, stree, gene2species, params)
    
    # total table
    row, headers = tree2row(totalTree)
    tab2 = tablelib.Table([row], headers=headers)
    
    # output
    tab.write(conf["output"] + ".tab")    
    tab2.write(conf["output"] + "-total.tab")
    
    drawEventTree(totalTree, out=file(conf["output"] + ".txt", "w"))


main(conf)
