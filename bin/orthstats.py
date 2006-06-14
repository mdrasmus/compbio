#!/usr/bin/env python

import util
import sys
import ensembl
import genomeio, genomeutil
import descriptions


options = [
    ("p:", "proper=", "proper", "[ -p <orth file>] [--proper <orth file>]"),
    ("e:", "extended=", "extended", "[ -e <orth file>] [--extended <orth file>]"),    
    ("a:", "ambig=", "ambiguous", "[ -a <orth file>] [--ambig <orth file>]"),
    ("g:", "genome=", "genome", "[ -g <genome>] [--genome <genome>]"),
    ("o:", "go=", "go", "[ -o <go file>] [--go <go file>]"),
    ("s:", "synteny=", "synteny", "[-s <synteny comps>] [ --synteny <synteny comps> ]"),
    ("c", "components", "components", "[ -c ] [ --components ]"),
    ("m", "members", "members", "[ -m ] [ --members ]"),
    ("d:", "desc=", "desc", "[ -d <desc>] [ --desc <desc> ]"),
    ("z", "sizesort", "sizesort", "[-z] [-sizesort]"),
    ["x", "exons", "exons", "[-e] [--exons]"]
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)



stats = {
    "proper": [],
    "extended": [],
    "ambiguous": []
}


# read exons
matching = genomeutil.Matching()
genomes = param["genome"]
if "exons" in param:
    for genome in ["human"]:
        matching.genomes[genome] = genomeutil.Genome()
        matching.genomes[genome].name = genome
    
    genomeio.readGenomeExons(matching.genomes)
    
    print matching.genomes["human"].genes.values()[0].trans.values()[0].exons


for f in param.get("proper", []):
    stats["proper"].extend(util.readDelim(f))
for f in param.get("extended", []):
    stats["extended"].extend(util.readDelim(f))
for f in param.get("ambiguous", []):
    stats["ambiguous"].extend(util.readDelim(f))

syncomps = []
for f in param.get("synteny", []):
    syncomps.extend(util.readDelim(f))


order = param["genome"]

desc = {}
if "desc" in param:
    for f in param["desc"]:
        desc.update(genomeio.readGeneDesc(f))

go = {}
if "go" in param:
    for f in param["go"]:
        for line in file(f):
            words = line.split()
            go[words[0]] = words[1:]




def printComp(param, header, comp, order, desc, props):
    counts = ensembl.genomeComposition(order, comp)
    
    print "%s: %d genes" % (header, len(comp))

    words = descriptions.summary(comp, desc)
    keys = words.keys()
    keys.sort(lambda a, b: cmp(words[b], words[a]))
    
    print
    print "%s: size %d [%s]" % \
        (header,
         len(comp), 
         ", ".join(map(lambda x: "%s(%d)" % (x, words[x]), keys[:5])))
    for genome in order:
        print "%s(%d) " % (genome, counts[genome]),
    print
    
    print "---------------------------------------"
    for gene in comp:
        # properties
        prefix = ""
        if gene in props.syntenic:
            prefix += "S"
        else:
            prefix += " "
        if gene in props.proper:
            prefix += "P"
        else:
            prefix += " "
        
        extra = ""
        if gene in go:
            extra += str(go[gene]) + " "
        if gene in desc:
            extra += desc[gene]
        print "%s %20s\t%s" % (prefix, gene, extra)
        
    print "---------------------------------------"            
    print
    

class GeneProperty:
    def __init__(self, syncomps, stats):
        self.syntenic = {}
        self.proper   = {}
        
        for comp in syncomps:
            for gene in comp:
                self.syntenic[gene] = 1
        for comp in stats["proper"]:
            for gene in comp:
                self.proper[gene] = 1
        
  

# build property dict
props = GeneProperty(syncomps, stats)

print "KEY"
print "P: the gene appears in a proper group"
print "S: the gene is also syntenic"
print 

# summary
for name in ["proper", "extended", "ambiguous"]:
    sizes = map(len, stats[name])
    if len(sizes) > 0: low = min(sizes)
    else: low = 0
    if len(sizes) > 0: top = max(sizes)
    else: top = 0
    
    print "%9s %7d genes, %6d groups, min %5d, max %5d" % \
        (name, sum(sizes), len(stats[name]), low, top)
print


# summary per group type
for name in ["proper", "extended", "ambiguous"]:
    sizes = map(len, stats[name])
    if len(sizes) > 0: low = min(sizes)
    else: low = 0
    if len(sizes) > 0: top = max(sizes)
    else: top = 0

    print "---------------------------------------"
    print "%s %d genes, %d groups, min %d, max %d" % \
        (name, sum(sizes), len(stats[name]), low, top)
    print 
    
    if len(sizes) == 0:
        continue
    
    print "histogram of group sizes [size(many)]:"
    hist = util.histInt(sizes)
    stride = 8
    items = [(i, hist[i]) for i in xrange(len(hist))]
    items = filter(lambda x: x[1] > 0, items)
    for i in xrange(0, len(items), stride):
        for j in xrange(i, min(i+stride, len(items))):
            print "%d(%d) " % items[j],
        print 
    print
    
    
    print "histogram of group compositions:"
    print order
    compositions = ensembl.componentCompositions(order, stats[name])
    keys = compositions.keys()
    keys.sort(lambda a,b: cmp(compositions[b], compositions[a]))
    for i in keys:
        print "%s: %d" % (str(i), compositions[i])        
    print


# break down
for name in ["proper", "extended", "ambiguous"]:
    print "======================================="
    
    comps = stats[name]
    
    ind = range(len(comps))
    ind.sort(lambda a, b: cmp(len(comps[b]), len(comps[a])))
    
    for i in ind:
        printComp(param, name, comps[i], order, desc, props)



