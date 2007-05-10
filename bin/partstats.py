#!/usr/bin/env python

from rasmus import util, genomeio, genomeutil, fasta, env
import sys



options = [
    ["p:", "part=", "part", "<part file>"],
    ["g:", "genomes=", "genomes", "<genome>"],
    ["o:", "go=", "go", "<go file>"],
    ["s", "summary", "summary", ""],
    ["c", "components", "components", ""],
    ["m", "members", "members", ""],
    ["t", "details", "details", ""],
    ["d:", "desc=", "desc", "<desc>"],
    ["z", "sizesort", "sizesort", ""],
    ["x:", "xml=", "xml", "<output dir>"],
    ["f:", "filter=", "filter", "<min part size>"],
    ["F:", "filtergenes=", "filtergenes", "<filtered genes>"],
    ["S:", "smap=", "smap", "<species map>"],
    ["P:", "paths=", "paths", "<data paths>",
        {"default": "."}]
    ]


param = util.parseOptions(sys.argv, options, quit=True)



def componentCompositions(order, comps, gene2species):
    compositions = util.Dict(default=0)
    for comp in comps:
        counts = util.Dict(util.histDict(map(gene2species, comp)), default=0)
        key = []
        for genome in order:
            key.append(counts[genome])
        compositions[tuple(key)] += 1
    return compositions


def cmpgenes(a, b):
    gene1 = genes[a]
    gene2 = genes[b]
    
    if gene1.chrom.genome != gene2.chrom.genome:
        return cmp(gene1.chrom.genome.name, gene2.chrom.genome.name)
    elif gene1.chrom != gene2.chrom:
        return cmp(gene1.chrom.name, gene2.chrom.name)
    else:
        return cmp(gene1.start, gene2.start)



def text2xml(text):
    text = text.replace("&", "&amp;")
    text = text.replace("<", "&lt;")
    text = text.replace(">", "&gt;")    
    return text


def writeMain(out):
    print >>out, """<?xml version="1.0" encoding="ISO-8859-1"?>
    <?xml-stylesheet type="text/xsl" href="partstats.xsl"?>"""
    
    # print summary    
    print >>out, "<partstats>"
    print >>out, "<summary>"
    
    print >>out, """
        <nitems>%d</nitems>
        <nparts>%d</nparts>
        <largestpart>%d</largestpart>
        <smallestpart>%d</smallestpart>
    """ % (sum(map(len, parts)), len(parts), max(sizes), min(sizes))
    
    # genome summary
    for genome in genomes:
        orthgenes = []
        for part in parts:
            orthgenes.extend(part)
        ngenes = len(filter(lambda x: genes[x].chrom.genome.name == genome, 
                            orthgenes))
        print >>out, """
            <genome>
                <name>%s</name>
                <totalgenes>%d</totalgenes>
                <orthgenes>%d</orthgenes>
            </genome>
            """ % (genome, len(genomes[genome].genes), ngenes)
    
    print >>out, """
        <genome>
            <name>total</name>"
            <totalgenes>%d</totalgenes>
            <orthgenes>%d</orthgenes>    
        </genome>
        """ % (sum([len(x.genes) for x in genomes.values()]), 
               sum(map(len, parts)))
    
    # part histogram
    print >>out, "<histogram>"
    hist = util.histInt(sizes)
    items = [(i, hist[i]) for i in xrange(len(hist))]
    items = filter(lambda x: x[1] > 0, items)
    
    for item in items:
        print >>out, "<part><size>%d</size><num>%d</num></part>" % (item[0], item[1])
    
    print >>out, "</histogram>"
    
    print >>out, "</summary>"

    
    print >>out, "<parts>"
    
    # components
    ind = range(len(parts))

    if "sizesort" in param:
        ind.sort(lambda a, b: cmp(len(parts[b]), len(parts[a])))

    for i in ind:
        writePart(out, parts[i], i)

    print >>out, "</parts>"
    print >>out, "</partstats>"    



def writePart(out, part, name):
    counts = util.Dict(util.histDict(map(gene2species, part)), default=0)

    words = {} #descriptions.summary(part, desc)
    keys = words.keys()
    keys.sort(lambda a, b: cmp(words[b], words[a]))

    print >>out, "<part>"
    print >>out, """
        <id>%d</id>
        <size>%d</size>
        <desc>%s</desc>
    """ % \
        (name, len(part), text2xml(", ".join(map(lambda x: "%s(%d)" % (x, words[x]), 
            keys[:5]))))
    for genome in order:
        print >>out, "<comp><genome>%s</genome><count>%d</count></comp>" % \
            (genome, counts[genome]),

    print >>out, "<members>"
    part.sort(cmpgenes)
    for gene in part:
        print >>out, """
            <member>
                <name>%s</name>""" % gene


        print >>out, """
            <genome>%s</genome>
            <chrom>%s</chrom>
            <start>%d</start>
            <end>%d</end>
            """ % (genes[gene].chrom.genome.name,
                   genes[gene].chrom.name,
                   genes[gene].start,
                   genes[gene].end)
        if gene in desc:
            print >>out, "<desc>%s</desc>" % text2xml(desc[gene])
        print >>out, """
            </member>"""
    print >>out, "</members>"

    print >>out, "</part>"    





# setup paths
env.addEnvPaths("DATAPATH")
env.addPaths(param["paths"])

# read species map
gene2species = genomeutil.readGene2species(* map(env.findFile, param["smap"]))


parts = []
for f in param["part"]:
    if f == "-":
        for line in sys.stdin:
            parts.append(line.rstrip().split())
    else:
        parts.extend(util.readDelim(f))

if "filter" in param:
    cutoff = int(param["filter"][-1])
    parts = filter(lambda x: len(x) >= cutoff, parts)
    

if "genomes" in param:
    order = param["genomes"][-1].split(",")
else:
    order = util.cget(maps, 1)

desc = {}
if "desc" in param:
    for f in param["desc"]:
        desc.update(genomeio.readGeneDesc(env.findFile(f)))

go = {}
if "go" in param:
    for f in param["go"]:
        for line in file(f):
            words = line.split()
            go[words[0]] = words[1:]

# default options
param["members"] = []
param["summary"] = []
param["components"] = []
#conf["z"] = []


sizes = map(len, parts)
comps = componentCompositions(order, parts, gene2species).items()
comps.sort(lambda a, b: cmp(b[1], a[1]))


genes = {}
if "details" in param or "xml" in param:
    genomes = {}
    for genome in order:
        util.tic("read %s genome" % genome)
        genomes[genome] = genomeutil.Genome(genome)
        genomeio.readEnsmartGenes(genomes[genome])
        genomes[genome].autoconf()
        genes.update(genomes[genome].genes)
        util.toc()


if "xml" in param:
    outdir = param["xml"][-1] + "/"
    
    writeMain(file(outdir + "main.xml", "w"))
    
    for i in xrange(len(parts)):
        out = file(outdir + "part%d.xml" % i, "w")
        print >>out, """<?xml version="1.0" encoding="ISO-8859-1"?>
            <?xml-stylesheet type="text/xsl" href="partstats.xsl"?>"""
        
        writePart(out, parts[i], i)
        out.close()
else:

    print "number of items:      %d" % sum(map(len, parts))
    print "number of partitions: %d" % len(parts)
    print "largest partition:    %d" % max(sizes)
    print "smallest partition:   %d" % min(sizes)
    print

    if "summary" in param:
        print "histogram of partition sizes [size(many)]:"
        hist = util.histInt(sizes)
        stride = 8
        items = [(i, hist[i]) for i in xrange(len(hist))]
        items = filter(lambda x: x[1] > 0, items)
        for i in xrange(0, len(items), stride):
            for j in xrange(i, min(i+stride, len(items))):
                print "%d(%d) " % items[j],
            print 
        print

        print "histogram of partition compositions:"
        print order
        
        coveringcomps = filter(lambda comp: 
                            len(filter(lambda x: x>0, comp[0])) == len(comp[0]), 
                            comps)
        totalcovering = sum(util.cget(coveringcomps, 1))
        
        print ("(" + ", ".join(["+"] * len(order)) + ")"), totalcovering
        
        for comp in comps:
            print comp[0], comp[1]
        print


    if "components" in param:
        ind = range(len(parts))

        if "sizesort" in param:
            ind.sort(lambda a, b: cmp(len(parts[b]), len(parts[a])))

        for i in ind:
            part = parts[i]

            counts = util.Dict(util.histDict(map(gene2species, part)), default=0)

            words = {} #descriptions.summary(part, desc)
            keys = words.keys()
            keys.sort(lambda a, b: cmp(words[b], words[a]))

            print
            print "part %d: size %d [%s]" % \
                (i, len(part), ", ".join(map(lambda x: "%s(%d)" % (x, words[x]), 
                    keys[:5])))
            for genome in order:
                print "%s(%d) " % (genome, counts[genome]),
            print

            if "members" in param:
                print "-----------------------------------------"
                part.sort()
                for gene in part:
                    extra = ""
                    if gene in genes:
                        extra += fasta.chromFormat(
                            genes[gene].chrom.name,
                            genes[gene].start,
                            genes[gene].end) + " " + \
                            "(%d) " % (genes[gene].end - genes[gene].start + 1)
                            
                    if gene in go:
                        extra += str(go[gene]) + " "
                    if gene in desc:
                        extra += desc[gene]
                    print "%20s\t%s" % (gene, extra)


                print "-----------------------------------------"


