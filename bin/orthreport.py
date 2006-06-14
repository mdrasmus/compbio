#!/usr/bin/env python

from rasmus import util, ensembl, genomeio, descriptions, genomeutil, fasta
import sys



options = [
    ("p:", "part=", "part", "[ -p <part file>] [--part=<part file>]"),
    ("g:", "genomes=", "genomes", "[ -g <genome>] [--genome=<genome>]"),
    ("d:", "desc=", "desc", "[ -d <desc>] [ --desc <desc> ]"),
    ("z", "sizesort", "sizesort", "[-z] [-sizesort]"),
    ["f:", "filter=", "filter", "AUTO<min part size>"],
    ["o:", "outdir=", "outdir", "AUTO<output directory>"]
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)



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


order = param["genomes"][-1].split(",")


desc = {}
if "desc" in param:
    for f in param["desc"]:
        desc.update(genomeio.readGeneDesc(f))


sizes = map(len, parts)
comps = ensembl.componentCompositions(order, parts).items()
comps.sort(lambda a, b: cmp(b[1], a[1]))


genes = {}
genomes = {}
for genome in order:
    util.tic("read %s genome" % genome)
    genomes[genome] = genomeutil.Genome(genome)
    genomeio.readEnsmartGenes(genomes[genome])
    genomes[genome].autoconf()
    genes.update(genomes[genome].genes)
    util.toc()



def text2xml(text):
    text = text.replace("&", "&amp;")
    text = text.replace("<", "&lt;")
    text = text.replace(">", "&gt;")    
    return text



def writeMain(out):
    print >>out, """
        <html>
        <head>
            <title>Orthologs</title>
        </head>
        <body>
        """


    # print summary    
    print >>out, """
        <table border='1'>
        <tr><td># items:</td><td>%d</td></tr>
        <tr><td># parts:</td><td>%d</td></tr>
        <tr><td>largest part</td><td>%d</td></tr>
        <tr><td>smallest part</td><td>%d</td></tr>
        </table>
    """ % (sum(map(len, parts)), len(parts), max(sizes), min(sizes))

    # genome summary
    print >>out, """
        <table border='1'>
            <tr>
                <td>genome</td>
                <td>matched genes</td>
                <td>total genes</td>
            </tr>
        """
    for genome in genomes:
        orthgenes = []
        for part in parts:
            orthgenes.extend(part)
        ngenes = len(filter(lambda x: genes[x].chrom.genome.name == genome, 
                            orthgenes))
        print >>out, """
            <tr>
                <td>%s</td>
                <td>%d</td>
                <td>%d (%2f)</td>
            </genome>
            """ % (genome, len(genomes[genome].genes), ngenes,
                      100.0 * len(genomes[genome].genes) / ngenes)


    totalmatch = sum([len(x.genes) for x in genomes.values()])
    totalgenes = sum(map(len, parts))

    print >>out, """
        <tr>
            <td>total</td>"
            <td>%d</td>
            <td>%d (%2f)</td>
        </genome>
        </table>
        """ % (totalmatch, totalgenes, 100.0 * totalmatch / totalgenes)


    # part histogram
    print >>out, """
        <table border='1'>
            <tr><td>part size</td><td>many</td></tr>
        """
    hist = util.histInt(sizes)
    items = [(i, hist[i]) for i in xrange(len(hist))]
    items = filter(lambda x: x[1] > 0, items)

    for item in items:
        print >>out, "<tr><td>%d</td><td>%d</td></tr>" % (item[0], item[1])

    print >>out, "</table>"


outdir = param["outdir"][-1] + "/"

writeMain(file(outdir + "main.html", "w"))


if False:
    # components
    ind = range(len(parts))

    if "sizesort" in param:
        ind.sort(lambda a, b: cmp(len(parts[b]), len(parts[a])))

    for i in ind:
        part = parts[i]

        counts = ensembl.genomeComposition(order, part)

        words = descriptions.summary(part, desc)
        keys = words.keys()
        keys.sort(lambda a, b: cmp(words[b], words[a]))

        print "<part>"
        print """
            <id>%d</id>
            <size>%d</size>
            <desc>%s</desc>
        """ % \
            (i, len(part), text2xml(", ".join(map(lambda x: "%s(%d)" % (x, words[x]), 
                keys[:5]))))
        for genome in order:
            print "<comp><genome>%s</genome><count>%d</count></comp>" % \
                (genome, counts[genome]),

        print "<members>"
        part.sort()
        for gene in part:
            print """
                <member>
                    <name>%s</name>""" % gene

            if gene in genes:
                print "<pos>%s</pos>" % fasta.chromFormat(
                    genes[gene].chrom.name,
                    genes[gene].start,
                    genes[gene].end)
            if gene in desc:
                print "<desc>%s</desc>" % text2xml(desc[gene])
            print """
                </member>"""
        print "</members>"

        print "</part>"

    print "</parts>"
    print "</partstats>"
