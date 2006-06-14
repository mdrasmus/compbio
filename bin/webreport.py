#!/usr/bin/env python

from rasmus import util, ensembl, genomeio, descriptions, genomeutil, fasta
import sys, os


#
# Parse arguments 
#
options = [
    ["p:", "part=", "part", "AUTO<part file>"],
    ["g:", "genomes=", "genomes", "AUTO<genomes>"],
    ["d:", "desc=", "desc", "AUTO<desc file>"],
    ["o:", "outdir=", "outdir", "AUTO<output dir>"],
    ["s", "summary", "summary", "AUTO"],
    ["a", "anno", "anno", "AUTO"],
    ["W:", "chromWindows=", "chromWindows", "AUTO<synteny vis index>"]
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)
    

assert "part" in param, "No partition file given"
assert "genomes" in param, "No genomes file given"
assert "outdir" in param, "No output directory given"


#
# Parse arguments
#
parts = []
for f in param["part"]:
    if f == "-":
        for line in sys.stdin:
            parts.append(line.rstrip().split())
    else:
        parts.extend(util.readDelim(f))


# load genome sorting order
order = param["genomes"][-1].split(",")


# load descriptions
desc = {}
if "desc" in param:
    for f in param["desc"]:
        desc.update(genomeio.readGeneDesc(f))


# generate genome compositions
sizes = map(len, parts)
comps = ensembl.componentCompositions(order, parts).items()
comps.sort(lambda a, b: cmp(b[1], a[1]))


# load genome information
genes = {}
genomes = {}
for genome in order:
    util.tic("read %s genome" % genome)
    genomes[genome] = genomeutil.Genome(genome)
    genomeio.readEnsmartGenes(genomes[genome])
    genomes[genome].autoconf()
    genes.update(genomes[genome].genes)
    util.toc()


# load annotations
annotations = genomeio.readGenes(order)
for row in annotations.values():
    if "gene" not in row:
        row["gene"] = ""


# load synteny
syntenyLookup = {}
if "chromWindows" in param:
    for line in file(param["chromWindows"][-1]):
        words = line.split("\t")
        filename = words[0]
        for gene in words[1:]:
            syntenyLookup[gene] = filename
    
    


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


def partDescription(part, desc):
    # determine description summary
    if "summary" in param:
        # use a summary of all descriptions
        words = descriptions.summary(part, desc)
        keys = words.keys()
        keys.sort(lambda a, b: cmp(words[b], words[a]))    
        description = ", ".join(
                              map(lambda x: "%s(%d)" % (x, words[x]), keys[:5]))
    else:
        # try to choose one description as representive
        common = [annotations[x]["gene"] for x in part]
        if max(map(len, common)) > 0:
            # use a gene with a common name (longest)
            gene =  part[util.argmaxfunc(len, common)]
            description = desc.setdefault(gene, "")
        else:
            # use longest description
            i = util.argmaxfunc(len, util.subdict(desc, part).values())
            if i != -1:
                description = util.subdict(desc, part).values()[i]
            else:
                description = ""    
    
    return description


def writeMain(out):
    print >>out, """<?xml version="1.0" encoding="ISO-8859-1"?>
    <?xml-stylesheet type="text/xsl" href="../partstats.xsl"?>"""
    
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
    if len(part) == 0:
        return

    counts = ensembl.genomeComposition(order, part)
    description = partDescription(part, desc)
    

    print >>out, "<part>"
    print >>out, """
        <id>%d</id>
        <size>%d</size>
        <desc>%s</desc>
    """ % \
        (name, len(part), text2xml(description))
    for genome in order:
        print >>out, "<comp><genome>%s</genome><count>%d</count></comp>" % \
            (genome, counts[genome]),

    
    # synteny visualization files
    visfiles = util.sort(util.unique(util.subdict(syntenyLookup, part).values()))
    for visfile in visfiles:
        print >>out, """
            <synteny>%s</synteny>
        """ % visfile
        
    
    

    print >>out, "<members>"
    part.sort(cmpgenes)
    for gene in part:
        # name
        print >>out, """
            <member>
                <name>%s</name>
        """ % gene
        
        # common name
        if annotations[gene]["gene"] != "":
            print >>out, """
                <commonname>%s</commonname>
            """ % annotations[gene]["gene"]
        

        # location info
        print >>out, """
            <genome>%s</genome>
            <chrom>%s</chrom>
            <start>%d</start>
            <end>%d</end>
            """ % (genes[gene].chrom.genome.name,
                   genes[gene].chrom.name,
                   genes[gene].start,
                   genes[gene].end)
        
        # description
        if gene in desc:
            print >>out, "<desc>%s</desc>" % text2xml(desc[gene])
         
        # Ensembl link
        if gene.startswith("ENS"):
            print >>out, """
                <link>http://www.ensembl.org/Homo_sapiens/textview?species=all&amp;idx=Gene&amp;q=%s</link>
            """ % gene
        
        # synteny visualization
        if gene in syntenyLookup:
            print >>out, """
                <synteny>%s</synteny>
            """ % syntenyLookup[gene]
        
        print >>out, """
            </member>"""
    print >>out, "</members>"

    print >>out, "</part>"    


def main():
    # determine paths
    outdir = param["outdir"][-1].rstrip("/")
    basedir = os.path.split(outdir)[0] + "/"
    outdir += "/"

    # write main page and convert to html
    writeMain(file(outdir + "index.xml", "w"))
    os.system("xsltproc %spartstats.xsl %sindex.xml > %sindex.html" % 
        (basedir, outdir, outdir))

    # write out part pages
    for i in xrange(len(parts)):
        out = file(outdir + "part%d.xml" % i, "w")
        print >>out, """<?xml version="1.0" encoding="ISO-8859-1"?>
            <?xml-stylesheet type="text/xsl" href="../partstats.xsl"?>"""

        writePart(out, parts[i], i)
        out.close()




main()
