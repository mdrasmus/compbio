#!/usr/bin/env python

from rasmus import util, ensembl, genomeio, descriptions, genomeutil, fasta
import sys, os


def mkcmpgenes(genes):
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


def partDescription(config, part, desc, annotations):
    # determine description summary
    if "report_part_summary" in config:
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


def writeMain(config, out, genomes, genes, parts, desc, 
        syntenyLookup, annotations):

    # generate genome compositions
    sizes = map(len, parts)
    comps = genomeutil.componentCompositions(config["genomes"], parts).items()
    comps.sort(lambda a, b: cmp(b[1], a[1]))
        
        
    print >>out, """<?xml version="1.0" encoding="ISO-8859-1"?>
    <?xml-stylesheet type="text/xsl" href="../wwwbase/partstats.xsl"?>"""
    
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
    for genome in config["genomes"]:
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

    
    ind.sort(lambda a, b: cmp(len(parts[b]), len(parts[a])))

    for i in ind:
        writePart(config, out, parts[i], i, genes, desc, syntenyLookup, annotations)

    print >>out, "</parts>"
    print >>out, "</partstats>"    



    

def writePart(config, out, part, name, genes, desc, syntenyLookup, annotations):
    if len(part) == 0:
        return

    counts = genomeutil.genomeComposition(config["genomes"], part)
    description = partDescription(config, part, desc, annotations)
    

    print >>out, "<part>"
    print >>out, """
        <id>%d</id>
        <size>%d</size>
        <desc>%s</desc>
    """ % \
        (name, len(part), text2xml(description))
    for genome in config["genomes"]:
        print >>out, "<comp><genome>%s</genome><count>%d</count></comp>" % \
            (genome, counts[genome]),

    
    # synteny visualization files
    visfiles = util.sort(util.unique(util.subdict(syntenyLookup, part).values()))
    for visfile in visfiles:
        print >>out, """
            <synteny>%s</synteny>
        """ % visfile
        

    if "align_file_func" in config:
        print >>out, """
            <align>%s</align>
        """ % config["align_file_func"](name)

    if "tree_file_func" in config:
        print >>out, """
            <tree>%s</tree>
        """ % config["tree_file_func"](name)
    

    print >>out, "<members>"
    part.sort(mkcmpgenes(genes))
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


def writeAll(config, matching, parts, desc):  
    # load genome information
    genes = matching.getGenes()
    
    # load annotations
    annotations = genomeio.readGenes(config["genomes"])
    for gene in genes:
        if gene not in annotations:
            annotations[gene] = {}
        if "gene" not in annotations[gene]:
            annotations[gene]["gene"] = ""

    # make sure there aren't any unknown genes
    for part in parts:
        for gene in part:
            assert gene in genes
    
    # load synteny
    syntenyLookup = {}
    if "vissynteny_index_file" in config:
        for line in file(config["vissynteny_index_file"]):
            words = line.split("\t")
            filename = words[0]
            for gene in words[1:]:
                syntenyLookup[gene] = filename
    
        
    # write main page and convert to html
    out = file(config["report_dir"] + "index.xml", "w")
    writeMain(config, out, matching.genomes, genes, parts, desc, 
        syntenyLookup, annotations)
    out.close()
    os.system("xsltproc %spartstats.xsl %sindex.xml > %sindex.html" % 
        (config["wwwbase_dir"], config["report_dir"], config["report_dir"]))
    
    
    # write out part pages
    for i in xrange(len(parts)):
        out = file(config["report_dir"] + "part%d.xml" % i, "w")
        print >>out, """<?xml version="1.0" encoding="ISO-8859-1"?>
            <?xml-stylesheet type="text/xsl" href="../wwwbase/partstats.xsl"?>"""

        writePart(config, out, parts[i], i, genes, desc, syntenyLookup, annotations)
        out.close()
