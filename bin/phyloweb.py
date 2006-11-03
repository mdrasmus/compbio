#!/usr/bin/env python

import os
import sys

from rasmus import util, env, genomeutil, treelib, progress, phylip, phyloutil
from rasmus.vis import treevis
from rasmus import svg, tablelib, fasta



options = [
    ["n:", "numtrees=", "numtrees", "<max number of trees>",
     {"parser": int,
      "single": True}],
    ["r", "resume", "resume", "",
     {"single": True}],
    
    "Options",
    ["", "separateByChrom=", "separateByChrom", "",
     {"default": True,
      "parser": util.str2bool}],
    ["", "useTopNames=", "useTopNames", "",
     {"default": True,
      "parser": util.str2bool}],
    ["", "dataTwoTiers=", "dataTwoTiers", "<bool>",
     {"single": True,
      "default": True,
      "parser": bool}],
    ["", "useGeneNames=", "useGeneNames", "<bool>",
     {"single": True,
      "default": True,
      "parser": bool}],
    ["", "useRfError=", "useRfError", "<bool>",
     {"single": True,
      "default": False,
      "parser": bool}],
    
    "Extentions",
    ["", "ntAlign=", "ntAlign", "<nt alignment extension>",
     {"default": ".nt.align"}],
    ["", "pepAlign=", "pepAlign", "<pep alignment extension>",
     {"default": ".pep.align"}],
    
    ["", "treeVisFormat=", "treeVisFormat", "<tree visualization format>",
     {"default": "svg"}],
    
    ["", "correctNtTree=", "correctNtTree", "<correct nt tree extension>",
     {"default": ".nt.tree"}],
    ["", "correctPepTree=", "correctPepTree", "<correct pep tree extension>",
     {"default": ".pep.tree"}],
    ["", "resultTree", "resultTree", "<result tree extension>",
     {"default": ".tree"}],     
    
    ["", "alignVisDir=", "alignVisDir", "<directory for alignment visulizations>",
     {"default": "alnvis"}],
    ["", "treeVisDir=", "treeVisDir", "<directory for tree visulizations>",
     {"default": "treevis"}],
    
    "Misc",
    ["", "topColors=", "topColors", "<topology colors>",
     {"default": ["#ccf", "#f40", "#f80", "#fb0", "#ff0", 
                          "#bf4", "#8f8", "#4fb", "#0ff", 
                          "#0bf", "#08f", #"#04f", #"#00f",
                          "#888"]}]
     
]



def color2html(r, g, b, a=1):
    txt = "#"
    digits = "0123456789abcdef"
        
    for val in [r, g, b]:
        txt += digits[int(val * 15)]
    
    return txt


def getGeneName(alignfile, mainspecies, gene2species):
    if conf["useGeneNames"]:
        for line in file(alignfile):
            if line.startswith(">"):
                gene = line[1:].rstrip()
                if gene2species(gene) == mainspecies:
                    return gene
    return None

######################################
# file and url naming conventions
#

def getAlignUrl(conf, dataurl, treename, seqtype = "nt"):
    if seqtype == "nt":
        ext = conf["ntAlign"]
    else:
        ext = conf["pepAlign"]
    
    if conf["dataTwoTiers"]:
        return os.path.join(dataurl, treename, treename + ext)
    else:
        return os.path.join(dataurl, treename + ext)
    

def getAlignFile(conf, datadir, treename, seqtype = "nt"):
    if seqtype == "nt":
        ext = conf["ntAlign"]
    else:
        ext = conf["pepAlign"]

    if conf["dataTwoTiers"]:
        return os.path.join(datadir, treename, treename + ext)
    else:
        return os.path.join(datadir, treename + ext)


def getAlignVisFile(conf, treename, seqtype = "nt"):
    return os.path.join(conf["webdir"], conf["alignVisDir"], 
                        treename + "." + seqtype + ".txt")

def getAlignVisUrl(conf, treename, seqtype = "nt"):
    return os.path.join(conf["weburl"], conf["alignVisDir"], 
                        treename + "." + seqtype + ".txt")


def getCorrectTreeFile(conf, datadir, treename, seqtype = "nt"):
    if seqtype == "nt":
        ext = conf["correctNtTree"]
    else:
        ext = conf["correctPepTree"]

    if conf["dataTwoTiers"]:
        return os.path.join(datadir, treename, treename + ext)
    else:
        return os.path.join(datadir, treename + ext)

def getCorrectTreeUrl(conf, dataurl, treename, seqtype = "nt"):
    if seqtype == "nt":
        ext = conf["correctNtTree"]
    else:
        ext = conf["correctPepTree"]

    if conf["dataTwoTiers"]:
        return os.path.join(dataurl, treename, treename + ext)
    else:
        return os.path.join(dataurl, treename + ext)


def getResultTreeFile(conf, resultdir, treename):
    return os.path.join(resultdir, treename + conf["resultTree"])

def getResultTreeUrl(conf, resulturl, treename):
    return os.path.join(resulturl, treename + conf["resultTree"])


def getSyntenyUrl(conf, filename):
    return conf['syntenyUrl'] + "/" + filename


def getTreeVisFile(conf, treename, kind):
    return os.path.join(conf["webdir"], conf["treeVisDir"], 
                        treename + "." + kind + "." + conf["treeVisFormat"])

def getTreeVisUrl(conf, treename, kind):
    return os.path.join(conf["weburl"], conf["treeVisDir"], 
                        treename + "." + kind + "." + conf["treeVisFormat"])




#####################################
# Web generation
#

def readSynteny(filename, genes):
    tab = tablelib.readTable(filename)
    
    geneLookup = {}
    
    for row in tab:
        rowGenes = row['genes'].split(',')
        
        for gene in rowGenes:
            if gene not in genes:
                continue
            
            dist = min(genes[gene]['start'] - row['start'],
                       genes[gene]['end'] - row['start'])
            if gene not in geneLookup or \
               dist > geneLookup[gene][1]:
                geneLookup[gene] = [row['file'], dist]
    
    # map directly to filename
    geneLookup = util.mapdict(geneLookup, valfunc=lambda x: x[0])
    
    return geneLookup


def readGeneCoords(filename):
    
    genes = tablelib.Table(headers=["gene", "chrom", "start", "end", "strand"])
    
    for line in file(filename):
        tokens = line.split()
        
        genes.add(gene=tokens[0],
                  chrom=tokens[1],
                  start=int(tokens[2]),
                  end=int(tokens[3]),
                  strand=int(tokens[4]))
    
    return genes


def readCommonNames(filename):
    common = util.Dict(default="-")
    
    for line in file(filename):
        tokens = line.rstrip().split("\t")   
        if len(tokens) == 2:
            common[tokens[0]] = tokens[1]
    
    return common


def readResults(filename):
    infile = file(filename)
    results = {}
    
    for line in infile:
        if len(line.rstrip()) == 0:
            break
    infile.next()
    
    for line in infile:
        tokens = line.split()
        
        if len(tokens) == 3 and tokens[0].isdigit():
            results[tokens[0]] = [util.str2bool(tokens[1]), float(tokens[2])]
        else:
            break
    
    return results
            


def getTreeNames(conf, datadir):
    if conf["dataTwoTiers"]:
        treenames = os.listdir(datadir)
    else:
        files = os.listdir(datadir)
        nums = []
        
        for f in files:
            if f[0].isdigit():
                i = f.index(".")
                nums.append(f[:i])
        treenames = sorted(util.unique(nums))
    
    if "numtrees" in conf:
        numtrees = conf["numtrees"]
    else:
        numtrees = len(treenames)
    
    if "evalSubset" in conf:
        treenames = filter(lambda x: x in conf["evalSubset"], treenames)
    
    treenames = treenames[:numtrees]
    
    return treenames


def generateGeneTreeTables(conf, datadir, resultdirs):
    genes = readGeneCoords(env.findFile(conf["mainspeciesCoords"]))
    
    chroms = genes.groupby(lambda x: x["chrom"])
    chromnames = sorted(chroms.keys())
    
    treeStats = conf["treeStatsLookup"]
    topNames = conf["topNamesLookup"]
    
    
    # main page
    out = file(os.path.join(conf["webdir"], "index.html"), "w")
    
    print >>out, "<html><head></head><body> <ul>"

    
    for chromname in chromnames:
        genelookup = chroms[chromname].lookup('gene')
        
        # get treenames from genes
        treenames = []        
        for row in treeStats.itervalues():
            if row['gene'] in genelookup:
                treenames.append(row['treeid'])        
        treenames.sort(key=lambda x: genelookup[treeStats[x]['gene']]['start'])
        treenames = map(str, treenames)
        
        
        if len(treenames) == 0:
            continue
        
        
        print >>out, "<li><a href='%s.html'>%s</a> %d trees" % \
            (chromname, chromname, len(treenames))
        
        generateGeneTreeTable(conf, chromname + ".html", treenames, 
                              datadir, resultdirs)
    
        
    print >>out, "</ul></body></html>"
    
    out.close()
    

def generateGeneTreeTable(conf, filename, treenames, datadir, resultdirs):
    util.tic("generate gene tree table '%s'" % filename)
    
    progs = util.cget(resultdirs, 0)    
    dataurl = conf["dataurl"]
    
    if treenames == None:
        treenames = getTreeNames(conf, datadir)
    
    
    treeStats = conf["treeStatsLookup"]
    topNames = conf["topNamesLookup"]
    genes = readGeneCoords(env.findFile(conf["mainspeciesCoords"])).lookup("gene")
    
    if "syntenyFile" in conf:
        syntenyIndex = readSynteny(conf["syntenyFile"], genes)
    else:
        syntenyIndex = None
    
    
    # create output
    os.system("mkdir -p %s" % conf["webdir"])
    out = file(os.path.join(conf["webdir"], filename), "w")    

    writeHeader(out)

    # header    
    print >>out, """<table cellspacing='1' cellpadding='5'>
                        <tr>
                            <td><b>Tree ID</b></td>
                            <td><b>Chrom pos</b></td>
                            <td><b>Gene</b></td> 
                            <td><b>Common</b></td>
                            <td><b>Gene len</b></td>
                            <td><b>Tree len</b></td>
                            <td><b>Links</b></td>"""
    for prog in progs:
        if "progNames" in conf:
            print >>out, "<td><b>%s</b></td>" % conf["progNames"][prog]
        else:
            print >>out, "<td><b>%s</b></td>" % prog
    print >>out, "</tr>"
    
    colormap = util.ColorMap([[0, (1, 1,1)],
                              [.3, (1, 0, 0)]])
    
    
    # table
    progbar = progress.ProgressBar(len(treenames))
    for treename in treenames:
        progbar.update()
        treenum = int(treename)
        
        assert treenum in treeStats, treenum
        
        gene = treeStats[treenum]['gene']        
        
        # get align files
        if os.path.exists(getAlignFile(conf, datadir, treename, "nt")):
            alnnturl = "http://compbio.mit.edu/spidir/browser/vis.cgi?vis=align&data=%s" % getAlignUrl(conf, dataurl, treename, "nt")
            alnntlink = "<a href='%s'>An</a> " % alnnturl
        else:
            alnntlink = ""
        
        if os.path.exists(getAlignFile(conf, datadir, treename, "pep")):
            alnpepurl = "http://compbio.mit.edu/spidir/browser/vis.cgi?vis=align&data=%s" % getAlignUrl(conf, dataurl, treename, "pep")
            alnpeplink = "<a href='%s'>Ap</a> " % alnpepurl
        else:
            alnpeplink = ""

        
        
        # get tree files
        treefile = getCorrectTreeUrl(conf, dataurl, treename)
        treeurl = "http://compbio.mit.edu/spidir/browser/vis.cgi?vis=tree&data=%s" % treefile
        
        row = treeStats[treenum]
        
        if gene == "":
            genename = "-"
            commonName = "-"
            pos = "-"
        else:
            genename = gene
            commonName = row["common"]
            chromName = genes[gene]['chrom']
            chromName = chromName.replace("chromosome", "chr")
            if chromName[0] != "c":
                chromName = "chr" + chromName
            pos = chromName + ":" + util.int2pretty(genes[gene]['start'])
        
        
        # gene name and tree ID
        if "geneOutLink" in conf and gene != "":
            genetext = "<a href='%s'>%s</a>" % \
                (conf['geneOutLink'] % genename, genename)
        else:
            genetext = genename
        
        if row["treelen"] == 0:
            treelen = "-"
        else:
            treelen = "%.3f" % row["treelen"]
        
        print >>out, "<tr><td><a href='%s'>%s</a></td>  <td>%s</td> <td>%s</td> " % \
            (os.path.join(dataurl, treename), treename, 
             pos, genetext)
        
        print >>out, "<td>%s</td><td>%d</td><td>%s</td>" % \
            (commonName, row["alignlen"], treelen)
        
        # synteny
        if syntenyIndex and gene in syntenyIndex:
            syntenyfile = getSyntenyUrl(conf, syntenyIndex[gene])
            syntenyurl = "http://compbio.mit.edu/spidir/browser/vis.cgi?vis=synteny&data=%s" % syntenyfile
            syntenylink = "<a href='%s'>S</a> " % syntenyurl
        else:
            syntenylink = ""
        
        
        # links
        print >>out, ("<td><nobr>" + \
                          alnntlink +
                          alnpeplink + \
                          syntenylink + \
                          "<a href='%s'>T</a></nobr></td>") % \
            (treeurl)
        
        
        for prog, resultdir, resulturl in resultdirs:
            treefile = getResultTreeUrl(conf, resulturl, treename)
            resulttree = "http://compbio.mit.edu/spidir/browser/vis.cgi?vis=tree&data=%s" % treefile
            
            if conf["useTopNames"]:
                top = treeStats[treenum][prog]
                topname = topNames[top]["name"]
                topnum = int(topname[1:]) - 1
            
                if topnum < len(conf["topColors"]) - 1:
                    color = conf["topColors"][topnum]
                else:
                    color = conf["topColors"][-1]
            else:
                if treeStats[treenum][prog + "_rferror"] == 0:
                    color = "#ccf"
                    if conf["useRfError"]:
                        topname = "0"
                    else:
                        topname = "True"
                else:
                    if conf["useRfError"]:
                        color = color2html(* colormap.get(
                                treeStats[treenum][prog + "_rferror"]))
                        topname = "%.3f" % treeStats[treenum][prog + "_rferror"]
                    else:
                        color =   "#f55"
                        topname = "False"
            
            print >>out, "<td style='background-color: %s'><a href='%s'>%s</a></td>" % \
                (color, resulttree, topname)
    
    
    # footer
    print >>out, "</table>"
    

    writeFooter(out)
    
    out.close()
    util.toc()


def generateTopologyHistograms(conf, datadir, resultdirs):
    treeStats = tablelib.Table(conf["treeStatsLookup"].values())
    topNames = conf["topNamesLookup"]
    dataurl = conf["dataurl"]
    
    numtops = 20
    
    out = file(os.path.join(conf["webdir"], "histogram.html"), "w")
    
    writeHeader(out)


    # table header    
    print >>out, """
        <table><tr>
        """

    examples = {}
    for row in treeStats:
        for prog, resultdir, resulturl in resultdirs:
            topname = topNames[row[prog]]['name']
            examples[topname] = [resulturl, str(row['treeid'])]

    tophists = {}
    for prog, resultdir, resulturl in resultdirs:
        tops = treeStats.cget(prog)
        tophists[prog] = tablelib.histTable(tops)
        
        if "progNames" in conf:
            progname = conf["progNames"][prog]
        else:
            progname = prog
        
        print >>out, "<td><b>%s</b></td>" % progname
    
    print >>out, "</tr>"
    
    # table data
    for i in range(numtops):
        print >>out, "<tr>"
        
        
        for prog, resultdir, resulturl in resultdirs:
            tab = tophists[prog]
            
            if i < len(tab):
                topname = topNames[tab[i]['item']]['name']
                topnum = int(topname[1:]) - 1

                if topnum < len(conf["topColors"]) - 1:
                    color = conf["topColors"][topnum]
                else:
                    color = conf["topColors"][-1]
                
                resulturl, treename = examples[topname]
                treefile = getResultTreeUrl(conf, resulturl, treename)
                treeurl = "http://compbio.mit.edu/spidir/browser/vis.cgi?vis=tree&data=%s" % treefile
        
                
                print >>out, "<td style='background-color: %s'><a href='%s'>%s</a> (%.1f%%)</td>" % (
                    color,
                    treeurl,
                    topname, 
                    tophists[prog][i]['percent'] * 100)
            else:
                print >>out, "<td></td>"
        
        print >>out, "</tr>"
    
    writeFooter(out)
    
        
    

################################
# Topology stats
#

def generateTopologyTable(conf, datadir, resultdirs):
    treenames = getTreeNames(conf, datadir)
    
    
    # skip is stats table already exists
    if os.path.exists(conf["treeStats"]) and conf["resume"]:
        treestats = tablelib.readTable(conf["treeStats"])
        if len(treestats) >= len(treenames):
            return
    
    
    util.tic("generate tree stats")
    
    if "commonNames" in conf:
        common = readCommonNames(conf["commonNames"])
    else:
        common = util.Dict(default="-")
    
    
    # load up correctness results
    results = {}
    for prog, resultdir, resulturl in resultdirs:
        results[prog] = readResults(os.path.join(resultdir, "results"))
        
    
    # init table
    headers = ["treeid", "gene", "common", "alignlen", "treelen"]
    for prog in util.cget(resultdirs, 0):
        headers.append(prog)
        headers.append(prog + "_rferror")
        
    
    treestats = tablelib.Table(headers=headers)
    counts = util.Dict(default=0)
    
    # populate table
    progbar = progress.ProgressBar(len(treenames))
    for treename in treenames:
        progbar.update()
        
        alnfile = getAlignFile(conf, datadir, treename)
        gene = getGeneName(alnfile, conf["mainspecies"], conf["gene2species"])
        
        if gene == None:
            gene = ""
            commonName = ""
        else:
            commonName = common[gene]
        
        row = {"treeid": int(treename), "gene": gene, "common": commonName}
        
        tree = treelib.readTree(getCorrectTreeFile(conf, datadir, treename))
        row["treelen"] = sum(x.dist for x in tree.nodes.values())
        
        aln = fasta.readFasta(getAlignFile(conf, datadir, treename, "nt"))
        row["alignlen"] = aln.alignlen()
        
        for prog, resultdir, resulturl in resultdirs:
            treefile = getResultTreeFile(conf, resultdir, treename)
            tree = treelib.readTree(treefile)
            thash = phyloutil.hashTree(tree, conf["gene2species"])
            
            row[prog] = thash
            row[prog + "_rferror"] = results[prog][treename][1]
            counts[thash] += 1
        treestats.append(row)
    
    # save table
    treestats.write(conf["treeStats"])
    
    # save treenames
    items = counts.items()
    items.sort(key=lambda x: x[1], reverse=True)
    topnames = tablelib.Table(headers=["name", "tree", "count"])
    
    for i, (tree, count) in enumerate(items):
        topnames.add(name="T%d" % (i+1), tree=tree, count=count)
    
    topnames.write(conf["treeNames"])
    
    util.toc()
    

def writeHeader(out):
    
    if "geneTableHeader" in conf:
        print >>out, conf["geneTableHeader"]
    else:
        print >>out, """
        <html>
        <head>
            <style>
                a {
                    color: black;
                }

                a:visited {
                    color: #00c;
                }

                body {
                    font-family: helvetica;
                    font-size: 10pt;
                }

                td {
                    font-size: 10pt;
                }
            </style>
        </head>
        <body>

        <a href="../index.html"> [home]</a>

        <h3>Key</h3>
        <ul>
            <li>An - nucleotide alignment
            <li>Ap - peptide alignment
            <li>S - Synteny (not avail)
            <li>T - Correct tree
            <li>T# - Topology number (T1 = correct topology)
        </ul>
        """

def writeFooter(out):
    print >>out, """</body></html>"""






def main(conf):
    # read data
    env.addEnvPaths("DATAPATH")
    conf['gene2species'] = genomeutil.readGene2species(env.findFile(conf['smap']))
    conf['stree'] = treelib.readTree(env.findFile(conf['stree']))
    
    
    generateTopologyTable(conf, conf["datadir"], conf["resultdirs"])    
    
    # load tree stats and names
    conf["treeStatsLookup"] = tablelib.readTable(conf["treeStats"]).lookup("treeid")
    conf["topNamesLookup"] = tablelib.readTable(conf["treeNames"]).lookup("tree")
    
    if conf["useTopNames"]:
        generateTopologyHistograms(conf, conf["datadir"], conf["resultdirs"])
    
    if conf["separateByChrom"]:
        generateGeneTreeTables(conf, conf["datadir"], conf["resultdirs"])
    else:
        generateGeneTreeTable(conf, "index.html", None, 
                              conf["datadir"], conf["resultdirs"])


if __name__ == "__main__":
    conf = util.parseOptions(sys.argv, options, quit=True)
    main(conf)
