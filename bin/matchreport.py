#!/usr/bin/env python

from rasmus import util, ensembl, genomeio, descriptions, genomeutil, fasta
import sys, os



options = [
    ["g:", "genomes=", "genomes", "AUTO<genome1>,<genome2>,..."],
    ["c:", "cols=", "cols", "AUTO<gene1 col>,<gene2 col>,<score col>"],
    ["o:", "outdir=", "outdir", "AUTO<output directory>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)

cols = map(int, param["cols"][-1].split(","))

genomes = {}
genes = {}
for genome in param["genomes"][-1].split(","):
    util.tic("read %s genome" % genome)
    genomes[genome] = genomeutil.Genome(genome)
    genomeio.readEnsmartGenes(genomes[genome])
    genomes[genome].autoconf()
    genes.update(genomes[genome].genes)
    util.toc()


def bgcolor(score):
    num = 255 - int(255 * min(score,1000) / 1000)
    text = hex(num).replace("0x","")
    text = ("0" * (2 - len(text))) + text
    return "#ff" + text + text


def cmpgenes(a, b):
    gene1 = genes[a]
    gene2 = genes[b]
    
    if gene1.chrom.genome != gene2.chrom.genome:
        return cmp(gene1.chrom.genome.name, gene2.chrom.genome.name)
    elif gene1.chrom != gene2.chrom:
        return cmp(gene1.chrom.name, gene2.chrom.name)
    else:
        return cmp(gene1.start, gene2.start)



# read in matches
for fn in rest:
    util.log(fn)
    
    mat = util.Dict(2, 0)
    
    for line in file(fn):
        tokens = line.rstrip().split()
        gene1 = tokens[cols[0]]
        gene2 = tokens[cols[1]]
        score = float(tokens[cols[2]])
        
        mat[gene1][gene2] = score
        mat[gene2][gene1] = score

    keys = mat.keys()
    keys.sort(cmpgenes)
    
    name = os.path.basename(fn).replace(".match", "")
    outfile = param["outdir"][-1] + "/match" + \
              os.path.basename(fn).replace(".match", ".html")
    out = file(outfile, "w")
    print >>out, """
    <html>
    <head><title>%s</title></head>
    <body>
    <table border="1" style="font-size: 10">
    <tr>
        <td></td>
    """ % (name)
    
    for i in xrange(len(keys)):
        print >>out, "<td>%d</td>" % (i+1)
    print >>out, "</tr>"
    
    for i in xrange(len(keys)):
        print >>out, "<tr><td>%d. %s</td>" % (i+1, keys[i])
    
        for j in xrange(len(keys)):
            color = bgcolor(mat[keys[i]][keys[j]])
            print >>out, "<td bgcolor='%s'><a title='%s,%s'>%d</a></td>" % \
                (color, keys[i], keys[j], int(mat[keys[i]][keys[j]]))
        print >>out, "</tr>"
    print >>out, "</table>"
    
    
    print >>out, """
    <body>
    </html>
    """
    
    out.close()



