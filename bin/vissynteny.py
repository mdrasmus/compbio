#!/usr/bin/env summon

# python libs
import os, sys

# third party libs
from reportlab import svglib

# rasmus libs
from rasmus import util, env

# rasmus bio libs
from rasmus.bio import genomeutil, genomeio, muscle, ensembl, fasta, alignlib
from rasmus.bio.genomeutil import *

# for easy user iteraction
from rasmus.common import *


# graphics libs
from summon.core import *
from rasmus.vis import syntenyvis as genomevis
from summon.colors import *

from summon import svg


options = [
    "Required options",
    ["g:", "genomes=", "genomes", "<genome>,<genome>,..."],
    ["s:", "synteny=", "synteny", "<synteny>"], 
    ["c:", "orthcomp=", "orthcomp", "<orth comps>"],
    ["S:", "smap=", "smap", "<gene2species map>",
        {"help": "mapping of gene names to species names"}],
    
    "Sequence loading",
    ["q", "sequence", "sequence", "", 
        {"help": "automatically read sequences"}],
    ["f:", "fasta=", "fasta", "<fasta file>",
        {"help": "loads a fasta file into sequence dictionary"}],
    
    "Automatic drawing",
    ["v:", "visual=", "visual", "<outputfile>"],
    ["w:", "winsize=", "winsize", "<window width>x<window height>"],
    ["p:", "part=", "part", "<part file>", 
        {"help": "draw every synteny group in a part file"}],
    ["W:", "chromWindows=", "chromWindows", "<window index filename>",
        {"help": "draw every window along a chromosome"}],
    ["r:", "refgenome=", "refgenome", "<reference genome>",
        {"help": "which genome should be drawn on top (reference)"}],
    ["e:", "gene=", "gene", "<gene to draw>"],
    ["x:", "context=", "context", "<context size (bp)>"],    
    
    "Gene label options",
    ["l", "scaleLabels", "scaleLabels", ""],
    ["L", "fixLabels", "fixLabels", ""],
    
    "Misc. options",
    ["P:", "paths=", "paths", "<data paths>", 
        {"default": ".",
         "single": True}]
]



# check args
param = util.parseOptions(sys.argv, options, quit=True, resthelp="<scripts> ...")


# global variables
markColor = color(0,0,1)
context = 1e6
selgene = None
selgenes = []

genes = {}
lookup = {}
comps = []
m = None
gene2species = None
vis = None
conf = {}
seqs = FastaDict()
winsize = (800, 400)



##################################################################
# Classes and functions
#

class SyntenyVis (genomevis.SyntenyVis):
    def __init__(self, conf, matching, **options):
        genomevis.SyntenyVis.__init__(self, conf, matching, **options)
        self.clickMode = "gene"
    
    def geneClick(self, gene):
        selgene = gene
        selgenes.append(gene)
        
        if self.clickMode == "gene":
            vis.printGene(gene)
        elif self.clickMode == "view":
            drawGene(gene.name, context)
        elif self.clickMode == "sequence":
            if gene.name in seqs:
                print
                print "%s\n%s" % (gene.name, seqs[gene.name])
                print
            else:
                print "%s has no sequence" % gene.name
                
        elif self.clickMode == "align":
            if gene.name not in lookup:
                print "gene %s has no matches" % gene.name
                return    
            
            comp = filter(lambda x: x in genes, comps[lookup[gene.name]])
            seqs2 = util.subdict(seqs, comp)
            aln = muscle.muscle(seqs2)
            
            keys = aln.keys()
            
            for key in keys:
                if "y" not in dir(genes[key]):
                    genes[key].y = -1e1000
            
            keys.sort(lambda a,b: cmp(genes[b].y, genes[a].y))
            alignlib.printAlign(aln, order=keys)


def draw(genome, chrom, start, end):
    util.tic("draw")    
    vis.draw(genome, chrom, start, end)
    util.toc()

def drawGene(geneName, context = context):
    gene = genes[geneName]
    draw(gene.chrom.genome.name, gene.chrom.name, 
         gene.start - context, gene.end + context)
    vis.mark(geneName, shape="box", col=markColor)

def drawAll(genome):
    vis.drawAll(genome)


def align(* names):
    if len(names) == 0:
        print "nothing to align"
    
    # get names from genes if they are not strings
    if type(names[0]) != str:
        names = [i.name for i in names]
    
    seqs2 = util.subdict(seqs, names)
    aln = muscle.muscle(seqs2)
    muscle.printAlign(aln)

def mark(shape="box", col=blue):
    names = []
    while True:
        line = sys.stdin.readline().rstrip()
        if line == "": break
        names.append(line)
    markGenes(names, shape, col)

def markGenes(names, shape="box", col=blue):
    for name in names:
        vis.mark(name, shape, col)

def markHoles(shape="box", col=blue):
    genes2 = filter(lambda x: len(x.matches) == 0, genes.values())
    names = [x.name for x in genes2]
    markGenes(names, "box", col)

def find(name):
    return vis.find(name)

printscreen = lambda *args, **kargs: svg.printScreen(vis.win, *args, **kargs)


def readFasta(filename):
    seqs.update(fasta.readFasta(env.findFile(f)))

def readAllSeqs():
    util.tic("read sequences")

    for genome in m.getGenomeOrder():
        try:
            seqfile = env.findFile("%s.fasta" % genome)
            util.tic("reading '%s'" % seqfile)
            seqs.read(seqfile)
            util.toc()
        except: 
            util.log("cannot read fasta '%s.fasta'" % genome)
    util.toc()


#
# Visualize ortholog sets across the genome
#
def visparts(parts, refGenome, outdir):
    for i in xrange(len(parts)):
        util.log("visualizing part %d" % i)
        # find the biggest group of gene from the reference genome from the same 
        # chromosome
        part = parts[i]
        refs = filter(lambda gene: genes[gene].chrom.genome.name == refGenome,
                      part)
        if len(refs) == 0:
            continue
        counts = util.histDict([genes[ref].chrom for ref in refs])
        keys = counts.keys()
        keys.sort(lambda a,b: cmp(counts[b], counts[a]))
        refChrom = keys[0]
        refgenes = filter(lambda gene: genes[gene].chrom == refChrom, part)
        
        start = min([genes[gene].start for gene in refgenes]) - context
        end = max([genes[gene].end for gene in refgenes]) + context
        
        draw(refGenome, refChrom.name, start, end)
        markGenes(part, markColor, "box")
        vis.win.set_visible(start, 2*conf["gene-size"], end, 
                            -conf["max-genome-sep"] * len(m.genomes))
                    
        # output svg 
        svg.printScreen(param["visual"][-1] + "/synteny%d.svg" % i)
        
        # conversion
        #os.system("convert %s %s" % (
        #    param["visual"][-1] + "/synteny%d.svg" % i,
        #    param["visual"][-1] + "/synteny%d.png" % i))
        svglib.svg2pdf(param["visual"][-1] + "/synteny%d.svg" % i)
        
        os.remove(param["visual"][-1] + "/synteny%d.svg" % i)

#
# Visualize windows across the genome
#
def viswindows(refGenome, windowSize, windowStep, outdir):
    chroms = m.genomes[refGenome].chroms.values()
    chroms.sort(lambda a,b: cmp(a.name, b.name))
    
    
    indexfile = file(param["chromWindows"][-1], "w")
    print >>indexfile, "##version:1.0"
    print >>indexfile, "##types:string\tstring\tstring\tint\tint\tstring"
    print >>indexfile, "file\tgenome\tchrom\tstart\tend\tgenes"
    
    
    for chrom in chroms:
        i = 0
        for start in xrange(0, int(chrom.size), windowStep):
            end = start + windowSize
            draw(refGenome, chrom.name, start, end)
            #visgroup = vis.draw(conf, m, refGenome, chrom.name, start, end)

            vis.win.set_visible(start, 10*conf["gene-size"], end, 
                                -conf["max-genome-sep"] * len(m.genomes))
            
            filename = ("%s_%s_%s-%s.svg" % 
                (refGenome, chrom.name, 
                 util.int2pretty(int(start)),
                 util.int2pretty(int(end))))
            
            # figure out which genes are in view
            visgenes = filter(lambda gene: 
                start < gene.x < end, vis.placedGenes.keys())
            
            # record in lookup filename->genes
            #pdffile = os.path.split(filename)[1].replace(".svg", ".pdf")
            #lookup.append([pdffile] + [x.name for x in visgenes])
            indexfile.write("\t".join([filename,
                                       refGenome,
                                       chrom.name,
                                       str(int(start)),
                                       str(int(end)),
                                       ",".join([x.name for x in visgenes])])
                            +"\n")
            
            # output svg 
            svg.printScreen(vis.win, outdir + "/" + filename) #, visgroup)
            
            # conversion
            #grid.execute("svg2pdf.py %s -r" % filename)
            
            i += 1
    indexfile.close()
    

# add key bindings
def press(mode):
    def func():
        print "mode is '%s'" % mode
        vis.clickMode = mode
    return func

def clear_selgenes():
    selgenes[:] = []
    print "selgenes cleared"

def align_selgenes():
    align(* selgenes)
    
def print_selgenes():
    print selgenes


def readData(genomes, compfile, syntenyfile, smapfile):
    global genes
    global lookup
    global comps
    global m
    global gene2species
    global conf
    global vis
    global winsize
    
    util.tic("read")
    
    # read genomes
    gene2species = genomeutil.readGene2species(env.findFile(smapfile))
    m = Matching()
    genomeio.readGenomes(m, genomes, gene2species)
    genes = m.getGenes()

    # read orthologs
    util.tic("read ortholog components")
    lookup = {}
    comps = util.readDelim(env.findFile(compfile))
    comps = map(lambda comp: filter(lambda gene: gene in genes, comp), comps)
    comps = filter(lambda comp: len(comp) > 0, comps)
    m.setGeneComponents(comps)
    util.toc()
    
    # read synteny blocks
    util.tic("read synteny blocks")
    genomeio.readBlockDimensions(m, env.findFile(syntenyfile), True)
    util.toc()
    
    util.toc()
    
    # setup visualization
    conf.update(genomevis.initConf({}, genomevis.calcGeneHeight(m)))
    vis = SyntenyVis(conf, m, winsize=winsize)



vis.win.set_binding(input_key("g"), press("gene"))
vis.win.set_binding(input_key("v"), press("view"))
vis.win.set_binding(input_key("s"), press("sequence"))
vis.win.set_binding(input_key("a"), press("align"))
vis.win.set_binding(input_key("d"), clear_selgenes)
vis.win.set_binding(input_key("w"), align_selgenes)
vis.win.set_binding(input_key("e"), print_selgenes)



#--------------------------------------------------------------------------------
# main function
#--------------------------------------------------------------------------------

def main(param):
    global conf
    global context
    global seqs
    global winsize    

    # init
    print genomevis.BINDINGS_HELP
    print """
Modes
g     display gene info mode
v     gene viewing mode (display gene with context)
s     display gene sequence mode
a     display gene alignment mode

Gene array
d     clear selgene array (selected genes)
e     print selgene array
w     print alignment of selgene array

    """

    env.addPaths(param["paths"])
    env.addEnvPaths("DATAPATH")

    # stop if no options are given
    if "smap" in param and \
       "genomes" in param and \
       "orthcomp" in param and \
       "synteny" in param:
        

        if "winsize" in param:
            winsize = map(int, param["winsize"][-1].split("x"))
            #set_window_size(* winsize)
        
        
        # read data
        readData(param["genomes"][-1].split(","), 
                 param["orthcomp"][-1],
                 param["synteny"][-1],
                 param["smap"][-1])

        seqs = fasta.FastaDict()
        globals()["seqs"] = seqs
        if "sequence" in param:
            readAllSeqs()
        
        
        # configuration    
        if "scaleLabels" in param:
            conf["show-gene-labels"] = "scale"
        
        if "fixLabels" in param:
            conf["show-gene-labels"] = "fix"

        if "context" in param:
            context = float(param["context"][-1])



        # automatic execution
        if "gene" in param:
            drawGene(param["gene"][-1], context)
            gene = genes[param["gene"][-1]]
            home()
            coords = get_visible()
            top = 2 * vis.conf["gene-size"] 
            bottom = -len(m.genomes) * vis.conf["max-genome-sep"]
            vis.win.set_visible(gene.start - context, top, gene.end + context, bottom)


        if "visual" in param:
            from summonlib import svg    
            if "part" in param and "refgenome" in param:
                parts = util.readDelim(param["part"][-1])
                refGenome = param["refgenome"][-1]
                visparts(parts, refGenome, param["visual"][-1])
                sys.exit(0)
            elif "chromWindows" in param and "refgenome" in param:
                refGenome = param["refgenome"][-1]
                viswindows(refGenome, context, param["visual"][-1])
                sys.exit(0)
            else:
                svg.printScreen(vis.win, param["visual"][-1])
                sys.exit(0)


    else:
        print "vissynteny: Not all required options are given.  Use commandline to load data"
    

    #
    # executing scripts
    #
    for f in param[""]:
        util.tic("executing '%s'" % f)
        execfile(f)
        util.toc()


main(param)



