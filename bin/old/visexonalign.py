#!/usr/bin/env summon

# python libs
import sys

# rasmus libs
from rasmus import fasta, util, genomeio, genomeutil

# graphics libs
from summon import *
from summonlib import *
from summonlib import shapes
from rasmus.vis import visual



options = [
 ["a:", "align=", "align", "AUTO<fasta alignment>"],
]

try:
    param, rest = util.parseArgs(sys.argv[1:], options)
except:
    sys.exit(1)




class ExonVis:
    def __init__(self, aln, genomes):
        self.aln = aln
        self.genomes = genomes
        self.conf = self.initConf()
        self.dnaGroup = None
        self.showDNA = True
        self.rulerGroup = None
        self.showRuler = True
        self.baseGroup = None
        
        self.install_bindings()
        

    def initConf(self, conf = None):
        if conf == None:
            conf = {}
            
        conf["color-seq"] = color(.9,.9,1)
        conf["color-unseq"] = color(.7,.7,.7)
        conf["color-exon"] = color(1,0,0,.7)
        conf["color-text"] = color(0,0,0)
        
        return conf
    
    def install_bindings(self):
        init_binding(input_key("d"), self.toggleDNA)
        init_binding(input_key("r"), self.toggleRuler)
    
    def toggleDNA(self):
        self.showDNA = not self.showDNA
    
    def toggleRuler(self):
        self.showRuler = not self.showRuler
    
    def draw(self, cur, conf = None):
        if conf != None:
            self.conf = conf   
        
        self.baseGroup = cur
        
        # draw sequence
        for i in xrange(len(self.aln)):
            row = self.aln[i]
            
            cur2 = push_transform(cur, translate, 0, -i)
            insert_group(cur2, group(
                self.conf["color-text"],
                text(row.genome, -1e10, 0, 0, -1, "right", "middle")))
            self.drawSequence(cur2, row)

        
        for i in xrange(len(self.aln)):
            row = self.aln[i]
            cur2 = push_transform(cur, translate, 0, -i)
            self.drawExons(cur2, row, self.genomes[row.genome])

        insert_group(cur, group(color(0,0,0), 
            shapes.boxStroke(0, 0, len(aln[0].seq), -len(aln))))


    def drawSequence(self, cur, row):   
        vis = []
        
        # draw sequence 
        test = [(c != "-" and c != "N") for c in row.seq]
        ranges = util.islands(test)[True]        
        vis.append(self.conf["color-seq"])
        vis.extend(self.makeRanges(ranges))
        
        # draw unsequenced
        test = [(c == "N") for c in row.seq]
        
        try:
            ranges = util.islands(test)[True]
            vis.append(self.conf["color-unseq"])
            vis.extend(self.makeRanges(ranges))
        except:
            pass
        
        return insert_group(cur, group(* vis))

            
    def drawExons(self, cur, row, genome):
        # find visible exons
        exons = findExons(genome.chroms[row.chrom].genes, row.start, row.end)
        lookup = genomeutil.mkAlignLookup(row.seq)

        vis = [self.conf["color-exon"]]
        
        def mkfunc(exon):
            def func():
                self.printExon(exon)
            return func
        
        for exon in exons.values():
            start = genomeutil.global2align(exon.start, 
                                            row.start, row.end, row.strand,
                                            lookup)
            end = genomeutil.global2align(exon.end, 
                                          row.start, row.end, row.strand,
                                          lookup)
            
            if end < start:
                tmp = start; start = end; end = tmp
            
            
            
            #vis.append(shapes.boxStroke(start, -.01, end+1, -.99))
            vis.append(shapes.box(start, -.01, end+1, -.99))
            vis.append(hotspot("click", start, 0, end, -1, mkfunc(exon)))
        
        insert_group(cur, group( * vis))


    def makeRanges(self, ranges):
        vis = []
        for r in ranges:
            vis.append(quads(r[0], 0, 
                             r[1], 0,
                             r[1], -1,
                             r[0], -1))
        return vis


    def printExon(self, exon):
        print exon.name, exon.trans.name, exon.trans.gene.name, \
              exon.start, exon.end


    def drawDNA(self):
        if not self.showDNA:
            if self.dnaGroup != None:
                remove_group(self.dnaGroup)
                self.dnaGroup = None
            return    
        
        
        x1,y1,x2,y2 = map(int, get_visible())
        width, height = get_window_size()
        
        # clip DNA if too small
        if (x2 - x1) > (width / 6.0):
            if self.dnaGroup:
                remove_group(self.dnaGroup)
                self.dnaGroup = None
            return
        
        vis = [color(0,0,0)]
        
        for i in xrange(len(self.aln)):
            row = self.aln[i]
            for j in xrange(max(0,x1), min(len(row.seq), x2)):
                vis.append(text_scale(row.seq[j], j, -i, j+1, -i-1, "center", "middle"))
        
        if self.dnaGroup:
            self.dnaGroup = replace_group(self.dnaGroup, group(* vis))
        else:
            self.dnaGroup = insert_group(self.baseGroup, group(* vis))
            
    
    def drawRuler(self):
        if not self.showRuler:
            if self.rulerGroup != None:
                remove_group(self.rulerGroup)
                self.rulerGroup = None
            return
        if self.rulerGroup == None:
            self.rulerGroup = insert_group(self.baseGroup, group())
        
        
        
        self.rulerGroup = replace_group(self.rulerGroup,
            group(translate(-1, 0,
                visual.drawRuler(1, len(self.aln[0].seq)+1))))
        
        
    def update(self):
        self.drawDNA()
        self.drawRuler()
    

def findExons(genes, start, end):
    exons = {}
    for gene in genes:
        for trans in gene.trans.values():
            for exon in trans.exons:
                if (start <= exon.start <= end) or \
                   (start <= exon.end <= end):
                    exons[exon.name] = exon
    return exons




# read alignment
util.tic("read input")
alignFile = param["align"][-1]
aln = fasta.readAlignment(alignFile)

# read exon information
genomes = {}
for row in aln:
    util.tic("read %s exons" % row.genome)
    genomes[row.genome] = genomeutil.Genome(row.genome)
 
    genomeio.readEnsmartExons(genomes[row.genome], 
                              genomeio.structfile(row.genome, row.chrom))
    genomes[row.genome].autoconf()
    util.toc()
util.toc()


# print stats
print "alignment %s:" % alignFile
print "%d sequences of length %d" % (len(aln), len(aln[0].seq))


# draw
set_bgcolor(1,1,1)
clear_groups()
set_antialias(False)
vis = ExonVis(aln, genomes)
vis.draw(push_transform(get_root_id(), scale, 1, 100))
home()


# continuous updating
add_update_func(vis.update)
begin_updating()
#def update():
#    vis.update()
#    timer_call(.5, update)
#update()

