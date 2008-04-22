# python libs
import math
import itertools

# summon libs
from summon.core import *
#from summon.simple import *
import summon
from summon import shapes, colors

# rasmus libs
from rasmus import algorithms
from rasmus import util
from rasmus import regionlib
from rasmus.bio import clustalw
from rasmus.bio import fasta



class SyntenyBlock:
    def __init__(self, region1, region2):
        self.region1 = region1
        self.region2 = region2
    
    def get_direction(self):
        """Returns 1 if region1.strand == region2.strand and
                  -1 if region1.strand != region2.strand
        """
        return self.region1 * self.region2
    

def readSyntenyBlocks(filename, feature="synteny"):
    infile = util.openStream(filename)
    blocks = []
    
    for line in infile:
        species1, chrom1, start1, end1, \
        species2, chrom2, start2, end2, direction = line.split("\t")
        
        blocks.append(SyntenyBlock(
            regionlib.Region(species1, chrom1, feature, 
                             int(start1), int(end1), 1),
            regionlib.Region(species2, chrom2, feature, 
                             int(start2), int(end2), int(direction))))
    return blocks




class Plot (object):
    def __init__(self, regions1, regions2, hits, hitnames=True, 
                 style="line", color=color(0, 0, 0),
                 selfhits=True):
        self.regions1 = regions1
        self.regions2 = regions2
        self.style = style
        self.color = color
        self.selfhits = selfhits
        
        if hitnames:
            self.hits = []
        
            # resolve hits to regions
            name2region = util.Dict(default=[])
            for region in itertools.chain(self.regions1, self.regions2):
                name2region[region.data["ID"]].append(region)
            
            for hit in hits:
                newhit = []
                for name in hit:
                    if name in name2region:
                        newhit.extend(name2region[name])
                if len(newhit) > 0:
                    self.hits.append(newhit)
                
        else:
            self.hits = hits


class SyntenyPlot (Plot):
    def __init__(self, blocks, **options):
        options.setdefault("style", "box")
        options.setdefault("color", color(0, 1, 0, .5))
        regions1 = []
        regions2 = []
        hits = []
        
        for block in blocks:
            regions1.append(block.region1)
            regions2.append(block.region2)
            regions1.append(block.region2)
            regions2.append(block.region1)
            hits.append([block.region1, block.region2])
        
        Plot.__init__(self, regions1, regions2, hits, hitnames=False, **options)



class DotplotMenu (summon.SummonMenu):
    """summatrix popup menu"""
    
    def __init__(self, viewer):
        summon.SummonMenu.__init__(self, viewer.win)
        
        self.viewer = viewer
        
        self.dotplot_menu = summon.Menu()
        self.dotplot_menu.add_entry("toggle chroms (c)", viewer.toggleChromDivs)
        
        
        self.plots = []
        def func(plot):
            return lambda: viewer.toggleTrace(plot)
        
        for i, plot in enumerate(viewer.plots):
            self.plots.append(summon.Menu())
            #self.plots[-1].add_entry("toggle", viewer.toggleRegions)
            self.plots[-1].add_entry("trace", func(plot))
            self.dotplot_menu.add_submenu("plot %i" %(i+1), self.plots[-1])
        self.insert_submenu(0, "Dotplot", self.dotplot_menu)
        
  

class Dotplot (object):
    def __init__(self, chroms1, chroms2):
        self.chroms1 = chroms1
        self.chroms2 = chroms2
        
        self.plots = []
        
        # visualization objects
        self.showChromDivs = 0
        self.chromDivs = None
        self.genomeDivs = None
        self.traces = {}
        
        # selection
        self.selgenes1 = []
        self.selgenes2 = []
    
        # colors
        self.colorChromDiv = color(1,0,0,.5)
        self.colorGenomeDiv = color(0, 0, 0)
        self.colorRegion = color(0, 1, 1, .2)
        self.colorRegionHit = color(0, 0, 1, .2)
    
        # layouts 
        self.chrom1Layout = {}
        self.chrom2Layout = {}        
        self.layout1 = {}
        self.layout2 = {}        
        self.plotSize = [0, 0]

        # create chrom lookup
        self.chrom1Lookup = {}
        self.chrom2Lookup = {}        
        for chrom in self.chroms1:
            self.chrom1Lookup[(chrom.species, chrom.seqname)] = chrom
        for chrom in self.chroms2:
            self.chrom2Lookup[(chrom.species, chrom.seqname)] = chrom
        


    def addPlot(self, plot):
        self.plots.append(plot)
                
    
    def show(self, winsize=(600,600), winpos=None):
        self.win = summon.Window("dotplot", size=winsize, position=winpos)
        self.win.set_bgcolor(1, 1, 1)
        self.win.add_group(self.draw())
        self.win.home()
        
        self.menu = DotplotMenu(self)
        self.win.set_menu(self.menu)
        
        self.win.set_binding(input_key("c"), self.toggleChromDivs)
        self.win.set_binding(input_key("t"), self.toggleTrace)        
    
    
    def draw(self):
        self.layout()
        
        # draw plots
        vis = group()
        for plot in self.plots:
            vis.append(self.drawPlot(plot))

        # draw genomes and chromosomes
        self.chromDivs = self.drawChromBorders()
        self.genomeDivs = self.drawGenomeBorders()
        vis.append(group(
            self.chromDivs,
            self.genomeDivs,
            hotspot("click",
                    0, 0,
                    self.plotSize[0], self.plotSize[1],
                    self.onClick,
                    give_pos=True)))
        return vis


    def layout(self):
        self.chrom1Layout = self.layoutChroms(self.chroms1)
        self.chrom2Layout = self.layoutChroms(self.chroms2)
        
        self.layout1.clear()
        self.layout2.clear()
        
        for plot in self.plots:
            self.layoutPlot(plot)
        
        self.plotSize = [self.chrom1Layout[self.chroms1[-1]] + self.chroms1[-1].length(),
                         self.chrom2Layout[self.chroms2[-1]] + self.chroms2[-1].length()]        
    
    
    def layoutPlot(self, plot):
        self.layoutRegions(self.layout1, self.chrom1Lookup, plot.regions1, self.chrom1Layout)
        self.layoutRegions(self.layout2, self.chrom2Lookup, plot.regions2, self.chrom2Layout)
        

    def layoutRegions(self, layout, chromLookup, regions, chromLayout):
        for region in regions:
            chrom = chromLookup.get((region.species, region.seqname), None)
            
            if chrom in chromLayout and \
               util.overlap(chrom.start, chrom.end, region.start, region.end):
                chrompos = chromLayout[chrom]
                layout[region] = chrompos + region.start - chrom.start
        
    
    def layoutChroms(self, chroms):
        # determine chrom layout
        chromLayout = {}
        x = 0
        
        for chrom in chroms:
            chromLayout[chrom] = x
            x += chrom.length()
                
        return chromLayout
    
    
    def drawPlot(self, plot):
        if plot.style in ("line", "box"):
            vis = []

            # draw hits
            for hit in plot.hits:
                set1 = []
                set2 = []

                # split genes into sets (possibily overlapping)
                for region in hit:
                    if region in self.layout1:
                        set1.append(region)
                    if region in self.layout2:
                        set2.append(region)
                
                
                # draw all pairs of hits
                for region1 in set1:
                    chrom1 = self.chrom1Lookup[(region1.species, 
                                                region1.seqname)]
                    
                    for region2 in set2:
                        if not plot.selfhits and \
                           region1.data["ID"] == region2.data["ID"]:
                            continue
                    
                        chrom2 = self.chrom2Lookup[(region2.species, 
                                                    region2.seqname)]
                        
                        s1 = max(self.chrom1Layout[chrom1], 
                                 self.layout1[region1])
                        e1 = min(self.chrom1Layout[chrom1] + chrom1.length(),
                                 self.layout1[region1] + region1.length())
                        s2 = max(self.chrom2Layout[chrom2], 
                                 self.layout2[region2])
                        e2 = min(self.chrom2Layout[chrom2] + chrom2.length(),
                                 self.layout2[region2] + region2.length())
                        
                        if plot.style == "line":
                            if region1.strand == region2.strand:
                                vis.extend([s1, s2, e1, e2])
                            else:
                                vis.extend([s1, e2, e1, s2])
                        elif plot.style == "box":
                            vis.append(shapes.box(s1, s2, e1, e2, fill=False))
                        else:
                            raise Exception("unknown plot style '%s'" % plot.style)
            
            if plot.style == "line":
                return group(plot.color, lines(*vis))
            elif plot.style == "box":
                return group(plot.color, *vis)
        else:
            return group()

    
    def drawChromBorders(self):
        vis = []
                
        # determine chrom layout
        divx = [0]
        divy = [0]
        
        for chrom in self.chroms1:
            divx.append(divx[-1] + chrom.length())
        for chrom in self.chroms2:
            divy.append(divy[-1] + chrom.length())
            
        maxx = divx[-1]
        maxy = divy[-1]
        
        # draw dividers
        for x in divx:
            vis.extend([x, 0, x, maxy])
        for y in divy:
            vis.extend([0, y, maxx, y])
        
        return group(self.colorChromDiv, lines(* vis))
    
    
    def drawGenomeBorders(self):
        vis = []
        
        assert len(self.chroms1) > 0 and len(self.chroms2) > 0
        
        # determine chrom layout
        divx = [0]
        divy = [0]
        x = self.chroms1[0].length()
        y = self.chroms2[0].length()
        
        for i in xrange(1, len(self.chroms1)):
            x += self.chroms1[i].length()
            if self.chroms1[i].species != self.chroms1[i-1].species:
                divx.append(x)
        for i in xrange(1, len(self.chroms2)):
            y += self.chroms2[i].length()
            if self.chroms2[i].species != self.chroms2[i-1].species:
                divy.append(y)
        
        divx.append(x)
        divy.append(y)
            
        maxx = x
        maxy = y
        
        # draw dividers
        for x in divx:
            vis.extend([x, 0, x, maxy])
        for y in divy:
            vis.extend([0, y, maxx, y])
        
        return group(self.colorGenomeDiv, lines(* vis))


    def drawTrace(self, plot):
        vis = []
        
        # get current view
        view = self.win.get_visible()
        
        for region in plot.regions1:
            if region in self.layout1:
                start = self.layout1[region]
            else:
                continue

            # only draw regions in view
            end = start + region.length()
            if util.overlap(view[0], view[2], start, end):
                vis.extend([self.colorRegion, 
                            shapes.box(start, 0, end, self.plotSize[1])])

        for region in plot.regions2:
            if region in self.layout2:
                start = self.layout2[region]
            else:
                continue

            # only draw regions in view
            end = start + region.length()
            if util.overlap(view[1], view[3], start, end):
                vis.extend([self.colorRegion, 
                            shapes.box(0, start, self.plotSize[0], end)])
        
        return group(*vis)
    
    


    def getGenesByPos(self, geneLayout, pos):
        gene_lst = []
        
        for gene, gene_pos in geneLayout.iteritems():
            if gene_pos <= pos <= gene_pos + gene.length():
                gene_lst.append(gene)
        
        return gene_lst

    def onClick(self, clickx, clicky):
        # determine genes that have been clicked
        self.selgenes1 = self.getGenesByPos(self.layout1, clickx)
        self.selgenes2 = self.getGenesByPos(self.layout2, clicky)
        
        print
        for gene in self.selgenes1:
            if "ID" not in gene.data:
                continue
            
            print "X %10s (%s %s:%s-%s)" % (gene.data["ID"],
                gene.species, 
                gene.seqname,
                util.int2pretty(gene.start),
                util.int2pretty(gene.end))
        for gene in self.selgenes2:
            if "ID" not in gene.data:
                continue
        
            print "Y %10s (%s %s:%s-%s)" % (gene.data["ID"],
                gene.species, 
                gene.seqname,
                util.int2pretty(gene.start),
                util.int2pretty(gene.end))


    def toggleChromDivs(self):
        self.showChromDivs = (self.showChromDivs + 1) % 3
        
        if self.showChromDivs == 0:
            self.chromDivs.set_visible(True)
            self.genomeDivs.set_visible(True)            
        if self.showChromDivs == 1:
            self.chromDivs.set_visible(False)
            self.genomeDivs.set_visible(True)            
        if self.showChromDivs == 2:
            self.chromDivs.set_visible(False)
            self.genomeDivs.set_visible(False)            


    def enableTrace(self, show, plot=None):
        if len(self.plots) == 0:
            return
        if plot == None:
            plot = self.plots[0]
    
        if show:
            vis = self.drawTrace(plot)

            # remove old regions
            if plot not in self.traces:
                self.traces[plot] = self.win.add_group(vis)
            else:
                self.traces[plot] = self.win.replace_group(self.traces[plot], vis)
        
        elif plot in self.traces:
            self.win.remove_group(self.traces[plot])
            del self.traces[plot]


    def toggleTrace(self, plot=None):
        if len(self.plots) == 0:
            return
        if plot == None:
            plot = self.plots[0]
    
        self.enableTrace(plot not in self.traces, plot)




