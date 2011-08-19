"""

    Dotplot visualization


"""



# python libs
import math
import itertools

# summon libs
from summon.core import *
#from summon.simple import *
import summon
from summon import shapes, colors

# rasmus libs
from rasmus import util
from compbio import regionlib
from compbio.synteny import SyntenyBlock
#from compbio.synteny import read_synteny_blocks as readSyntenyBlocks


def make_chroms(genes):
    """Make default chromosomes for a set of genes"""
    
    chroms = {}
    for gene in genes:
        if gene.seqname in chroms:
            chrom = chroms[gene.seqname]
            chrom.end = max(chrom.end, gene.end)
        else:
            chrom = regionlib.Region(gene.species, gene.seqname, "chromosome",
                                     1, gene.end)
            chroms[gene.seqname] = chrom

    return sorted(chroms.values(), key=lambda x: x.end, reverse=True)
   



class Plot (object):
    def __init__(self, regions1, regions2, hits, hitnames=True, 
                 style="line", color=(0, 0, 0),
                 fill_color=None,
                 trace_color=(0, 1, 1, .2),
                 selfhits=True, name=None):
        self.name = name
        self.regions1 = regions1
        self.regions2 = regions2
        self.style = style
        self.color = color
        self.fill_color = fill_color
        self.trace_color = trace_color
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

    
    def draw_plot(self, dotplot):
        
        if self.style not in ("line", "box") and not callable(self.style):
            return group()

        vis = group()
        line_pts = []

        # draw hits
        for hit in self.hits:
            set1 = []
            set2 = []

            # split genes into sets (possibily overlapping)
            for region in hit:
                if region in dotplot.layout1:
                    set1.append(region)
                if region in dotplot.layout2:
                    set2.append(region)


            # draw all pairs of hits
            for region1 in set1:
                chrom1 = dotplot.chrom1lookup[(region1.species, 
                                            region1.seqname)]

                for region2 in set2:
                    if not self.selfhits and \
                       region1.data["ID"] == region2.data["ID"]:
                        continue

                    chrom2 = dotplot.chrom2lookup[(region2.species, 
                                                   region2.seqname)]

                    s1 = dotplot.layout1[region1]
                    e1 = dotplot.layout1[region1] + region1.length()
                    s2 = dotplot.layout2[region2]
                    e2 = dotplot.layout2[region2] + region2.length()

                    # styles
                    if self.style == "line":
                        if region1.strand == region2.strand:
                            line_pts.extend([s1, s2, e1, e2])
                        else:
                            line_pts.extend([s1, e2, e1, s2])

                    elif self.style == "box":
                        vis.append(shapes.box(s1, s2, e1, e2, fill=False))

                        if self.fill_color:
                            vis.append(color(*self.fill_color))
                            vis.append(shapes.box(s1, s2, e1, e2,
                                                  fill=True))
                            vis.append(color(*self.color))

                    elif callable(self.style):
                        vis.append(self.style(
                                region1, s1, s2, region2, e1, e2))

                    else:
                        raise Exception("unknown plot style '%s'" % self.style)

        if self.style == "line":
            return group(color(*self.color), lines(*line_pts))
        elif self.style == "box":
            return group(color(*self.color), vis)
        elif callable(self.style):
            return vis


    def draw_trace(self, dotplot):
        vis = []
        
        # get current view
        view = dotplot.win.get_visible()
        
        for region in self.regions1:
            if region in dotplot.layout1:
                start = dotplot.layout1[region]
            else:
                continue

            # only draw regions in view
            end = start + region.length()
            if util.overlap(view[0], view[2], start, end):
                vis.extend([color(* self.trace_color),
                            lines(start, 0, start, dotplot.plot_size[1]),
                            shapes.box(start, 0, end, dotplot.plot_size[1])])

        for region in self.regions2:
            if region in dotplot.layout2:
                start = dotplot.layout2[region]
            else:
                continue

            # only draw regions in view
            end = start + region.length()
            if util.overlap(view[1], view[3], start, end):
                vis.extend([color(* self.trace_color),
                            lines(0, start, dotplot.plot_size[0], start),
                            shapes.box(0, start, dotplot.plot_size[0], end)])
        
        return group(*vis)
    



class SyntenyPlot (Plot):
    def __init__(self, blocks, **options):
        options.setdefault("style", "box")
        options.setdefault("color", (0, 1, 0, .5))
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
        self.dotplot_menu.add_entry("toggle chroms (c)",
                                    viewer.toggle_chrom_divs)
        
        
        self.plots = []
        def func(plot):
            return lambda: viewer.toggle_trace(plot)
        
        for i, plot in enumerate(viewer.plots):
            self.plots.append(summon.Menu())
            self.plots[-1].add_entry("trace", func(plot))

            if plot.name is None:
                plot_name = "plot %i" %(i+1)
            else:
                plot_name = plot.name
            
            self.dotplot_menu.add_submenu(plot_name, self.plots[-1])
        self.insert_submenu(0, "Dotplot", self.dotplot_menu)
        
  

class Dotplot (object):
    
    def __init__(self, chroms1, chroms2, labels=True, get_label=None,
                 chrom_div_color=(1, 0, 0, .5),
                 genome_div_color=(0, 0, 0)):
        self.chroms1 = chroms1
        self.chroms2 = chroms2
        self.show_labels = labels

        if get_label is None:
            self.get_label = lambda x: x.seqname
        else:
            self.get_label = get_label
        
        self.plots = []
        
        # visualization objects
        self.show_chrom_divs = 0
        self.chrom_divs = None
        self.genome_divs = None
        self.traces = {}
        
        # selection
        self.selgenes1 = []
        self.selgenes2 = []
    
        # colors
        self.chrom_div_color = chrom_div_color
        self.genome_div_color = genome_div_color
    
        # layouts 
        self.chrom1layout = {}
        self.chrom2layout = {}        
        self.layout1 = {}
        self.layout2 = {}        
        self.plot_size = [0, 0]

        # create chrom lookup
        self.chrom1lookup = {}
        self.chrom2lookup = {}
        for chrom in self.chroms1:
            self.chrom1lookup[(chrom.species, chrom.seqname)] = chrom
        for chrom in self.chroms2:
            self.chrom2lookup[(chrom.species, chrom.seqname)] = chrom
        


    def add_plot(self, plot):
        self.plots.append(plot)
        return plot
    addPlot = add_plot # back-compatiable
                
    
    def show(self, winsize=(600,600), winpos=None):
        self.win = summon.Window("dotplot", size=winsize, position=winpos)
        self.win.set_bgcolor(1, 1, 1)
        self.win.add_group(self.draw())
        self.win.home()
        
        self.menu = DotplotMenu(self)
        self.win.set_menu(self.menu)
        
        self.win.set_binding(input_key("c"), self.toggle_chrom_divs)
        self.win.set_binding(input_key("t"), self.toggle_trace)        
    
    
    def draw(self):
        self.layout()
        
        # draw plots
        vis = group()
        for plot in self.plots:
            vis.append(plot.draw_plot(self))

        # draw genomes and chromosomes
        self.chrom_divs = self.draw_chrom_borders()
        self.genome_divs = self.draw_genome_borders()
        vis.append(group(
            self.chrom_divs,
            self.genome_divs,
            hotspot("click",
                    0, 0,
                    self.plot_size[0], self.plot_size[1],
                    self.on_click,
                    give_pos=True)))
        
        if self.show_labels:
            vis.append(self.draw_labels())
        
        return vis


    def layout(self):
        self.chrom1layout = self.layout_chroms(self.chroms1)
        self.chrom2layout = self.layout_chroms(self.chroms2)
        
        self.layout1.clear()
        self.layout2.clear()
        
        for plot in self.plots:
            self.layout_plot(plot)
        
        self.plot_size = [self.chrom1layout[self.chroms1[-1]] +
                          self.chroms1[-1].length(),
                          self.chrom2layout[self.chroms2[-1]] +
                          self.chroms2[-1].length()]
    
    
    def layout_plot(self, plot):
        self.layout_regions(self.layout1, self.chrom1lookup,
                            plot.regions1, self.chrom1layout)
        self.layout_regions(self.layout2, self.chrom2lookup,
                            plot.regions2, self.chrom2layout)
        

    def layout_regions(self, layout, chromLookup, regions, chromLayout):
        for region in regions:
            chrom = chromLookup.get((region.species, region.seqname), None)
            
            if chrom in chromLayout and \
               util.overlap(chrom.start, chrom.end, region.start, region.end):
                chrompos = chromLayout[chrom]
                layout[region] = chrompos + region.start - chrom.start
        
    
    def layout_chroms(self, chroms):
        # determine chrom layout
        chromLayout = {}
        x = 0
        
        for chrom in chroms:
            chromLayout[chrom] = x
            x += chrom.length()
                
        return chromLayout


    
    def draw_chrom_borders(self):
        vis = [color(*self.chrom_div_color)]
                
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
        
        return lines(* vis)
        



    def draw_labels(self):
        vis = group(color(*self.genome_div_color))
        thick = 1000000
    
        # determine chrom layout
        divx = [0]
        divy = [0]
        labelsx = []
        labelsy = []
        
        for chrom in self.chroms1:
            labelsx.append(self.get_label(chrom))
            divx.append(divx[-1] + chrom.length())
        for chrom in self.chroms2:
            labelsy.append(self.get_label(chrom))
            divy.append(divy[-1] + chrom.length())
            
        maxx = divx[-1]
        maxy = divy[-1]
        
        last = 0
        for x, label in zip(divx[1:], labelsx):
            vis.append(text_clip(label, last, 0, x, -thick, 4, 20, "left", "top"))
            last = x

        last = 0
        for y, label in zip(divy[1:], labelsy):            
            vis.append(translate(0, last,
                rotate(90,
                    text_clip(label, 0, 0, thick, y-last, 4, 20, "left", "bottom"))))
            last = y
        return vis
            
    
    def draw_genome_borders(self):
        vis = [color(*self.genome_div_color)]
        
        assert len(self.chroms1) > 0 and len(self.chroms2) > 0
        
        # determine chrom layout
        divx = [0]
        divy = [0]
        x = self.chroms1[0].length()
        y = self.chroms2[0].length()
        
        for i in xrange(1, len(self.chroms1)):
            if self.chroms1[i].species != self.chroms1[i-1].species:
                divx.append(x)                
            x += self.chroms1[i].length()
        for i in xrange(1, len(self.chroms2)):
            if self.chroms2[i].species != self.chroms2[i-1].species:
                divy.append(y)
            y += self.chroms2[i].length()
        
        divx.append(x)
        divy.append(y)
        
        maxx = x
        maxy = y
        
        # draw dividers
        for x in divx:
            vis.extend([x, 0, x, maxy])
        for y in divy:
            vis.extend([0, y, maxx, y])

        return lines(* vis)
    


    def get_regions_by_pos(self, geneLayout, pos):
        gene_lst = []
        
        for gene, gene_pos in geneLayout.iteritems():
            if gene_pos <= pos <= gene_pos + gene.length():
                gene_lst.append(gene)
        
        return gene_lst

    def on_click(self, clickx, clicky):
        # determine regions that have been clicked
        self.selgenes1 = self.get_regions_by_pos(self.layout1, clickx)
        self.selgenes2 = self.get_regions_by_pos(self.layout2, clicky)
        
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


    def toggle_chrom_divs(self):
        self.show_chrom_divs = (self.show_chrom_divs + 1) % 3
        
        if self.show_chrom_divs == 0:
            self.chrom_divs.set_visible(True)
            self.genome_divs.set_visible(True)            
        if self.show_chrom_divs == 1:
            self.chrom_divs.set_visible(False)
            self.genome_divs.set_visible(True)            
        if self.show_chrom_divs == 2:
            self.chrom_divs.set_visible(False)
            self.genome_divs.set_visible(False)            


    def enable_trace(self, show, plot=None):
        if len(self.plots) == 0:
            return
        if plot is None:
            plot = self.plots[0]
    
        if show:
            vis = plot.draw_trace(self)

            # remove old regions
            if plot not in self.traces:
                self.traces[plot] = self.win.add_group(vis)
            else:
                self.traces[plot] = self.win.replace_group(self.traces[plot], vis)
        
        elif plot in self.traces:
            self.win.remove_group(self.traces[plot])
            del self.traces[plot]


    def toggle_trace(self, plot=None):
        if len(self.plots) == 0:
            return
        if plot == None:
            plot = self.plots[0]
        
        self.enable_trace(plot not in self.traces, plot)




