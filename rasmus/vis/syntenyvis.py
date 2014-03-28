
# python libs
import copy

# rasmus libs
from rasmus import stats
from rasmus import util

# rasmus bio libs
from compbio import fasta
from compbio import regionlib
import compbio.regionlib
from compbio.regionlib import iter_chrom


# graphics libs
from summon.core import *
from summon import shapes
import summon

# globals
gene_label_types = [False, "fix", "scale", "vertical"]


#
# Bindings
#

BINDINGS_HELP = """
==============================
  Synteny Visualization Keys
==============================
c     toggle controls
l     toggle gene labels
"""



def invcolor(c, alpha = 1):
    return color(1 - c[0], 1 - c[1], 1 - c[2], alpha)

   
light_colors   = [color(1, 0, 0, .6),
                  color(1, .2, 0, .6),
                  color(1, .5, 0, .6),
                  color(1, 1, 0, .6),
                  color(.5, 1, 0, .6),
                  color(.2, 1, 0, .6),
                  color(0, 1, 0, .6),
                  color(0, 1, .2, .6),
                  color(0, 1, .5, .6),
                  color(0, 1, 1, .6),
                  color(0, .5, 1, .6),
                  color(0, .2, 1, .6),
                  color(0, 0, 1, .6),
                  color(.2, 0, 1, .6),
                  color(.5, 0, 1, .6),
                  color(1, 0, 1, .6),
                  color(1, 0, .5, .6),
                  color(1, 0, .2, .6),
                  color(.5, .5, .5, .6),
                  color(.5, 0, 0, .6),
                  color(.5, .5, 0, .6),
                  color(0, .5, 0, .6),
                  color(0 , .5, .5, .6),
                  color(0, 0, .5, .6),
                  color(.5, 0, .5, .6)
                  ]

def draw_gene(vis, gene, direction, col):
    "draw a single gene"
    
    length = gene.length()
    height = 1
    steep = .1

    mid = max(length - (height * steep), 0)
    
    return group(polygon(col,  
                         0, 0,
                         mid, 0,
                         length, height/2,
                         mid, height,
                         0, height),
                     lines(0, 0, 0, height))

def draw_gene_box(vis, gene, direction, col):
    "draw a single gene"
    
    length = gene.length()
    height = 1
    
    return group(quads(col,  
                       0, 0,
                       length, 0,
                       length, height,
                       0, height),
                 lines(0, 0, 0, height))

        # determine which row gene is in (for color)
def gene_color_alt_species(vis, gene, eff_dir):
    order = 0
        
    if eff_dir == 1:
        if order % 2 == 0:
            col = vis.colors['gene_pos']
        else:
            col = vis.colors['gene2_pos']
            
    else:
        if order % 2 == 0:
            col = vis.colors['gene_neg']
        else:
            col = vis.colors['gene2_neg']
    return col

def effective_dir(direction, geneDirection):
    return (-1, 1)[(direction == 1) == (geneDirection == 1)]


#
# main visualization 
#
class Frag:
    "represents fragment of DNA"
    
    def __init__(self, genome=None, chrom=None, start=0, end=0, strand=0,
                 x=0, y=0):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.end = end
        self.direction = strand
        self.x = x
        self.y = y
        self.genes = []

    def length(self):
        return self.end - self.start + 1


class Layout (object):
    def __init__(self, x, y, orient):
        self.x = x
        self.y = y
        self.orient = orient



class SyntenyVisBase:    
    def __init__(self, genomes, chroms, regions, blocks, orths,
                 rootid=None,
                 winsize=(800, 400),

                 # graphics
                 gene_label=lambda x: x.data["ID"],
                 draw_gene=draw_gene_box,
                 gene_color=gene_color_alt_species,

                 # misc options
                 fat_matches  = True,
                 use_controls = False,
                 min_frag     = 0,
                 show_gene_labels = False,
                 font_size = 8,
 
                 # layout
                 genome_sep = 3,
                 max_genome_sep = 6,
                 frag_sep = 1.5,

                 # colors
                 color_genome_div = color(0, 0, 0, .5),
                 color_gene_pos   = color(1, .6, 0, .95),
                 color_gene_neg   = color(1, .6, 0, .95),
                 color_gene2_pos  = color(1, .2, 0, .95),
                 color_gene2_neg  = color(1, .2, 0, .95),
                 color_matches    = color(.8, .8, 1, .8),
                 color_arrow      = color(1, .8, 0, .5),
                 color_frag       = color(0, 0, 0, .8),
                 color_blocks     = [color(.8,.8,1,.5)]):

        self.win = None
        self.winsize = winsize   
        self.rootid = rootid
        
        self.ref_genome     = None
        self.frags         = set()
        self.controlids    = []
        self.markids       = []
        self.labelids      = []
        self.groupid       = 0
        self.visid = None
        
        self.gene_label = gene_label
        self.draw_gene = draw_gene
        self.gene_color = gene_color
        self.font_size = font_size
        
        self.region2frags = {}
        self.region_layout = {}

        
        self.genomes = genomes
        self.chroms = chroms
        self.chroms_lookup = dict(((x.species, x.seqname), x) for x in chroms)
        self.db = regionlib.RegionDb(regions)
        self.blocks = blocks


        # filter orths for only regions we have in our db
        self.orths = [[x for x in orth
                       if self.db.has_region(x)]
                      for orth in orths]
        
        # make ortholog lookup
        self.orth_lookup = {}
        for orth in self.orths:
            for region in orth:
                self.orth_lookup[region] = orth

        # options
        self.fat_matches = fat_matches
        self.use_controls = use_controls
        self.min_frag = min_frag
        self.show_gene_labels = show_gene_labels

        # layout
        self.genome_sep = genome_sep
        self.max_genome_sep = max_genome_sep
        self.frag_sep = frag_sep

        self.colors = {
            "genome_div" : color_genome_div,
            "gene_pos"   : color_gene_pos,
            "gene_neg"   : color_gene_neg,
            "gene2_pos"  : color_gene2_pos,
            "gene2_neg"  : color_gene2_neg,
            "matches"    : color_matches,
            "arrow"      : color_arrow,
            "frag"       : color_frag,
            "blocks"     : color_blocks
            }
                
        
    
    
    def show(self):
        if self.win == None or not self.win.is_open():
            self.win = summon.Window("synteny")
        
        if self.rootid == None:
            self.rootid = self.win.get_root()
        
        if self.visid == None:
            self.win.set_size(* self.winsize)
            self.win.set_bgcolor(1,1,1)
            self.win.set_binding(input_key("c"), self.toggle_controls)
            self.win.set_binding(input_key("l"), self.toggle_labels)
            self.visid = self.win.insert_group(self.rootid, group())

    def home(self):
        self.win.set_visible(*(self.win.get_root().get_bounding() + ("exact",)))
        self.win.zoom(.9, .9)

        
    def clear_drawing(self):
        self.frags = set()
        self.region_layout = {}
        self.win.remove_group(self.visid)
        self.clear_marks()
        self.visid = self.win.insert_group(self.rootid, group())
        
    
    
    def draw(self, genome, chrom, start, end, direction=1):
        self.show()
        self.clear_drawing()
        self.win.insert_group(self.visid, 
            self.draw_chrom(genome, chrom, 
                            start, end, direction=direction))
        self.show_controls(self.use_controls)
    
    '''
    def drawAll(self, refGenomeName):
        self.show()
        self.clear_drawing()
        
        y = 0
        chroms = self.matching.genomes[refGenomeName].chroms.values()
        chroms.sort(lambda a,b: cmp(b.size, a.size))
        
        for chrom in chroms:
            util.tic("drawing chrom %s" % chrom.name)
            self.win.insert_group(self.visid, group(translate(0, y, 
                self.draw_chrom(refGenomeName, chrom.name, 0, chrom.size))))
            self.show_controls(self.conf["use-controls"])
            util.toc()
            y -= self.conf['max-genome-sep'] * (len(self.matching.genomes) + 1)
    '''
    
    def draw_chrom(self, genome_name, chrom_name, start, end, direction=1):
        """Draw the synteny for a region of a chromosome"""
        
        self.ref_genome = genome_name

        self.layout_frags(genome_name, chrom_name, start, end, direction)
        return self.draw_placed()
        

    def layout_frags(self, genome_name, chrom_name, start, end, direction=1):

        ref_chrom  = self.chroms_lookup[(genome_name, chrom_name)]

        # setup genome display order
        order = {}
        for i, genome in enumerate(self.genomes):
            order[genome] = i
        
        # swap the genome with order 0 and the reference genome
        j = order[self.ref_genome]
        order[self.genomes[0]] = j
        order[self.ref_genome] = 0                
        
        # init reference fragment
        ref_frag = Frag(genome=genome_name,
                        chrom=chrom_name, 
                        start=max(start, 0),
                        end=min(end, ref_chrom.end),
                        strand=direction,
                        x=max(start,0),
                        y=0)
        self.frags.add(ref_frag)
        self.layout_frag_contents(ref_frag)
        
        
        # find all synteny blocks in this region
        # sort blocks by appearance in ref_chrom
        blocks = list(self.filter_blocks(self.blocks, ref_chrom, start, end))
        def blocksort(a):
            if a[1] == 0:
                starta = a[0].region1.start
            else:
                starta = a[0].region2.start
        blocks.sort(key=blocksort)
        
        
        # make lookup for genes to block and block to fragment
        block_lookup = {}
        frag_lookup = {}
        for block, flip in blocks:            
            if flip == 0:
                other = block.region2
            else:
                other = block.region1
                
            frag = Frag()
            frag.genome = other.species
            frag.chrom = other.seqname
            frag_lookup[block] = frag

            for gene2 in iter_chrom(self.db.get_regions(frag.genome, 
                                                        frag.chrom),
                                    other.start, other.end):
                block_lookup[gene2] = block
                
        self.block_lookup = block_lookup
        
        
        # find all genes that will be drawn
        # walk along ref_chrom and store drawn genes into fragments
        refLookup = {}
        for gene in iter_chrom(self.db.get_regions(genome_name, chrom_name),
                               start, end):
            for name2 in self.orth_lookup.get(gene.data["ID"], []):
                gene2 = self.db.get_region(name2)
                if gene2 in block_lookup:
                    frag_lookup[block_lookup[gene2]].genes.append(gene2)
                    refLookup[gene2] = gene
        self.refLookup = refLookup
        
        # determine fragment dimensions
        for frag in frag_lookup.itervalues():
            if len(frag.genes) == 0:
                frag.x = None
                continue
            frag.genes.sort(key=lambda a: a.start)
            
            # set fragment start and end
            frag.start = frag.genes[0].start
            frag.end = frag.genes[-1].end
            
            # find fragment direction
            vote = 0
            last = None
            
            for gene2 in frag.genes:
                pos = refLookup[gene2].start
                
                if last != None and pos != last:
                    if last < pos:
                        vote += 1
                    else:
                        vote -= 1
                last = pos
            
            if vote > 0:
                frag.direction = direction
            else:
                frag.direction = -direction
            
            # find fragment x-coordinate
            diffs = []
            for gene2 in frag.genes:
                if direction == 1:
                    offset1 = refLookup[gene2].start - ref_frag.start
                else:
                    offset1 = ref_frag.end - refLookup[gene2].end
                
                if frag.direction == 1:
                    offset2 = gene2.start - frag.start
                else:
                    offset2 = frag.end - gene2.end
                diffs.append(offset2 - offset1)
            frag.x = ref_frag.x - stats.median(diffs)
        
        # place blocks
        fragY = util.Dict(default=-self.genome_sep)
        for block, flip in blocks:
            frag = frag_lookup[block]
            otherGenome = frag.genome
            
            if frag.x == None:
                # fragment could not be placed
                continue
            
            frag.y = fragY[otherGenome] - \
                     ((order[otherGenome] - 1) * 
                       self.max_genome_sep)

            # re-get all genes between those coordinates
            #frag.genes = list(iter_chrom(self.db.get_regions(frag.genome, 
            #                                                 frag.chrom),
            #                             frag.start, frag.end))
            
            # store and lyaout frag
            self.frags.add(frag)
            self.layout_frag_contents(frag)

            # stagger fragments
            fragY[otherGenome] -= self.frag_sep
            if fragY[otherGenome] < -self.max_genome_sep:
                fragY[otherGenome] = -self.genome_sep
        

    def filter_blocks(self, blocks, ref_chrom, start, end):
        """Filter blocks for those that contain the chromosome.
           Also the matching side of the block is returned
        """

        for block in blocks:
            if ((block.region1.species,
                 block.region1.seqname) not in self.chroms_lookup or
                (block.region2.species,
                 block.region2.seqname) not in self.chroms_lookup):
                continue
            
            if (block.region1.seqname == ref_chrom.seqname and 
                util.overlap(block.region1.start,
                             block.region1.end,
                             start, end)):
                yield (block, 0)
            elif (block.region2.seqname == ref_chrom.seqname  and 
                  util.overlap(block.region2.start,
                               block.region2.end,
                               start, end)):
                yield (block, 1)


    def layout_frag_contents(self, frag):
        """Layout the contents of a fragment"""

        for gene in iter_chrom(self.db.get_regions(frag.genome,
                                                   frag.chrom),
                               frag.start, frag.end):
            if frag.direction == 1:
                x = frag.x + gene.start - frag.start
            else:
                x = frag.x + frag.end - gene.end
            self.region_layout[gene.data["ID"]] = Layout(x, frag.y, frag.direction)
            self.region2frags[gene.data["ID"]] = frag
        

    def draw_placed(self):
        vis = []
        
        util.tic("create draw code")
        
        # draw frags
        for frag in self.frags:
            vis.append(self.frag_widget(frag))

        # draw genes
        for reg, l in self.region_layout.iteritems():
            vis.append(translate(l.x, l.y, 
                                 self.gene_widget(self.db.get_region(reg))))
        
        # draw matches
        drawn = set()
        for frag in self.frags:
            vis.append(self.draw_matches(frag.genome, frag.chrom,
                                         frag.start, frag.end, drawn))
        
        util.toc()

        self.groupid = group(*vis)
        return self.groupid


    def draw_matches(self, sp, chrom, start, end, drawn=None):
        vis = []

        if drawn is None:
            drawn = set()
        
        # build list of matches in order of drawing
        
        for gene in iter_chrom(self.db.get_regions(sp, chrom), start, end):
            # need to sort matches by genome order so that mult-genome synteny
            # is drawn top-down

            # get orthologs
            genes2 = [x for x in self.orth_lookup.get(gene.data["ID"], [])
                      if x in self.region_layout]
            if len(genes2) == 0:
                continue
            
            rows = util.groupby(lambda x: self.region_layout[x].y, genes2)
            keys = util.sort(rows.keys(), reverse=True)
            rows = util.mget(rows, keys)

            l = self.region_layout
            
            for i in range(1, len(rows)):
                for botGene in rows[i]:
                    gene1 = self.db.get_region(botGene)
                    for topGene in rows[i-1]:

                        if (botGene, topGene) in drawn:
                            continue

                        drawn.add((botGene, topGene))
                        
                        gene2 = self.db.get_region(topGene)
                        y1 = l[topGene].y 
                        y2 = l[botGene].y + 1
                        x1 = l[topGene].x
                        x2 = l[topGene].x + gene2.length()
                        x3 = l[botGene].x + gene1.length()
                        x4 = l[botGene].x
                        
                        if self.fat_matches:
                            vis.append(quads(
                                    self.colors["matches"],
                                    x1, y1,
                                    x2, y1,
                                    x3, y2,
                                    x4, y2))

                        vis.append(lines(self.colors["matches"],
                                         x1, y1,
                                         x4, y2))
        return group(* vis)


    def frag_widget(self, frag):
        '''
        def arrow(direction, width, height, func):
            return group(
                triangles(conf['color-arrow'],
                          0, height/2,
                          direction * width, 0,
                          0, height/-2),
                hotspot("click",
                        0, height/-2,
                        direction * width, height/2,
                        func))
        
        def leftArrowFunc():
            if frag.direction == 1:
                frag.x = x - min(frag.start, conf['frag-extend'])
            else:
                frag.x = x - min(frag.chrom.size - frag.end, 
                                 conf['frag-extend'])
            if frag.direction == 1:
                frag.start = max(0, frag.start - conf['frag-extend'])
            else:
                frag.end = min(frag.chrom.size, frag.end + conf['frag-extend'])
            self.assignFragPos(conf, frag)
            self.redraw()

        def rightArrowFunc():
            if frag.direction == 1:
                frag.end = min(frag.chrom.size, frag.end + conf['frag-extend'])
            else:
                frag.start = max(0, frag.start - conf['frag-extend'])
            self.assignFragPos(conf, frag)
            self.redraw()
            '''
        
        # calculate coordinates from drawing
        x   = frag.x
        y   = frag.y
        x2  = x + frag.length()
        mid = y + .5
        top = y + 1
        vis = []
        

        # backbone
        vis.append(lines(self.colors['frag'],
                         x, mid, x2, mid,
                         x, y, x, top,
                         x2, y, x2, top))
        # hotspot
        vis.append(hotspot("click", x, y, x2, top,
                           lambda: self.frag_click(frag)))
        
        # controls
        '''
        if True: # conf['use-controls']:
            # build controls
            controls = group(
                # left arrow
                translate(x + conf['arrow-width'] / -2, y + conf['gene-size']/2,
                    arrow(-1, conf['arrow-width'], conf['arrow-height'],
                          leftArrowFunc)),
                # right arrow
                translate(x2 + conf['arrow-width'] / 2, y + conf['gene-size']/2,
                    arrow(1, conf['arrow-width'], conf['arrow-height'],
                          rightArrowFunc)))

            # add controls to controls list
            self.controlids.append(controls)
            
            # add controls to vis
            vis.append(controls)
        '''
        
        return group(* vis)
    
    #
    # gene functions
    #
    

    def gene_widget(self, gene):
        def func():
           self.gene_click(gene)

        
        eff_dir = effective_dir(self.region_layout[gene.data["ID"]].orient,
                               gene.strand)

        length = gene.length()
        if length == 0:
            return group()

        col = self.gene_color(self, gene, eff_dir)

        if eff_dir == 1:
            g = self.draw_gene(self, gene, eff_dir, col)
        else:
            g = translate(length, 0,
                          flip(0, 1, self.draw_gene(self, gene, eff_dir, col)))
        
        return group(g, hotspot("click", 0, 0, length, 1, func),
                     self.draw_gene_label(gene))
    
    
    def draw_gene_label(self, gene):

        name = self.gene_label(gene)
        length = gene.length()

        if self.show_gene_labels == "scale":
            label = group(color(0, 0, 0),
                          text_clip(name, 0, 0, length, 1,
                                    5, 20, "middle", "center"))
        elif self.show_gene_labels == "fix":
            label = group(color(0, 0, 0),
                          text(name, 0, 0, length, 1, 
                               "middle", "center"))
        elif self.show_gene_labels == "vertical":
            if gene.species == self.ref_genome:
                top = 10000000
            else:
                top = self.max_genome_sep
            label = rotate(90, color(0, 0, 0),
                           text_clip(name, 0, -100000,
                                     top, 0,
                                     self.font_size, self.font_size,
                                     "left", "top"))
        elif self.show_gene_labels == 'main_only':
            if gene.species == self.ref_genome:
                label = rotate(90, color(0, 0, 0),
                               text_clip(name, 0, -100000,
                                         10000000, 0,
                                         self.font_size, self.font_size,
                                         "left", "top"))
            else:
                label = group()
        else:
            label = group()

        return label

    
    def gene_click(self, gene):
        self.print_gene(gene)
    
    def frag_click(self, frag):

        x, y = self.win.get_mouse_pos("world")
        if frag.direction == 1:
            pos = int(frag.start + x - frag.x)
        else:
            pos = int(frag.end - (x - frag.x))
        
        print "%s:%s:%s" % (frag.genome, frag.chrom,
                            util.int2pretty(pos))

    
    def print_gene(self, gene):
        print "%s %s:%s:%s-%s (%s)" % (
               gene.data["ID"], gene.species, 
               gene.seqname, 
               util.int2pretty(gene.start), 
               util.int2pretty(gene.end), 
               util.int2pretty(gene.length()))
        print ";".join("%s=%s" % (a, b) for a,b in gene.data.items())
        print 
        
    def redraw(self):
        if self.groupid != 0:
            self.win.replace_group(self.groupid, self.draw_placed())
        else:
            self.groupid = self.win.add_group(self.draw_placed())
        self.show_controls()
    
    
    def get_gene_coords(self, gene):
        l = self.region_layout
        name = gene.data["ID"]
        return (l[name].x, l[name].y, 
                l[name].x + gene.end - gene.start,
                l[name].y + 1)
    
        
    def find(self, name):
        try:
            region = self.db.get_region(name)
            
            if region in self.region_layout: 
                self.win.set_visible(* self.get_gene_coords(region))
            else:
                print "gene '%s' is not shown" % name
        except KeyError:
            print "cannot find gene '%s'" % name


    def mark(self, name, shape="box", col=color(0, 0, 1)):
        if not (name in self.region_layout): #placedGenes):
            print "gene '%s' is not shown" % name
            return
        gene = self.db.get_region(name)
        coords = self.get_gene_coords(gene)
        
        gid = self.win.add_group(self.draw_marking(shape, col, coords[1],
                                                   coords[0], coords[2]))
        self.markids.append(gid)
    
    
    def show_marks(self, visible):
        for gid in self.markids:
            self.win.show_group(gid, visible)
    
    def clear_marks(self):
        for gid in self.markids:
            self.win.remove_group(gid)
        self.markids = []
    
    
    def show_controls(self, visible = None):
        if visible == None:
            visible = self.use_controls
        for gid in self.controlids:
            self.win.show_group(gid, visible)
    
    
    #===================================================================
    # regions
    #
    '''
    def add_regions(self, regions, shape=None, col=None, height=None):
    
        for region in regions:
            # set default visualizatios attributes
            if "shape" not in region.data:
                if shape == None:
                    region.data["shape"] = "fill"
                else:
                    region.data["shape"] = shape
            
            if "color" not in region.data:
                if col == None:
                    region.data["color"] = color(0,1,0)
                else:
                    region.data["color"] = col
            else:
                if isinstance(region.data["color"], str):
                    region.data["color"] = eval(region.data["color"])
                else:
                    region.data["color"] = region.data["color"]
            
            if "height" not in region.data:
                if height == None:
                    region.data["height"] = 1.0
                else:
                    region.data["height"] = height
            else:
                if isinstance(region.data["height"], str):
                    region.data["height"] = float(region.data["height"])
                else:
                    region.data["height"] = region.data["height"]
            
            # ensure species is specified
            assert "species" in region.data
            
            # force stand to +1 or -1
            if region.strand not in [1, -1]:
                region.strand = 1
            
            chrom = self.matching.genomes[region.data["species"]].chroms[region.seqname]
            self.regions[chrom].append(region)
        
        for lst in self.regions.itervalues():
            lst.sort(key=lambda x: x.start)
    '''
    
    # TODO: use regions as markings
        
    def draw_marking(self, shape, col, y, x1, x2, direction=1, height=1.0):
        mid = y + .5
        
        y1 = mid - height / 2.0
        y2 = mid + height / 2.0
        
        
        if shape == "box":
            return group(col, shapes.box(x1, y1, x2, y2, fill=False))
        
        elif shape == "half_box":
            if direction == 1:
                return group(col, shapes.box(x1, mid, x2, y2, fill=False))
            else:
                return group(col, shapes.box(x1, y1, x2, mid, fill=False))
        
        elif shape == "fill":
            return group(col, quads(
                x1, y1, 
                x1, y2,
                x2, y2,
                x2, y1), lines(x1, y1, x1, y2))
        
        elif shape == "half_fill":
            if direction == 1:
                return group(col, quads(
                    x1, mid, 
                    x1, y2,
                    x2, y2,
                    x2, mid), lines(x1, mid, x1, y2))
            else:
                return group(col, quads(
                    x1, y1, 
                    x1, mid,
                    x2, mid,
                    x2, y1), lines(x1, y1, x1, mid))
            
        elif shape == "cross":
            return group(lines(col, 
                x1, y1, 
                x2, y2,
                x1, y2, 
                x2, y1))
            
        elif shape == "flag":
            x = min(x1, x2)
            return group(lines(col, 
                x, y1, x, 6))
            
        elif shape == "hash":
            return group(col, lines(x1, y1, x1, y2))
        
        elif shape == "half_hash":
            if direction == 1:
                return group(col, lines(x1, mid, x1, y2))
            else:
                return group(col, lines(x1, y1, x1, mid))
            
        else:
            raise "unknown shape '%s'" % shape
    
    '''
    def drawMark(self, genome, chrom, start, end, strand=1, shape="box", col=color(0,0,1)):
        y, x1, x2 = self.getRegionDrawCoords(genome, chrom, start, end)
        
        if y == None:
            print "region not shown"
        else:
            if x1 < x2:
                direction = strand
            else:
                direction = -strand
            
            self.drawMarking(shape, col, y, x1, x2, direction)
    
    
    def getRegionDrawCoords(self, genome, chrom, start, end):
        """Returns (y, x1, x2) or (None, None, None)"""
    
        chrom = self.matching.genomes[genome].chroms[chrom]
        gene1, gene2 = self.matching.findNearGenes(chrom, start)
        gene3, gene4 = self.matching.findNearGenes(chrom, end)
        
        frags = []
        
        for gene in [gene1, gene2, gene3, gene4]:
            if gene in self.region_layout: #placedGenes:
                frags.append(self.region2frag[gene]) #gene.frag)
        
        for frag in frags:        
            if util.overlap(start, end, frag.start, frag.end):
                if frag.direction == 1:
                    return frag.y, \
                           start - frag.start + frag.x, \
                           end - frag.start + frag.x
                else:
                    return frag.y, \
                           frag.x + frag.end - start, \
                           frag.x + frag.end - end
        
        return None, None, None
    '''
    

    def toggle_controls(self):
        self.use_controls = not self.use_controls
        self.show_controls(self.use_controls)

    def toggle_labels(self):
        i = gene_label_types.index(self.show_gene_labels)
        i = (i + 1) % len(gene_label_types)
        self.show_gene_labels = gene_label_types[i]
        self.redraw()




        
        
        
#=============================================================================



# global variables
markColor = color(0,0,1)
context = 1e6
selgene = None
selgenes = []

genes = {}




##################################################################
# Classes and functions
#

class SyntenyVis (SyntenyVisBase):
    def __init__(self, genomes, chroms, regions, blocks, orths,
                 **options):
        SyntenyVisBase.__init__(self, genomes, chroms,
                             regions, blocks, orths, **options)
        self.click_mode = "gene"
        self.selgenes = []
        self.seqs = fasta.FastaDict()
                        
            
    
    def show(self):
        SyntenyVisBase.show(self)
        
        self.win.set_binding(input_key("g"), self.press("gene"))
        self.win.set_binding(input_key("v"), self.press("view"))
        self.win.set_binding(input_key("s"), self.press("sequence"))
        self.win.set_binding(input_key("a"), self.press("align"))
        self.win.set_binding(input_key("d"), self.clear_selgenes)
        self.win.set_binding(input_key("w"), self.align_selgenes)
        self.win.set_binding(input_key("e"), self.print_selgenes)

    
    def gene_click(self, gene):
        global selgene
        selgene = gene
        selgenes.append(gene)
        
        if self.click_mode == "gene":
            self.print_gene(gene)
        elif self.click_mode == "view":
            drawGene(gene.name, context)
        elif self.click_mode == "sequence":
            if gene.name in self.seqs:
                print
                print "%s\n%s" % (gene.name, self.seqs[gene.name])
                print
            else:
                print "%s has no sequence" % gene.name

        '''
        elif self.click_mode == "align":
            if gene.name not in lookup:
                print "gene %s has no matches" % gene.name
                return    
            
            orth = filter(lambda x: x in self.genes, 
                          self.comps[self.lookup[gene.name]])
            seqs2 = util.subdict(self.seqs, self.comp)
            aln = muscle.muscle(seqs2)
            
            keys = aln.keys()
            
            for key in keys:
                if not hasattr(self.genes[key], "y"):
                    self.genes[key].y = -1e1000
            
            keys.sort(lambda a,b: cmp(self.genes[b].y, self.genes[a].y))
            alignlib.printAlign(aln, order=keys)
        '''


    # add key bindings
    def press(self, mode):
        def func():
            print "mode is '%s'" % mode
            self.click_mode = mode
        return func

    def clear_selgenes(self):
        self.selgenes[:] = []
        print "selgenes cleared"

    def align_selgenes(self):
        self.align(* self.selgenes)

    def print_selgenes(self):
        print self.selgenes


    def align(self, * names):
        if len(names) == 0:
            print "nothing to align"

        # get names from genes if they are not strings
        if type(names[0]) != str:
            names = [i.name for i in names]
        
        seqs2 = util.subdict(self.seqs, names)
        aln = muscle.muscle(seqs2)
        muscle.printAlign(aln)









#=============================================================================

'''
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



def mark(shape="box", col=util.blue):
    names = []
    while True:
        line = sys.stdin.readline().rstrip()
        if line == "": break
        names.append(line)
    markGenes(names, shape, col)

def markGenes(names, shape="box", col=util.blue):
    for name in names:
        vis.mark(name, shape, col)

def markHoles(shape="box", col=util.blue):
    genes2 = filter(lambda x: len(x.matches) == 0, genes.values())
    names = [x.name for x in genes2]
    markGenes(names, "box", col)

def find(name):
    return vis.find(name)

printscreen = lambda *args, **kargs: svg.printScreen(vis.win, *args, **kargs)


def read_fasta(filename):
    vis.seqs.update(fasta.read_fasta(env.findFile(f)))

def readAllSeqs():
    util.tic("read sequences")

    for genome in m.getGenomeOrder():
        try:
            seqfile = env.findFile("%s.fasta" % genome)
            util.tic("reading '%s'" % seqfile)
            vis.seqs.read(seqfile)
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
        counts = util.hist_dict([genes[ref].chrom for ref in refs])
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
'''    


"""    
    def placeSyntenyBlock(self, conf, block, refGenome, refChrom, 
                          start, end, direction, y):
        if refGenome == block.genome2:
            otherGenome = block.genome1
            otherChrom  = block.chrom1
        else:
            otherGenome = block.genome2
            otherChrom  = block.chrom2           
        
        
        # find gene index of first gene in block
        low, startIndex = util.binsearch(refChrom.genes, start-1, 
                                        lambda a,b: cmp(a.start, b))            
        if startIndex == None:
            # quit if no gene is found
            return
        
        # interpoloate start, end, and direction to other chrom
        vote = 0
        last = None
        genes2 = {}
        for gene in GeneIter(refChrom.genes, start, end, startIndex):
            matches = self.findGeneSynteny(block, gene)
            
            if len(matches) > 0:
                pos = stats.mean(map(lambda x: x.start, matches))
                if last != None and pos != last:
                    if last < pos:
                        vote += 1
                    else:
                        vote -= 1
                last = pos
            
            for gene2 in matches:
                genes2[gene2] = 1
        genes2 = genes2.keys()        
        
        # quit if no matches are found
        if len(genes2) == 0:
            return
        
        otherStart = min(map(lambda x: x.start, genes2))
        otherEnd = max(map(lambda x: x.end, genes2))
        
        if vote > 0:
            otherDir = direction
        else:
            otherDir = -direction
        
        
                    
        # create frament
        frag = Frag()
        frag.chrom     = otherChrom 
        frag.start     = otherStart 
        frag.end       = otherEnd
        frag.direction = otherDir
        frag.x         = self.findFragX(conf, block, startIndex, refChrom, 
                            start, end, direction, otherChrom, otherStart, 
                            otherEnd, otherDir)
        frag.y         = y - ((self.order[otherGenome.name] - 1) *
                              conf['max-genome-sep'])
        
        # store frag
        self.frags[frag] = True
    
        # assign genes (x,y)-coords
        self.assignFragPos(conf, frag)
    

    
    def findOtherDir(self, conf, block, refChrom, refIndex, start, end, 
                     direction, otherChrom, otherStart, otherEnd):
        last = 0
        vote = 0
        
        # for every increase in refChrom is there a general dec or inc?
        for gene in GeneIter(refChrom.genes, start, end, refIndex):
            genes2 = self.findGeneSynteny(block, gene)
            if len(genes2) > 0:
                pos = stats.mean(map(lambda x: x.start, genes2))
                if last < pos:
                    vote += 1
                else:
                    vote -= 1
                last = pos

        # return dir based on majority vote
        if vote > 0:
            return direction
        else:
            return -direction
    
        
    
    def findFragX(self, conf, block, index1,
                  chrom1, start1, end1, direction1,
                  chrom2, start2, end2, direction2):
        diffs = []
        initStart = chrom1.genes[index1].start
        initX     = chrom1.genes[index1].x
        
        # find the average difference between start of frag and 'start' of gene
        for gene in GeneIter(chrom1.genes, start1, end1, index1):
            diff1 = gene.start - initStart
            
            for gene2 in self.findGeneSynteny(block, gene):
                diff2 = None
                if direction2 == 1:
                    diff2 = gene2.start - start2
                else:
                    diff2 = end2 - gene2.end
                diffs.append(diff2 - diff1)
        
        if len(diffs) == 0:
            return None
        else:    
            return initX - stats.mean(diffs)
    
    def findGeneSynteny(self, block, gene):
        genes2 = []
        
        if gene in self.matching.complookup:
            for gene2 in self.matching.comps[self.matching.complookup[gene]]:
                
                # find the genes that match in this block
                #   not on same chrom, and with block start and end
                if gene2.chrom != gene.chrom:
                    if block == None:
                        genes2.append(gene2)
                    elif gene2.chrom == block.chrom1:
                        if gene2.start >= block.start1 and \
                           gene2.end <= block.end1:
                            genes2.append(gene2)
                    elif gene2.chrom == block.chrom2:
                        if gene2.start >= block.start2 and \
                           gene2.end <= block.end2:
                            genes2.append(gene2)
        return genes2
    
    """
