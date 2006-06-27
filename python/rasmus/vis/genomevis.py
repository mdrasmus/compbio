
# python libs
import copy

# rasmus libs
from rasmus import gff
from rasmus import fff
from rasmus import stats
from rasmus import util
from rasmus.genomeutil import *


# graphics libs
from summon import *
from summonlib import shapes


# globals
defaultVis = None
geneLabelTypes = [False, "fix", "scale", "vertical"]



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


def toggle_controls():
    if defaultVis != None:
        defaultVis.conf['use-controls'] = not defaultVis.conf['use-controls']
        defaultVis.showControls(defaultVis.conf['use-controls'])

def toggle_labels():
    if defaultVis != None:
        i = geneLabelTypes.index(defaultVis.conf['show-gene-labels'])
        i = (i + 1) % len(geneLabelTypes)
        defaultVis.conf['show-gene-labels'] = geneLabelTypes[i]
        defaultVis.redraw()


set_binding(input_key("c"), toggle_controls)
set_binding(input_key("l"), toggle_labels)





def invcolor(c, alpha = 1):
    return color(1 - c[0], 1 - c[1], 1 - c[2], alpha)

   

def initConf(conf = None, size = 1000):
    "initialize a configuration table"
    
    if conf == None:
        conf = {}
    
    bgcolor = get_bgcolor()
    
    setConfScale(conf, size)
    conf['color-genome-div'] = invcolor(bgcolor, .5)
    conf['color-gene-pos'] =  color(1, .6, 0, .95)
    conf['color-gene-neg'] =  color(1, .6, 0, .95)
    conf['color-gene2-pos'] =  color(1, .2, 0, .95)
    conf['color-gene2-neg'] =  color(1, .2, 0, .95)
    conf['color-matches']  =  color(1, 1, 0, .9)
    conf['color-arrow']    =  color(1, .8, 0, .5)
    conf['color-frag']     =  invcolor(bgcolor, .8)
    conf['color-blocks']   = [color(.8,.8,1,.5)]
    
    conf['fat-matches']    = True
    conf['use-controls']   = False
    conf['min-frag']       = 0
    conf['show-gene-labels'] = False
    
    return conf


def setLightColors(conf = None):
    if conf == None:
        conf = {}
        
    bgcolor = get_bgcolor()
    
    conf['color-genome-div'] = invcolor(bgcolor, .5)
    conf['color-gene-pos'] =  color(1, .9, .5, .95)
    conf['color-gene-neg'] =  color(1, .5, .9, .95)
    conf['color-gene2-pos'] =  color(1, .6, 0, .95)
    conf['color-gene2-neg'] =  color(1, .6, 0, .95)    
    conf['color-matches']  =  color(1, 1, 0, .9)
    conf['color-arrow']    =  color(1, .8, 0, .5)
    conf['color-frag']     =  invcolor(bgcolor, .8)
    conf['color-blocks']   = [color(1, 0, 0, .6),
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
    
    base = .5
    for i in xrange(len(conf["color-blocks"])):
        c = conf["color-blocks"][i]
        newcolor = []
        for j in [1,2,3]:
            newcolor.append((c[j] * (1 - base)) + base)
        newcolor.append(.6)
        conf["color-blocks"][i] = color(* newcolor)
    
    return conf


def calcGeneHeight(matching):
    genes = matching.getGenes().values()
    lens = []
    for gene in genes:
        lens.append(gene.end - gene.start + 1)
    meanlen = sum(lens) / float(len(lens))
    
    return meanlen / 2.0


def setConfScale(conf, size):
    conf['gene-size']      = size
    conf['genome-sep']     = 3 * size
    conf['max-genome-sep'] = 6 * size
    conf['frag-sep']       = 1.5 * size
    conf['frag-extend']    = 50 * size
    conf['arrow-width']    =  1.5 * size
    conf['arrow-height']   =  2 * size


def drawChromLegend(conf):
    vis = []
    i = 1
    for c in conf['color-blocks'][1:]:
        vis.append(c)
        vis.append(translate(i, 0,
            shapes.box(0, 0, .9, .7),
            color(0,0,0),
            text_scale(str(i), .1, 0, .8, -.5, "center", "bottom")))
        i += 1
    vis.append(color(0,0,0))
    vis.append(text_scale("chromosome colors", 1, 1, i+1, 2, "center"))
    return list2group(vis)



def effectiveDir(direction, geneDirection):
    return (direction == 1) == (geneDirection == 1)


#
# main visualization 
#
class Frag:
    "represents fragment of DNA"
    
    def __init__(self):
        self.chrom = None
        self.start = 0
        self.end = 0
        self.direction = 0
        self.x = 0
        self.y = 0
        self.genes = []

    

class SyntenyVis:    
    def __init__(self, conf, matching, rootid=get_root_id()):
        self.conf = None
        self.matching = None
        self.genes         = {}
        self.frags         = {}
        self.placedGenes   = {}
        self.placedFrags   = {}
        self.placedRegions = {}
        self.controlids    = []
        self.markids       = []
        self.labelids      = []
        self.order         = {}
        self.groupid       = 0
        
        self.conf = conf
        self.matching = matching
        self.genes = matching.getGenes()
        self.regions = util.Dict(1, [])
        self.findex = fff.FeatureIndex()
        
        self.rootid = rootid
        self.visid = insert_group(self.rootid, group())
        
        
    def clearDrawing(self):
        remove_group(self.visid)
        self.clearMarks()
        self.visid = insert_group(self.rootid, group())
        
    
    
    def draw(self, refGenomeName, refChromName, start, end, direction=1):
        self.clearDrawing()
        insert_group(self.visid, 
            self.drawChromosome(refGenomeName, 
                                refChromName, 
                                start, end, direction=direction))
        self.showControls(self.conf["use-controls"])
    
    
    def drawAll(self, refGenomeName):
        self.clearDrawing()
        
        y = 0
        chroms = self.matching.genomes[refGenomeName].chroms.values()
        chroms.sort(lambda a,b: cmp(b.size, a.size))
        
        for chrom in chroms:
            util.tic("drawing chrom %s" % chrom.name)
            insert_group(self.visid, group(translate(0, y, 
                self.drawChromosome(refGenomeName, chrom.name, 0, chrom.size))))
            self.showControls(self.conf["use-controls"])
            util.toc()
            y -= self.conf['max-genome-sep'] * (len(self.matching.genomes) + 1)
    
    
    def drawChromosome(self, refGenomeName, refChromName, start, end, direction=1):
        refGenome = self.matching.genomes[refGenomeName]
        refChrom  = refGenome.chroms[refChromName]
    
        # init reference fragment
        refFrag = Frag()
        refFrag.chrom     = refChrom
        refFrag.start     = max(start,0)
        refFrag.end       = min(end, refChrom.size)
        refFrag.direction = direction
        refFrag.x         = max(start,0)
        refFrag.y         = 0


        # init visualization
        self.frags = { refFrag : True }
        self.placedGenes = {}
        self.placedFrags = []
        self.placedRegions = {}
        self.controlids = []
        self.order = {}
        self.start = start
        self.end = end
        
        
        # swap the genome with order 0 and the reference genome
        for genome in self.matching.genomes.values():
            if genome.order == 0:
                self.order[genome.name] = refGenome.order
            elif genome == refGenome:
                self.order[genome.name] = 0
            else:
                self.order[genome.name] = genome.order
        
        
        util.tic("placing fragments")
        
        
        # store ref chrom as seed for matches
        self.placedFrags.append(refFrag)
        
        # assign genes (x,y)-coords from reference fragment
        self.assignFragPos(self.conf, refFrag)
        
        
        # find all synteny blocks in this region
        drawBlocks = []
        for block in self.matching.blocks:
            if block.chrom1 == refChrom:
                drawBlocks.append((block, 0))
            elif block.chrom2 == refChrom:
                drawBlocks.append((block, 1))
        
        # sort blocks by appearance in refChrom
        def blocksort(a, b):
            if a[1] == 0:
                starta = a[0].start1
            else:
                starta = a[0].start2
            
            if b[1] == 0:
                startb = b[0].start1
            else:
                startb = b[0].start2
            return starta - startb
        drawBlocks.sort(cmp=blocksort)
        
        
        # make lookup for genes to block and block to fragment
        blockLookup = {}
        fragLookup = {}
        for block, flip in drawBlocks:
            if flip == 0:
                otherChrom = block.chrom2
                otherStart = block.start2
                otherEnd = block.end2
            else:
                otherChrom = block.chrom1
                otherStart = block.start1
                otherEnd = block.end1
                
            frag = Frag()
            frag.chrom = otherChrom
            fragLookup[block] = frag
            
            for gene2 in GeneIter(otherChrom.genes, otherStart, otherEnd):
                blockLookup[gene2] = block

        
        
        # find all genes that will be drawn
        # walk along refChrom and store drawn genes into fragments
        refLookup = {}
        for gene in GeneIter(refChrom.genes, start, end):
            genes2 = self.matching.getComp(gene)
            for gene2 in genes2:
                if gene2 in blockLookup:
                    fragLookup[blockLookup[gene2]].genes.append(gene2)
                    refLookup[gene2] = gene
        
        
        # determine fragment dimensions
        for frag in fragLookup.values():
            if len(frag.genes) == 0:
                frag.x = None
                continue
            frag.genes.sort(lambda a,b: a.start - b.start)
            
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
                    offset1 = refLookup[gene2].start - refFrag.start
                else:
                    offset1 = refFrag.end - refLookup[gene2].end
                
                if frag.direction == 1:
                    offset2 = gene2.start - frag.start
                else:
                    offset2 = frag.end - gene2.end
                diffs.append(offset2 - offset1)
            frag.x = refFrag.x - stats.median(diffs)
        
        # place blocks
        fragY = util.Dict(1, -self.conf['genome-sep'])
        for block, flip in drawBlocks:
            frag = fragLookup[block]
            otherGenome = frag.chrom.genome
            
            if frag.x == None:
                # fragment could not be placed
                continue
            
            frag.y = fragY[otherGenome] - \
                     ((self.order[otherGenome.name] - 1) * 
                       self.conf['max-genome-sep'])
            
            # store frag
            self.frags[frag] = 1

            # assign genes (x,y)-coords
            self.assignFragPos(self.conf, frag)

            # stagger fragments
            fragY[otherGenome] -= self.conf['frag-sep']
            if fragY[otherGenome] < -self.conf['max-genome-sep']:
                fragY[otherGenome] = -self.conf['genome-sep']
        
        util.toc()
        
        # draw placed objects
        return self.drawPlaced()
    

    def assignFragPos(self, conf, frag):
        for gene in GeneIter(frag.chrom.genes, frag.start, frag.end):
            if frag.direction == 1:
                gene.x = frag.x + gene.start - frag.start
            else:
                gene.x = frag.x + frag.end - gene.end
            gene.y = frag.y
            gene.visDir = frag.direction
            gene.frag = frag
            
            # record gene as "placed"
            self.placedGenes[gene] = True
        
        
        for region in gff.RegionIter(self.regions[frag.chrom], frag.start, frag.end):
            if frag.direction == 1:
                region.x = frag.x + region.start - frag.start
            else:
                region.x = frag.x + frag.end - region.end
            region.y = frag.y
            region.direction = frag.direction
            region.frag = frag
            
            # record region as "placed"
            self.placedRegions[region] = True


    def drawPlaced(self):
        global defaultVis
        
        vis = []
        
        util.tic("create draw code")

        defaultVis = self
        
        # draw frags
        for frag in self.frags:
            vis.append(self.fragWidget(self.conf, frag))

        # draw genes
        for gene in self.placedGenes:
            vis.append(translate(gene.x, gene.y, 
                                 self.geneWidget(self.conf, gene)))
        
        # draw regions
        for region in self.placedRegions:
            vis.append(self.drawRegion(region))
            
        for frag in self.placedFrags:
            vis.append(self.drawFrequentFeatures(frag))
        
        
        # draw matches
        for frag in self.placedFrags:
            vis.append(self.drawMatches(self.conf,
                                     frag.chrom, frag.start, frag.end))
        
        util.toc()

        g = list2group(vis)
        self.groupid = get_group_id(g)
        return g


    def drawMatches(self, conf, chrom, start, end):
        def getBlockColor(chrom):
            if not chrom.isdigit():
                num = 0
            else:
                num = int(chrom)
            return conf['color-blocks'][num % len(conf['color-blocks'])]

        vis = []
        
        # build list of matches in order of drawing
        for gene in GeneIter(chrom.genes, start, end):
            # need to sort matches by genome order so that mult-genome synteny
            # is drawn top-down
            
            #genes2 = self.findGeneSynteny(None, gene)
            
            if gene in self.matching.complookup:
                genes2 = self.matching.comps[self.matching.complookup[gene]]
            else:
                continue
            
            genes2 = filter(lambda gene: "y" in dir(gene) and
                                         gene in self.placedGenes, genes2)
            
            rows = util.groupby(lambda x: x.y, genes2)
            keys = util.sort(rows.keys(), compare=util.invcmp)
            rows = util.sublist(rows, keys)
            
            for i in range(1, len(rows)):
                for botGene in rows[i]:
                    for topGene in rows[i-1]:
                        y1 = topGene.y 
                        y2 = botGene.y + conf['gene-size']
                        x1 = topGene.x
                        x2 = topGene.x + topGene.end - topGene.start
                        x3 = botGene.x + botGene.end - botGene.start
                        x4 = botGene.x

                        if conf['fat-matches']:
                            vis.append(quads(
                                getBlockColor(botGene.chrom.name),
                                vertices(
                                  x1, y1,
                                  x2, y1,
                                  x3, y2,
                                  x4, y2)))

                        vis.append(lines(
                            getBlockColor(botGene.chrom.name),
                            vertices(
                              x1, y1,
                              x4, y2)))
        return group(* vis)


    def fragWidget(self, conf, frag):
        def arrow(direction, width, height, func):
            return group(
                triangles(conf['color-arrow'],
                    vertices(0, height/2,
                             direction * width, 0,
                             0, height/-2)),
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
        
        # calculate coordinates from drawing
        x   = frag.x
        y   = frag.y
        x2  = x - frag.start + min(frag.end, frag.chrom.size)
        mid = y + conf['gene-size'] / 2
        top = y + 1.2 * conf['gene-size']
        vis = []
        

        # backbone
        vis.append(lines(conf['color-frag'],
                   vertices(x, mid, x2, mid,
                            x, y, x, top,
                            x2, y, x2, top)))
        
        
        # controls
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
            self.controlids.append(get_group_id(controls))
            
            # add controls to vis
            vis.append(controls)
        
        return group(* vis)
    
    #
    # gene functions
    #
    def drawGene(self, conf, length, direction, order):
        "draw a single gene"
        height = conf['gene-size']
        steep = .1
        
        if direction == 1:
            if order % 2 == 0:
                col = conf['color-gene-pos']
            else:
                col = conf['color-gene2-pos']
        
            mid = max(length - (height * steep), 0)
            return group(polygon(col, vertices( 
                0, 0, mid, 0, length, height/2, mid, height, 0, height)),
                lines(0, 0, 0, height))
        else:
            if order % 2 == 0:
                col = conf['color-gene-neg']
            else:
                col = conf['color-gene2-neg']
        
            mid = min(height * steep, length)
            return group(polygon(col, vertices( 
                length, 0, length, height, mid, height, 0, height/2, mid, 0)),
                lines(length, 0, length, height))


    def geneWidget(self, conf, gene):
        def func():
           self.geneClick(gene)

        length = gene.end - gene.start
        effDir = effectiveDir(gene.visDir, gene.direction)

        if conf['show-gene-labels'] == "scale":
            label = text_clip(gene.name, 0, 0, length, conf['gene-size'], 
                         5, 20, "middle", "center")
        elif conf['show-gene-labels'] == "fix":
            label = text(gene.name, 0, 0, length, conf['gene-size'], 
                         "middle", "center")
        elif conf['show-gene-labels'] == "vertical":
            label = rotate(90, text_clip(gene.name, 0, -conf["gene-size"],
                                          conf['max-genome-sep'] / 2.0, 0,
                                          5, 20,
                                          "left", "top"))
        else:
            label = group()
        
        # determine which row gene is in (for color)
        order = self.order[gene.chrom.genome.name]
    
        return group(
            self.drawGene(conf, length, effDir, order),
            hotspot("click", 0, 0, length, conf['gene-size'], func),
            color(0,0,0),
            label)

    
    def geneClick(self, gene):
        self.printGene(gene)
    
    
    def printGene(self, gene):
        print "%s:%s:%s (%s-%s %s)" % (
               gene.chrom.genome.name,
               gene.chrom.name, gene. name, \
               util.int2pretty(gene.start), 
               util.int2pretty(gene.end), 
               util.int2pretty(gene.end - gene.start))
        
    def redraw(self):
        if self.groupid != 0:
            replace_group(self.groupid, self.drawPlaced())
        else:
            self.groupid = add_group(self.drawPlaced())
        self.showControls()
    
    
    def getGeneCoords(self, gene):
        return (gene.x, gene.y, 
                gene.x + gene.end - gene.start, gene.y + self.conf['gene-size'])
    
        
    def find(self, name):
        if name in self.genes:
            gene = self.genes[name]
            if gene in self.placedGenes:
                apply(set_visible, self.getGeneCoords(gene))
            else:
                print "gene '%s' is not shown" % name
        else:
            print "cannot find gene '%s'" % name
    
    
    def mark(self, name, shape="box", col=color(1, 1, 0)):
        if name in self.genes:
            gene = self.genes[name]
            self.markGene(gene, shape, col)
        else:
            print "cannot find gene '%s'" % name


    def markGene(self, gene, shape="box", col=color(0, 0, 1)):
        if not (gene in self.placedGenes):
            print "gene '%s' is not shown" % gene.name
            return
        coords = self.getGeneCoords(gene)
        
        gid = add_group(self.drawMarking(shape, col, coords[1], coords[0], coords[2]))
        self.markids.append(gid)
    
    
    def showMarks(self, visible):
        for gid in self.markids:
            show_group(gid, visible)
    
    def clearMarks(self):
        for gid in self.markids:
            remove_group(gid)
        self.markids = []
    
    
    def showControls(self, visible = None):
        if visible == None:
            visible = self.conf['use-controls']
        for gid in self.controlids:
            show_group(gid, visible)
    
    
    #
    # regions
    #
    
    def addRegions(self, regions, shape=None, col=None, height=None):
    
        for region in regions:
            # set default visualizatios attributes
            if "shape" not in region.attrs:
                if shape == None:
                    region.attrs["shape"] = "fill"
                else:
                    region.attrs["shape"] = shape
            
            if "color" not in region.attrs:
                if col == None:
                    region.attrs["color"] = color(0,1,0)
                else:
                    region.attrs["color"] = col
            else:
                region.attrs["color"] = eval(region.attrs["color"])
            
            if "height" not in region.attrs:
                if height == None:
                    region.attrs["height"] = 1.0
                else:
                    region.attrs["height"] = height
            else:
                region.attrs["height"] = float(region.attrs["height"])
            
            # ensure species is specified
            assert "species" in region.attrs
            
            # force stand to +1 or -1
            if region.strand not in [1, -1]:
                region.strand = 1
            
            chrom = self.matching.genomes[region.attrs["species"]].chroms[region.seqname]
            self.regions[chrom].append(region)
        
        for lst in self.regions.itervalues():
            lst.sort(key=lambda x: x.start)
    
    
    
    def addFff(self, filename):
        self.findex.read(filename)
    
    
    def drawRegion(self, region):
        
        return self.drawMarking(region.attrs["shape"], region.attrs["color"], 
                                region.y, region.x, 
                                region.x + region.end - region.start,
                                region.direction * region.strand,
                                region.attrs["height"])
    
    
    def drawFrequentFeatures(self, frag):
        
        colors = [color(1, 0, 0, .6),
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
        
        features = self.findex.getFeatures(frag.chrom.genome.name,
                                           frag.chrom.name,
                                           frag.start,
                                           frag.end)
        
        vis = []
        i = 0
        for feature, sites in features.iteritems():
            col = colors[i]
            length = self.findex.features[feature].length
            
            for site in sites:
                if frag.direction == 1:
                    x1 = frag.x + abs(site) - frag.start
                else:
                    x1 = frag.x + frag.end - abs(site) - length
                
                vis.append(self.drawMarking(
                                "half_fill", col, 
                                frag.y, x1, 
                                x1 + length,
                                util.sign(site) * frag.direction,
                                1.5))
            
            i = (i + 1) % len(colors)
            
        
        return group(* vis)
    
    
    def drawMarking(self, shape, col, y, x1, x2, direction=1, height=1.0):
        mid = y + self.conf["gene-size"] / 2.0
        
        y1 = mid - height * self.conf["gene-size"] / 2.0
        y2 = mid + height * self.conf["gene-size"] / 2.0
        
        
        if shape == "box":
            return group(col, shapes.boxStroke(x1, y1, x2, y2))
        
        elif shape == "half_box":
            if direction == 1:
                return group(col, shapes.boxStroke(x1, mid, x2, y2))
            else:
                return group(col, shapes.boxStroke(x1, y1, x2, mid))
        
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
            return group(lines(col, vertices(
                x1, y1, 
                x2, y2,
                x1, y2, 
                x2, y1)))
            
        elif shape == "flag":
            x = min(x1, x2)
            return group(lines(col, vertices(
                x, y1, x, self.conf["gene-size"] * 6)))
            
        elif shape == "hash":
            return group(col, lines(x1, y1, x1, y2))
        
        elif shape == "half_hash":
            if direction == 1:
                return group(col, lines(x1, mid, x1, y2))
            else:
                return group(col, lines(x1, y1, x1, mid))
            
        else:
            raise "unknown shape '%s'" % shape
    
    
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
            if gene in self.placedGenes:
                frags.append(gene.frag)
        
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
        low, startIndex = algorithms.binsearch(refChrom.genes, start-1, 
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

        
        
        
