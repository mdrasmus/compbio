# python libs
import math

# summon libs
from summon.core import *
#from summon.simple import *
import summon
from summon import shapes

# rasmus libs
from rasmus import algorithms
from rasmus import util
from rasmus.bio import clustalw
from rasmus.bio import fasta




def initConf(conf = {}):
    "Initialize a configuration table"
       
    conf['color-gene']         = color(0,0,1, .2)
    conf['color-gene-nomatch'] = color(0,1,1,.2)
    conf['color-genome-div'] = color(0,0,0)
    conf['color-chrom-div']  = color(1,0,0,.5)
    conf['color-match']      = color(0,0,0)
    conf['color-block']      = color(0,1,0)
    conf['scale']            = .0000001

    conf['color-match-by-score'] = False
    
    return conf


def displayBindings():
    """Print the default bindings"""

    print "key bindings"
    print "------------------"
    print "middle click    display gene info"
    print "g               'gene' mode"
    print "a               'align' mode"
    print "r               show gene regions"
    print "alt+r           hide gene regions"
    print "c               toggle chrom divisions"


class Dotplot (summon.VisObject):
    """A dot plot visualization implemented for SUMMON"""

    def __init__(self, conf=None, matching=None, genomes=None):
        self.win = None
        
        self.genomes = genomes
        self.genomeOrder = []
        self.geneRegionsId = -1
        self.mode = "gene"
        self.showchrom = True
        self.seqs = fasta.FastaDict()
        
        self.initMatching(conf, matching, genomes)
    
    def initMatching(self, conf, matching, genomes=None):
        self.conf = conf
        self.matching = matching
        self.matching.autoconf(genomes)
        self.matching.genomeOrder = util.mget(self.matching.genomes,
                                                 genomes)
        
        # determine positioning and size
        pos = 0
        for genome in self.matching.genomeOrder:
            genome.position = pos
            pos += genome.size
        
            # determine chrom order and position
            pos2 = 0
            genome.chromOrder = util.mget(genome.chroms,
                                             genome.getChromOrder())
            for chrom in genome.chromOrder:
                chrom.position = pos2
                pos2 += chrom.size
        
        # determine plot size
        self.plotSize = self.getPlotSize()

    
    def show(self):
        self.win = summon.Window("dotplot")
        self.win.set_bgcolor(1, 1, 1)
        self.win.add_group(self.draw(self.conf, self.matching, self.genomes))
        self.win.home()
    
    
    def setSeqs(self, seqs):
        self.seqs = seqs
            
    
    def installBindings(self):
        # set bindings
        self.win.set_binding(input_key("g"), self.setModeFunc("gene"))        
        self.win.set_binding(input_key("a"), self.setModeFunc("align"))
        self.win.set_binding(input_key("r"), self.showGeneRegions)
        self.win.set_binding(input_key("r", "alt"), self.hideGeneRegions)
        self.win.set_binding(input_key("c"), self.toggleChrom)
    
    
    def setModeFunc(self, mode):
        def func():
            print "Dotplot mode is now '%s'" % mode
            self.mode = mode
        return func
    
    
    def showGeneRegions(self):
            # remove old regions
            if self.geneRegionsId != -1:
                self.win.remove_group(self.geneRegionsId)
            
            # add new regions
            geneRegions = self.drawGeneRegions(self.conf, self.matching)
            self.geneRegionsId = get_group_id(geneRegions)
            self.win.add_group(geneRegions)
            del geneRegions


    def hideGeneRegions(self):
        # remove old regions
        if self.geneRegionsId != -1:
            self.win.remove_group(self.geneRegionsId)


    def toggleChrom(self):
        self.showchrom = not self.showchrom    

        if not self.showchrom:
            self.conf['color-chrom-div-old'] = self.conf['color-chrom-div']
            self.conf['color-chrom-div'] = color(0,0,0,0)
        else:
            self.conf['color-chrom-div'] = self.conf['color-chrom-div-old']

        util.tic("drawing")
        self.win.clear_groups()
        self.win.add_group(self.draw(self.conf, self.matching, self.genomeOrder))
        util.toc()

    
    #
    # coordinate functions
    #
    
    def getChromPosition(self, chrom):
        """Returns the position of a chromosome"""
        
        return chrom.genome.position + chrom.position

    
    def getGeneCoords(self, gene):
        """Returns the start and end coordinates of gene"""
        
        chrom = gene.chrom
        offset = chrom.genome.position + chrom.position
        return (offset + gene.start, offset + gene.end)


    def getPlotSize(self):
        """Returns the total size of all genomes"""
        
        allSize = 0
        for genome in self.matching.genomes.values():
            allSize += genome.size
        return allSize
    
    

    
    #
    # drawing functions
    #
    
    def drawMatches(self, conf, matching):
        vis = []
        
        for match in matching.matches:
            gene1   = match.genes[0]
            gene2   = match.genes[1]
            
            coord1 = self.getGeneCoords(gene1)
            coord2 = self.getGeneCoords(gene2)
            
            # swap genes so that match is in upper triangle
            if coord1[0] > coord2[0]:
                gene1, gene2 = gene2, gene1
                coord1, coord2 = coord2, coord1
            
            # determine match color
            if "color" in dir(match):
                vis.append(match.color)
            else:
                vis.append(conf['color-match'])

            
            # draw match based on direction
            if gene1.direction == gene2.direction:
                vis.append(lines(
                    conf['scale'] * coord1[0],
                    conf['scale'] * coord2[0],
                    conf['scale'] * coord1[1],
                    conf['scale'] * coord2[1]))
            else:
                vis.append(lines(
                    conf['scale'] * coord1[0],
                    conf['scale'] * coord2[1],
                    conf['scale'] * coord1[1],
                    conf['scale'] * coord2[0]))
        
        return group(* vis)
    
        
    
    def drawSyntenyBlocks(self, conf, matching):
        vis = [conf['color-block']]
        
        for block in matching.blocks:
            if block == None:
                continue
            if block.genome1.order < block.genome2.order:
                offset1 = block.genome1.position + block.chrom1.position
                offset2 = block.genome2.position + block.chrom2.position
                vis.append(shapes.boxStroke(
                                     conf['scale'] * (offset1 + block.start1),
                                     conf['scale'] * (offset2 + block.start2),
                                     conf['scale'] * (offset1 + block.end1),
                                     conf['scale'] * (offset2 + block.end2)))
            else:
                offset2 = block.genome1.position + block.chrom1.position
                offset1 = block.genome2.position + block.chrom2.position
                vis.append(shapes.boxStroke(
                                     conf['scale'] * (offset1 + block.start2),
                                     conf['scale'] * (offset2 + block.start1),
                                     conf['scale'] * (offset1 + block.end2),
                                     conf['scale'] * (offset2 + block.end1)))
        return group(* vis)


    def drawChromBorders(self, conf, matching):
        vis = []
        allSize = conf['scale'] * self.plotSize
        
        vis.append(lines(conf['color-chrom-div'], 
                         0, allSize, allSize, allSize, 
                                  allSize, 0, allSize, allSize))
        for genome in matching.genomes.values():
            for chrom in genome.chroms.values():
                pos = conf['scale'] * (genome.position + chrom.position)
                vis.append(lines(0, pos, 
                                 allSize, pos,
                                 pos, 0,
                                 pos, allSize))
        return group(* vis)
    

    def drawGenomeBorders(self, conf, matching):
        vis = []
        allSize = conf['scale'] * self.plotSize
        
        vis.append(lines(conf['color-genome-div'], 
                         0, allSize, allSize, allSize, 
                                  allSize, 0, allSize, allSize))
        for genome in matching.genomes.values():
            vis.append(lines(0, conf['scale'] * genome.position, 
                                      allSize, conf['scale'] * genome.position,
                                      conf['scale'] * genome.position, 0,
                                      conf['scale'] * genome.position, allSize))
        
        # hotspot
        vis.append(hotspot("click", 0, 0, allSize, allSize, self.click))
        
        return group(* vis) 
    
    
    
    def drawGeneRegions(self, conf, matching):
        vis = []
        allSize = conf["scale"] * self.plotSize
        
        # get current view
        view = self.win.get_visible()
                
        for genome in matching.genomes.itervalues():
            for chrom in genome.chroms.itervalues():
                for gene in chrom.genes:
                    pos = conf['scale'] * (genome.position + chrom.position)

                    # calculate coords
                    start = pos + conf['scale'] * gene.start
                    end   = pos + conf['scale'] * gene.end

                    # only draw regions in view
                    if (view[0] <= start and view[2] >= start) or \
                       (view[0] <= end and view[2] >= end) or \
                       (view[1] <= start and view[3] >= start) or \
                       (view[1] <= end and view[3] >= end):


                        if len(gene.matches) > 0:
                            geneColor = conf["color-gene"]
                        else:
                            geneColor = conf["color-gene-nomatch"]

                        vis.extend(self.drawGeneRegion(conf, gene, geneColor, allSize))
        
        return group(*vis)
    
    
    def drawGeneRegion(self, conf, gene, geneColor, allSize = None):
        vis = []
        if allSize == None:
            allSize = conf["scale"] * self.plotSize
        
        pos = conf['scale'] * (gene.chrom.genome.position + gene.chrom.position)

        # calculate coords
        start = pos + conf['scale'] * gene.start
        end   = pos + conf['scale'] * gene.end        

        # draw shaded box
        vis.append(geneColor)
        vis.append(lines(0, start, allSize, start))
        vis.append(shapes.box(0, start, allSize, end))
        vis.append(lines(start, 0, start, allSize))
        vis.append(shapes.box(start, 0, end, allSize))

        # create hotspots
        #spots = makeHotspots(pos, allSize, gene)
        #vis.append(spots[0])
        #vis.append(spots[1])
        
        return vis
    
    
    def draw(self, conf, matching, genomes=None):
        self.installBindings()
        
        return group(self.drawMatches(conf, matching),
                     self.drawSyntenyBlocks(conf, matching),
                     self.drawChromBorders(conf, matching),
                     self.drawGenomeBorders(conf, matching))
    
    
    #
    # interactive functions
    #
    
    def click(self):
        pos = self.win.get_mouse_pos('world')
        
        if self.mode == "gene":
            self.printGeneInfo(pos[0], pos[1])
        elif self.mode == "align":
            self.printAlignment(pos[0], pos[1])
    
    
    def printGeneInfo(self, x, y):
        gene1 = self.getGeneFromPos(x / self.conf['scale'])
        gene2 = self.getGeneFromPos(y / self.conf['scale'])
        
        stat = [["SIDE", "GENOME", "CHROM", "GENE", "START", "END", "LEN"]]
        if gene1:
            stat.append(["bottom", 
                         gene1.chrom.genome.name,
                         gene1.chrom.name,
                         gene1.name,
                         gene1.start,
                         gene1.end,
                         gene1.length()])
        else:
            stat.append(["bottom", "-", "-", "-", "-", "-", "-"])
        
        if gene2:
            stat.append(["left", 
                         gene2.chrom.genome.name,
                         gene2.chrom.name,
                         gene2.name,
                         gene2.start,
                         gene2.end,
                         gene2.length()])
        else:
            stat.append(["bottom", "-", "-", "-", "-", "-", "-"])
        
        print
        util.printcols(stat, spacing=2)
    
    
    def printAlignment(self, x, y):
        gene1 = self.getGeneFromPos(x / self.conf['scale'])
        gene2 = self.getGeneFromPos(y / self.conf['scale'])
        
        if gene1 == None or gene2 == None:
            print "could not find two genes for alignment"
            return
        
        seqs2 = self.seqs.get([gene1.name, gene2.name])
        
        if len(seqs2) < 2:
            print "not enough sequences available for alignment"
            return
        
        aln = clustalw.clustalw(seqs2)
        clustalw.printAlign(aln)
        
    
    #
    # searching functions
    #

    def markGenes(self, names, markColor):
        vis = []
        allSize = self.conf["scale"] * self.plotSize
        for name in names:
            gene = self.matching.genes[name]
            vis.extend(self.drawGeneRegion(self.conf, gene, markColor, allSize))
        self.win.add_group(group(* vis))
    
    
    def getGenomeFromPos(self, pos):
        if pos < 0 or pos > self.plotSize:
            return None
        
        low, top = algorithms.binsearch(self.matching.genomeOrder, pos, 
            lambda genome, p: cmp(genome.position, p))

        return self.matching.genomeOrder[low]
    
    
    def getChromFromPos(self, pos):
        genome = self.getGenomeFromPos(pos)
        
        if genome == None:
            return None
        
        if pos < genome.position or pos > genome.position + genome.size:
            return None
        
        low, top = algorithms.binsearch(genome.chromOrder, 
                                        pos - genome.position, 
                                        lambda chrom, p: cmp(chrom.position, p))
        
        return genome.chromOrder[low]
    
    
    def getGeneFromPos(self, pos):
        chrom = self.getChromFromPos(pos)
        
        if chrom == None:
            return None
        
        if pos < chrom.genome.position + chrom.position or \
           pos > chrom.genome.position + chrom.position + chrom.size:
            return None
        
        low, top = algorithms.binsearch(chrom.genes, 
                                        pos - chrom.genome.position 
                                            - chrom.position, 
                                        lambda gene, p: cmp(gene.start, p))
        
        gene = chrom.genes[low]
        
        if chrom.genome.position + chrom.position + gene.start <= pos <= \
           chrom.genome.position + chrom.position + gene.end:    
            return gene
        else:
            return None
        

    def findBlock(self, blockid):
        block = self.matching.blocks[blockid]
        
        pos1 = block.genome1.position + block.chrom1.position
        pos2 = block.genome2.position + block.chrom2.position
        
        if block.genome1.order < block.genome2.order:
            self.win.set_visible(self.conf['scale'] * (block.start1 + pos1),
                                 self.conf['scale'] * (block.start2 + pos2),
                                 self.conf['scale'] * (block.end1 + pos1),
                                 self.conf['scale'] * (block.end2 + pos2))
        else:
            self.win.set_visible(self.conf['scale'] * (block.start2 + pos2),
                                 self.conf['scale'] * (block.start1 + pos1),
                                 self.conf['scale'] * (block.end2 + pos2),
                                 self.conf['scale'] * (block.end1 + pos1))

    def findChroms(self, genomeName1, chromName1, genomeName2, chromName2):
        pos1 = self.matching.genomes[genomeName1].position + \
               self.matching.genomes[genomeName1].chroms[chromName1].position
        size1 = self.matching.genomes[genomeName1].chroms[chromName1].size
        pos2 = self.matching.genomes[genomeName2].position + \
               self.matching.genomes[genomeName2].chroms[chromName2].position
        size2 = self.matching.genomes[genomeName2].chroms[chromName2].size
        
        if self.matching.genomes[genomeName1].order < \
           self.matching.genomes[genomeName2].order:
            self.win.set_visible(self.conf['scale'] * (pos1),
                                 self.conf['scale'] * (pos2),
                                 self.conf['scale'] * (pos1 + size1),
                                 self.conf['scale'] * (pos2 + size2))
        else:
            self.win.set_visible(self.conf['scale'] * (pos2),
                                 self.conf['scale'] * (pos1),
                                 self.conf['scale'] * (pos2 + size2),
                                 self.conf['scale'] * (pos1 + size1))

        
        
