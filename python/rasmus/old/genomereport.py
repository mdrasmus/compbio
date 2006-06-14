from genomeutil import *
import copy
import math
import matrix
import util
import stats
import sys


def reportHomology(matching, filename):
    groups = copy.copy(matching.homology)

    out = file(filename, "w")
    groups.sort(lambda a,b: len(b.genes) - len(a.genes))

    print >>out, "# homology groups:", len(groups)

    # 1:...:1
    numorths = 0
    for group in groups:
        ones = 0
        for genome in matching.genomes.values():
            lst = filter(lambda x: x.chrom.genome == genome, group.genes)
            if len(lst) == 1:
                ones += 1
        if ones == len(matching.genomes):
            numorths += 1
    print >>out, "# of 1:..:1 ", numorths

    # singletons
    singles = {}
    for genome in matching.genomes.values():
        singles[genome] = 0

    for group in groups:
        ones = 0
        loner = None
        for genome in matching.genomes.values():
            lst = filter(lambda x: x.chrom.genome == genome, group.genes)
            if len(lst) == 1 and len(group.genes) == 1:
                ones += 1
                loner = genome
        if ones == 1:
            singles[loner] += 1
    for genome in matching.genomes.values():
        print >>out, genome.name, "singletons:", singles[genome]
    print >>out, "# of singletons ", sum(singles.values())
    
    
    # average correspondence
    mapping = {}
    for genome in matching.genomes.values():
        mapping[genome] = 0
    for group in groups:
        nums = {}
        for genome in matching.genomes.values():
            lst = filter(lambda x: x.chrom.genome == genome, group.genes)
            nums[genome] = len(lst)
        low = min(nums.values())
        if low == 0: low = 1
        for i in nums:
            mapping[i] += nums[i] / float(low)
    for i in mapping:
        mapping[i] /= float(len(groups))
    print >>out, "avg mapping:"
    for i in mapping:
        print >>out, " ", i.name, mapping[i]


    # report for each group
    for group in groups:
        print >>out, "----------------------------------"
        print >>out, "group id:", group.groupid
        print >>out, "total size:", len(group.genes)

        for genome in matching.genomes.values():
            lst = filter(lambda x: x.chrom.genome == genome, group.genes)
            print >>out, len(lst), "genes from", genome.name
        print >>out 

        for gene in group.genes:
            print >>out, "%s\t%s\t%s" % (gene.chrom.genome.name, gene.name, \
                                         gene.description)
        print >>out


def reportSynteny(matching, filename):
    blocks = copy.copy(matching.blocks)

    out = file(filename, "w")
    blocks.sort(lambda a,b: max(b.end1 - b.start1, b.end2 - b.start2) - \
                            max(a.end1 - a.start1, a.end2 - a.start2))

    print >>out, "# synteny blocks:", len(blocks)

    # width stats
    widths = []
    lengths = []    
    for block in blocks:
        widths.append(block.width())
        lengths.append(min(block.length1(), block.length2()))

    print >>out, "max width:", max(widths)
    print >>out, "avg width:", stats.mean(widths)
    print >>out, "sdev width:", stats.sdev(widths)

    # questionable synteny
    bads = []
    for i in range(len(widths)):
        if widths[i] > lengths[i] * .5:
            bads.append(blocks[i])
    print >>out, "# questionable blocks:", len(bads)
    for bad in bads:
        print >>out, "block id:", bad.blockid, ", length:", \
            min(bad.length1(), bad.length2()), ", width:", bad.width()


    # report for each block
    for block in blocks:
        holes = findBlockHoles(block)
        ngenes = len(block.getGenes())
        
        print >>out, "----------------------------------"
        print >>out, "block id:", block.blockid
        print >>out, "# genes: ", ngenes
        print >>out, "slope: ", block.getSlope()
        print >>out, "width:", block.width()
        print >>out, "holes in", block.genome1.name ,":", len(holes[0]), \
                      len(holes[0]) / float(ngenes)
        print >>out, "holes in", block.genome2.name ,":", len(holes[1]), \
                      len(holes[1]) / float(ngenes)
        print >>out, block.genome1.name, block.chrom1.name, \
                    block.start1, block.end1, block.length1()
        print >>out, block.genome2.name, block.chrom2.name, \
                    block.start2, block.end2, block.length2()
        print >>out

def findBlockHoles(block):
    holes = [[], []]
    geneIter = GeneIter(block.chrom1.genes, block.start1, block.end1)
    while geneIter.more():
        gene = geneIter.get()
        if not block.hasGene(gene):
            holes[0].append(gene)
    
    geneIter = GeneIter(block.chrom2.genes, block.start2, block.end2)
    while geneIter.more():
        gene = geneIter.get()
        if not block.hasGene(gene):
            holes[1].append(gene)
    return holes


def reportMultiBlocks(matching, refGenome, filename):
    out = file(filename, "w")
    
    def makeKey(blocks):
        keys = map(lambda x: x.otherGenome(refGenome), blocks.keys())
        keys.sort()
        return tuple(keys)
    
    # coverage stats
    coverage = {}
    for chrom in refGenome.chroms.values():
        for mblock in chrom.multiBlocks:
            key = makeKey(mblock.blocks)
            if not(key in coverage):
                coverage[key] = 0
            coverage[key] += mblock.end - mblock.start
    keys = coverage.keys()
    keys.sort(lambda a,b: coverage[b] - coverage[a])
    size = float(refGenome.size)
    
    print >>out, "multiblock: coverage"    
    for i in range(min(20, len(keys))):
        print >>out, map(lambda x: x.name, keys[i]), coverage[keys[i]] / size
    print >>out
    
    # print actual mutliblocks
    for chrom in refGenome.chroms.values():
        print >>out, refGenome.name, chrom.name        
        for mblock in chrom.multiBlocks:
            print >>out, " ", mblock.start, mblock.end, len(mblock.blocks)
        


#
# matching statistics
#
def numMatches(gene, otherGenome = None):
    degree = 0
    for match in gene.matches:
        if match.otherGene(gene).chrom.genome == otherGenome or otherGenome == None:
            degree += 1
    return degree


def getOrths(matching, genomeName1, genomeName2):
    genome1 = matching.genomes[genomeName1]
    genome2 = matching.genomes[genomeName2]
    orths = []
    
    for gene in genome1.genes:
        if numMatches(genome1.genes[gene], genome2) == 1:
            orths.append(genome1.genes[gene])
    return orths

def getOrphans(matching):
    genes = []
    
    for genome in matching.genomes.values():
        for geneName in genome.genes:
            gene = genome.genes[geneName]
            if numMatches(gene) == 0:
                genes.append(gene)
    return genes

def getRelativeScores(genome1, genome2):
    allScores = []
    for gene in genome1.genes.values():
        # collect scores for gene
        scores = []
        for match in gene.allMatches():
            if match.otherGene(gene).chrom.genome == genome2:
                scores.append(match.score)
        
        if len(scores) == 0:
            continue
        
        # normalize scores to max
        top = max(scores)
        for i in range(len(scores)):
            scores[i] /= top
        
        # add normalized scores to total list
        allScores += scores
    return allScores
        

#
# verification functions
#
def verifyMatches(matching):
    status = True
    coords = {}
    
    for match in matching.matches:
        if match.pruned:
            continue    # skip dups
    
        coord = None
        if match.genes[0].geneid < match.genes[1].geneid:
            coord = (match.genes[0].geneid, match.genes[1].geneid)
        else:
            coord = (match.genes[1].geneid, match.genes[0].geneid)
            
        if coords.has_key(coord):
            coords[coord].append(match)
            print "dups on lines ", \
                  map(lambda x: (x.genes[0].name, x.genes[1].name), coords[coord])
            print
            status = False
        else:
            coords[coord] = [match]
    
    for label in matching.labels:
        gene = matching.genomes[label[0]].genes[label[1]]
        matches = {}
        
        for match in gene.allMatches():
            if matches.has_key(match.otherGene(gene)):
                print "dup ", gene.name, match.otherGene(gene).name
                status = False
            matches[match.otherGene(gene)] = True
    return status


def verifyHomology(matching):
    # check whether a gene and homology group are connected
    def check(hgroup, other):
        if other.homology != hgroup:
            print "error in gene %s:%s" % (
                other.chrom.genome.name,
                other.name)
            return False
        if not (other in hgroup.genes):
            print "error in hgroup %d" % (hgroup.groupid)
            return False
        return True
    
    status = True
    
    for genome in matching.genomes.values():
        for geneName in genome.genes:
            gene = genome.genes[geneName]
            hgroup = gene.homology
            
            # ensure gene and homology group are connected
            if not(gene in hgroup.genes):
                print "ERROR: gene %s is not in its group" % (gene.name)
                status = False
            
            # ensure all matches are in homology group
            for match in gene.allMatches():
                if not check(hgroup, match.otherGene(gene)):
                    print "ERROR: gene ", match.otherGene(gene).name, \
                          " not in hgroup"
                    status = False
            
        # ensure all homology group genes are matches of each other
        #for hgroup in matching.homology:
        #    for gene in hgroup.genes:
        #        for gene2 in hgroup.genes:            
        #            if gene == gene2:
        #                continue
        #            if not gene2 in gene.neighbors():
        #                print "ERROR: ", gene.name, " not adj to ", gene2.name
        #                status = False
    return status
        

def syntenyCoverage(matching, refGenome):
    print "alignment coverage by at least one alignment"
    
    cov = 0
    tot = 0
    
    chroms = refGenome.chroms.values()
    chroms.sort(lambda a,b: cmp(a.name, b.name))
    
    for chrom in chroms:
        if chrom.name.find("NT") != -1:
            continue
        size = 0    
        for mblock in chrom.multiBlocks:
            if len(mblock.blocks) > 0:
                size += mblock.end - mblock.start
        
        cov += size
        tot += chrom.size
        
        print "chrom %s: %dbp out of %dbp covered (%.2f %%)" % \
            (chrom.name, size, chrom.size, size / float(chrom.size))
    
    print "total: %dbp out of %dbp covered (%.2f %%)" % \
        (cov, tot, cov / float(tot))

    
    print
    print "alignment coverage by all species"
    
    cov = 0
    tot = 0
    
    for chrom in chroms:
        if chrom.name.find("NT") != -1:
            continue
        size = 0
        for mblock in chrom.multiBlocks:
            if len(mblock.getGenomes()) == len(matching.genomes):
                size += mblock.end - mblock.start
        
        cov += size
        tot += chrom.size
        
        print "chrom %s: %dbp out of %dbp covered (%.2f %%)" % \
            (chrom.name, size, chrom.size, size / float(chrom.size))
    print "total: %dbp out of %dbp covered (%.2f %%)" % \
        (cov, tot, cov / float(tot))

    
    
    for genome in matching.genomes.values():
        if genome == refGenome:
            continue
    
        print
        print "alignment coverage by ", genome.name
        
        cov = 0
        tot = 0
        
        for chrom in chroms:
            if chrom.name.find("NT") != -1:
                continue
            
            size = 0
            for mblock in chrom.multiBlocks:
                if genome in mblock.getGenomes():
                    size += mblock.end - mblock.start
            
            cov += size
            tot += chrom.size
            
            print "chrom %s: %dbp out of %dbp covered (%.2f %%)" % \
                (chrom.name, size, chrom.size, size / float(chrom.size))
        print "total: %dbp out of %dbp covered (%.2f %%)" % \
            (cov, tot, cov / float(tot))


        
        
        
    
    
    
    

