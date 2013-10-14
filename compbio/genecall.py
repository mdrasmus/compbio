from rasmus import graph
from rasmus import util

from compbio import regionlib
from . import alignlib

    
def findFragments(regiondb, aln, overlapCutoff=.10):
    """Determine if alignment has gene fragments"""
    
    aln_genes = util.mget(regiondb.regions, aln.keys())
    nbrs = findNeighbors(regiondb, aln_genes)
    frags = []
    
    # are there any neighbors?
    if max(map(len, nbrs)) > 1:
        # do neighbors overlap in alignment?
        for nbr in nbrs:
            if len(nbr) > 1:
                aln2 = aln.get(x.data['ID'] for x in nbr)
                frags.extend(findMerges(aln2, overlapCutoff=overlapCutoff))
    return frags


def findNeighbors(regiondb, genes):
    """Determine which genes in 'genes' are neighboring genes in the 
       chromosomes 'chroms'
       
       returns a clustering of the genes in clusters of neighboring streaks
    """

    geneset = set(genes)
    neighbors = util.Dict(dim=2)    
    
    for gene in genes:
        chrom = regiondb.species[gene.species][gene.seqname]
        ind = regionlib.findRegion(chrom, gene)
        
        # look for neighboring genes on same strand
        if ind > 0:
            left = chrom[ind-1]
            if left in geneset and left.strand == gene.strand:
                neighbors[gene][left] = 1
                neighbors[left][gene] = 1
        
        if ind < len(chrom)-1:
            right = chrom[ind+1]
            if right in geneset and right.strand == gene.strand:
                neighbors[gene][right] = 1
                neighbors[right][gene] = 1
    
    comps = graph.connectedComponents(genes, lambda x: neighbors[x].keys())
    
    # sort neighbors in order of appearance along strand
    for comp in comps:
        comp.sort(key=lambda x: x.start, 
                  reverse=(comp[0].strand == -1))
    
    return comps


def findMerges(aln, overlapCutoff=.10):
    """Given several neighboring genes in a sub-alignment, determine whether 
       they are good candidates for merging
       
       overlapCutoff -- what percentage of a sequence is allowed to overlap 
                        with another and still be considered a candidate for
                        merging. This allows for merging in the presence of
                        sloppy alignments

                        ############---#-#-#----- 
                        -------#-----############ 
                        note: this is still mergebale
    """
    
    # find sequence ranges
    ranges = {}
    for name, seq in aln.items():
        lookup = alignlib.getAlignLookup(seq)
        ranges[name] = ((lookup[0], lookup[-1]))
    
    # names should already be sort in order along strand
    names = aln.keys()
    
    merges = [[names[0]]]
    for name in names[1:]:
        last = merges[-1][-1]
        
        overlapSize = max(ranges[last][1] - ranges[name][0], 0)
        totalSize = max(ranges[name][1], ranges[last][1]) - \
                    min(ranges[name][0], ranges[last][0])
        
        if overlapSize / float(totalSize) < overlapCutoff:
            merges[-1].append(name)
        else:
            merges.append([name])
    
    return filter(lambda x: len(x) > 1, merges)

            

"""
        
def findAllFragments(parts, genes, seqs):
    fraggenes = []
    fragparts = []
    
    for i, part in enumerate(parts[:500]):
        genes2 = []
        for item in part:
            if item in genes:
                genes2.append(genes[item])
        seqs2 = seqs.get(part)

        frags = findFragments(genes2, seqs2)
        if len(frags) > 0:
            util.log("%d Fragments found:" % i)
            util.log("   genes: %s", str(frags))
            fraggenes.append(frags)
            fragparts.append(i)

    return fragparts, fraggenes
    

"""


