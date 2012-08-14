




from rasmus import util

from compbio import regionlib

from . import SyntenyBlock, make_orth



def is_contig(db, genes):
    """Returns True if genes are contiguous along chromosome"""
    
    if len(genes) > 1:
        pos = [db.get_region_pos_full(i) for i in genes if i in db.regions]

        # ensure hits are on same chromosome
        if not util.equal(* util.cget(pos, 1)):
            return False

        ind = util.cget(pos, 2)
        ind.sort()

        # check that each position is present
        i = ind[0]
        for j in ind[1:]:
            if j != i+1:
                return False
            i += 1

    return True


def is_orth_contig(db, orth):
    """Returns True if both sets of genes are contiguous"""
    return is_contig(db, orth[0]) and is_contig(db, orth[1])


def orth_regions(db, orth):
    """Returns the pair of regions spanning an ortholog pair"""

    # get boundaries of spanning regions
    start1 = min(db.get_region(i).start for i in orth[0])
    end1 = max(db.get_region(i).end for i in orth[0])    
    start2 = min(db.get_region(i).start for i in orth[0])
    end2 = max(db.get_region(i).end for i in orth[0])
    
    gene1 = db.get_region(orth[0][0])
    gene2 = db.get_region(orth[1][0])

    # get strands
    if len(orth[0]) == 1:
        strand1 = gene1.strand
    else:
        strand1 = 0
    
    if len(orth[1]) == 1:
        strand2 = gene2.strand
    else:
        strand2 = 0

    # return spanning regions
    return regionlib.Region(gene1.species, gene1.seqname, "synteny",
                            start1, end1, strand1), \
           regionlib.Region(gene2.species, gene2.seqname, "synteny",
                            start2, end2, strand2)



def find_next_gene(func, genes, index, step):
    while True:            
        if index >= len(genes):
            return None
        gene = genes[index].data["ID"]
        if func(gene):
            return gene
        index += step



def can_append_orth(db, block_orth, block_dir, orth, orthdb={}):
    """Try to add an ortholog pair to the synteny block.
       returns direction of append:
           1  = parallel
           -1 = anti-parallel
           0  = no append was possible
    """

    current_orth1, current_orth2 = block_orth

    # get current2 left and right
    current_species2, current_chrom2, current_left2 = \
        db.get_region_pos_full(current_orth2[0])
    current_right2 = db.get_region_pos(current_orth2[-1])


    o1, o2 = orth
    species2, chrom2, left2 = db.get_region_pos_full(o2[0])
    right2 = db.get_region_pos(o2[-1]) #orth[-1][0])


    # get the regions on current_chrom2
    lst2 = db.get_regions(current_species2, current_chrom2)

    # must be on same species and chrom
    if species2 != current_species2 or \
       chrom2 != current_chrom2:
        return 0

    # first non-loss must be in our orth
    def is_parallel(current_right2, orth):
        gene2 = find_next_gene(lambda x: x in orthdb,
                               lst2, current_right2 + 1, 1)
        if gene2 is None or gene2 not in orth[1]:
            return 0
        else:
            return 1

    def is_antiparallel(current_left2, orth):
        # first non-loss must be in our orth
        gene2 = find_next_gene(lambda x: x in orthdb,
                               lst2, current_left2 - 1, -1)
        if gene2 is None or gene2 not in orth[1]:
            return 0
        else:
            return -1


    if block_dir == 0:
        # block has no direction yet, check both

        if left2 > current_right2:
            return is_parallel(current_right2, orth)
        elif right2 < current_left2:
            return is_antiparallel(current_left2, orth)
        else:
            return 0

    elif block_dir == 1:
        # block has direction 1

        # ortholog pair of single genes must match direction perfectly
        if len(orth[0]) == len(orth[1]) == 1:
            s1 = db.get_region(orth[0][0]).strand
            s2 = db.get_region(orth[1][0]).strand
            if s1 * s2 != 1:
                return 0

        return is_parallel(current_right2, orth)        

    elif block_dir == -1:
        # block has direction -1

        # check direction matches
        if len(orth[0]) == len(orth[1]) == 1:
            s1 = db.get_region(orth[0][0]).strand
            s2 = db.get_region(orth[1][0]).strand
            if s1 * s2 != -1:
                return 0

        return is_antiparallel(current_left2, orth)
    else:
        raise Exception("unknown direction '%d'" % block_dir)



def find_synteny(species1, species2, regions1, regions2, orths):

    # ortholog db {gene1 -> orthologs of gene1}
    orthdb = util.Dict(default=set())
    for row in orths:
        orthdb[row[0]].add(row[1])
        orthdb[row[1]].add(row[0])


    # TODO: generalize
    # inparalogs db {gene1 -> gene1a | gene1 and gene1a both have same orthologs}
    inpardb = util.Dict(default=set())
    for gene, others in orthdb.iteritems():
        inpardb[gene] = orthdb[iter(others).next()]


    # make region db
    regiondb = regionlib.RegionDb(regions1 + regions2)

    # get chromosome sets
    chroms1 = regiondb.get_chroms(species1)
    chroms2 = regiondb.get_chroms(species2)


    blocks = []
    for chname1, chrom1 in chroms1.iteritems():
        # skip empty chromosomes
        if len(chrom1) == 0:
            continue
        
        # start a new block
        need_new_block = True
        loss_streak = []
        for i, gene1 in enumerate(chrom1):
            names2 = orthdb[gene1.data["ID"]]

            # no orthologs, start a loss streak
            if len(names2) == 0:
                #need_new_block = True
                loss_streak.append(make_orth(regiondb, [gene1.data["ID"]], []))
                continue

            # make ortholog cluster
            names1 = inpardb[gene1.data["ID"]]            
            orth = make_orth(regiondb, names1, names2)
            
            # orthologs are not contiguous, stop block
            if not is_orth_contig(regiondb, orth):
                loss_streak = []
                need_new_block = True
                continue

            # try to add to existing block
            if not need_new_block:
                block = blocks[-1]

                # just continue if we are still in the last ortholog pair
                #   i.e. gene1 is a paralog in a tandem set
                if orth == block.orths[-1]:
                    continue

                # try to append
                direction = can_append_orth(regiondb, block.orths[-1],
                                            block.dir, orth, orthdb)
                if direction == 0:
                    loss_streak = []
                    need_new_block = True
                else:
                    for loss in loss_streak:
                        block.add_orth(loss, direction)
                    loss_streak = []
                    block.add_orth(orth, direction)

            # start a new block
            if need_new_block:
                loss_streak = []
                if len(blocks) > 0:
                    blocks[-1].recalc_regions(regiondb)
                blocks.append(SyntenyBlock(* orth_regions(regiondb, orth)))
                blocks[-1].add_orth(orth)
                need_new_block = False

        if len(blocks) > 0:
            blocks[-1].recalc_regions(regiondb)

    return blocks
