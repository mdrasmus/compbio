"""

   Implements a fuzzy definition of synteny

"""

from itertools import chain

import rasmus
from rasmus import util
from rasmus.linked_list import LinkedList
from rasmus.sets import UnionFind

from compbio.regionlib import Region

from . import SyntenyBlock


def iter_windows(hits, radius):
    """Iterate through blast hits using a window with a radius in the
       query genome"""

    hits = util.PushIter(hits)
    cache = LinkedList()

    upstream = set()
    downstream = set()

    try:
        center = hits.next()
    except StopIteration:
        return


    while True:
        # discard anyone in the upstream that is not within radius distance
        for hit in list(upstream):
            if hit[0].end + radius < center[0].start:
                upstream.remove(hit)

        # populate downstream with all regions within in radius
        for hit in hits:
            if hit[0].start - radius > center[0].end:
                hits.push(hit)
                break
            downstream.add(hit)
            cache.append(hit)

        yield (center, upstream, downstream)

        # populate upstream
        upstream.add(center)
        
        # move center to next hit
        try:
            center = cache.pop_front()
        except IndexError:
            break

        # remove new center from downstream
        downstream.remove(center)
            

def print_window(center, upstream, downstream):

    print center[0]
    print 

    print "\n".join([str(x[0]) for x in upstream])
    print 
    print "\n".join([str(x[0]) for x in downstream])
        
    print "-"*70


def iter_chroms(hits):
    """
    Returns an iterator of iterators it, such that each it iterates over
    hits from the same species and chromosome.
    """

    hits = util.PushIter(hits)

    try:
        hit = hits.next()
    except StopIteration:
        # no hits to iterate
        return

    # initiate which species and chrom we are going to start with
    last_sp = [hit[0].species]
    last_chrom = [hit[0].seqname]
    hits.push(hit)

    def inner_iter(hits):
        """An iterator of hits from only one species, chromome"""
        for hit in hits:
            if hit[0].species != last_sp[0] or hit[0].seqname != last_chrom[0]:
                # if species,chrom changes, push hit back and return
                last_sp[0] = hit[0].species
                last_chrom[0] = hit[0].seqname
                hits.push(hit)
                return
            yield hit
    
    while hits.peek(None) != None:
        yield inner_iter(hits)
            
        
    

def find_syntenic_neighbors(hits, radius, radius2=None):
    """
    For each hit find the neighboring hits that are syntenic.

    hits -- iterable of tuples (region1, region2, extra)
    radius -- radius of window in query genome
    radius2 -- radius of window in subject genome (default=radius)

    hits must be sorted by query region species, chrom, and start
    """

    if radius2 is None:
        radius2 = radius

    for hits2 in iter_chroms(hits):
        for center, upstream, downstream in iter_windows(hits2, radius):

            start = center[1].start - radius2
            end = center[1].end + radius2
            syntenic = []

            for hit in chain(upstream, downstream):
                # determine which subjects are wihtin the window of 
                if (hit[1].species == center[1].species and 
                    hit[1].seqname == center[1].seqname and 
                    util.overlap(start, end, hit[1].start, hit[1].end)):
                    syntenic.append(hit)

            yield (center, syntenic)


def samedir_hits(hit1, hit2):
    dir1 = hit1[0].strand * hit1[1].strand
    dir2 = hit2[0].strand * hit2[1].strand

    if dir1 != dir2:
        return False

    if dir1 > 0:
        return ((hit2[0].end >= hit1[0].start and
                 hit2[1].end >= hit1[1].start) or
                (hit2[0].start <= hit1[0].end and
                 hit2[1].start <= hit1[1].end))
    elif dir1 < 0:
        return ((hit2[0].start <= hit1[0].end and
                 hit2[1].end >= hit1[1].start) or
                (hit2[0].end >= hit1[0].start and
                 hit2[1].start <= hit1[1].end))

    return True

        
def cluster_hits(hits, radius1, radius2=None, samedir=False):
    """
    Cluster hits using windows

    hits -- iterable of tuples (region1, region2, extra)
    radius -- radius of window in query genome
    radius2 -- radius of window in subject genome (default=radius)
    samdir -- whether or not to require genes in same direction

    hits must be sorted by query region species, chrom, and start
    """

    # connected components set
    comps = {}
    
    for hit, syntenic in find_syntenic_neighbors(hits, radius1, radius2):

        # get block of hit
        block = comps.get(hit, None)
        if block is None:
            block = UnionFind([hit])
            comps[hit] = block

        # union block with syntenic hits
        for hit2 in syntenic:
            block2 = comps.get(hit2, None)

            # check whether hits are in the same direction
            if samedir and not samedir_hits(hit, hit2):
                if hit2 not in comps:
                    comps[hit2] = UnionFind([hit2])
                continue
            
            if block2 is None:
                comps[hit2] = block
                block.add(hit2)
            else:
                block2.union(block)

    # get the set of blocks
    comps = set(b.root() for b in comps.itervalues())

    return comps


def hits2synteny_block(hits):
    """
    Create a Synteny block from a cluster of hits

    hits -- list of tuples (region1, region2, extra)

    """

    # find containing regions within each genome
    start1 = util.INF
    end1 = -util.INF
    start2 = util.INF
    end2 = -util.INF

    for hit in hits:
        a, b = hit[:2]
        start1 = min(start1, a.start)
        end1 = max(end1, a.end)
        start2 = min(start2, b.start)
        end2 = max(end2, b.end)

    return SyntenyBlock(Region(a.species, a.seqname, "synreg", start1, end1),
                        Region(b.species, b.seqname, "synreg", start2, end2),
                        data={"hits": hits})
    
    
