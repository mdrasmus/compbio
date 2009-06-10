"""

   Implements a fuzzy definition of synteny

"""

from itertools import chain

import rasmus
from rasmus import util
from rasmus.linked_list import LinkedList

        


def iter_windows(hits, radius):

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

    hits = util.PushIter(hits)
    
    running = [True]
    last_sp = [None]
    last_chrom = [None]

    def inner_iter():
        for hit in hits:
            if hit[0].species != last_sp[0] or hit[0].seqname != last_chrom[0]:
                last_sp[0] = hit[0].species
                last_chrom[0] = hit[0].seqname
                hits.push(hit)
                return
            yield hit

        running[0] = False
                
    list(inner_iter())
    while running[0]:
        yield inner_iter()
        
    

def find_syntenic_neighbors(hits, radius, radius2=None):
    """First regions must be sorted by chrom and start"""

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
        
        
