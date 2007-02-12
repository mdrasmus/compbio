#!/usr/bin/env python
#

"""
COORD format
<gene name> <chrom name> <start> <end> <strand>


GFF format
<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

AKA

<chrom> . gene <start> <end> . <strand> 0 ID=<gene>


"""

import sys

if (sys.argv) < 3:
    print >>sys.stderr, "coord2gff.py <species> <coord file>"
    sys.exit(1)

species = sys.argv[1]

for line in file(sys.argv[2]):
    gene, chrom, start, end, strand = line.rstrip().split("\t")
    
    if strand == "1":
        strand = "+"
    else:
        strand = "-"
    
    print "\t".join([chrom, ".", "gene", start, end, ".", strand, "0", 
                     "ID=%s;species=%s;" % (gene, species)])
    
