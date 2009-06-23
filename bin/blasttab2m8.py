#!/usr/bin/env python

import sys, os

from rasmus import util


# NCBI BLASTALL 2.2.10 -m 8 tab-delimited output
# Fields: 
# 0. Query id, 
# 1. Subject id, 
# 2. % identity, 
# 3. alignment length, 
# 4. mismatches, 
# 5. gap openings, 
# 6. q. start, 
# 7. q. end, 
# 8. s. start, 
# 9. s. end, 
# 10. e-value, 
# 11. bit score


# Manolis's tab-delimited blast format
# 0. query
# 1. target
# 2. idperc
# 3. length
# 4. qstart
# 5. qend
# 6. tstart
# 7. tend
# 8. gap_open
# 9. mismatches
# 10. subhits
# 11. score
# 12. probability



options = [
    ["o:", "oldext=", "oldext", "<old extension>",
        {"single": True}],
    ["n:", "newext=", "newext", "<new extension>",
        {"single": True}],
    
]

conf = util.parseOptions(sys.argv, options, quit=True, 
                         resthelp="<blast tab files> ...")


conversion = [0, 1, 2, 3, 9, 8, 4, 5, 6, 7, 12, 11]


for tabfile in conf["REST"]:
    m8file = util.replace_ext(tabfile, conf["oldext"], conf["newext"])
    
    print tabfile, "==>", m8file
    
    infile = file(tabfile)
    out = file(m8file, "w")
    
    try:
        infile.next()
        infile.next()
    
        for line in infile:
            fields = line.rstrip().split("\t")
        
            print >>out, "\t".join(util.mget(fields, conversion))
    except StopIteration:
        pass


