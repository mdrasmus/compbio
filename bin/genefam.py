#!/usr/bin/env python
# Wed Apr 25 14:22:01 EDT 2007
# Matt Rasmussen
#
# Make a more generalized format for Family IDs
#

import os
import sys

from rasmus import env
from rasmus import fasta
from rasmus import genomeutil
from rasmus import util
from rasmus import tablelib
from rasmus import genecluster


options = [
    ["", "make=", "make", "<partfile>", {"single": True}],
    ["", "to_part", "to_part", ""],
    
    ["f:", "famid=", "famid", "<famtab file>", {"single": True}],
    ["", "startid=", "startid", "<starting famid>", 
        {"single": True,
         "parser": int,
         "default": 0}]
]


conf = util.parseOptions(sys.argv, options, quit=True)


def main(conf):
    # setup paths
    env.addEnvPaths("DATAPATH")
    #env.addPaths(conf["paths"])
    
    if "make" in conf:
        # convert parts file to famid file
        parts = util.readDelim(conf["make"])
        famtab = genecluster.makeFamtab(parts, famid=conf["startid"])
        famtab.write()
    
    elif "to_part" in conf:
        # convert famid file to parts file
        
        partstab = tablelib.readTable(conf["famid"])
        
        for row in partstab:
            print row["genes"].replace(",", "\t")
        

main(conf)
