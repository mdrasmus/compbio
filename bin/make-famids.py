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


options = [
    ["p:", "part=", "part", "<part file>", {"single": True}],
    ["f:", "famid=", "famid", "<famid file>", {"single": True}],
]


conf = util.parseOptions(sys.argv, options, quit=True)


def main(conf):
    # setup paths
    env.addEnvPaths("DATAPATH")
    #env.addPaths(conf["paths"])
    
    if "part" in conf:
        # convert parts file to famid file
        parts = util.readDelim(conf["part"])

        partstab = tablelib.Table(headers=["famid", "genes"])

        for i, part in enumerate(parts):
            partstab.add(famid=str(i),
                         genes=",".join(part))

        partstab.write()
    
    elif "famid" in conf:
        # convert famid file to parts file
        
        partstab = tablelib.readTable(conf["famid"])
        
        for row in partstab:
            print row["genes"].replace(",", "\t")
        

main(conf)
