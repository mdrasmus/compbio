#!/usr/bin/env python

import sys

from rasmus import util

from compbio import phylip
from compbio import fasta

options = [
    ["d:", "dist=", "dist", "<distance matrix>", 
        {"single": True}],
    ["l:", "label=", "label", "<alignment or label file>",
        {"single": True}],
    ["m:", "margin=", "margin", "<margin>",
        {"single": True,
         "default": 100,
         "parser": int}]
]


def main(conf):
    label, distmat = phylip.read_dist_matrix(conf["dist"])
    
    if "label" in conf:
        if conf["label"].endswith(".align") or \
           conf["label"].endswith(".aln") or \
           conf["label"].endswith(".afa"):
            label = fasta.read_fasta(conf["label"]).keys()
        else:
            label = util.read_strings(conf["label"])
    
    util.heatmap(distmat, rlabels=label, clabels=label, 
                 xmargin=conf["margin"], ymargin=conf["margin"], 
                 width=12, height=12)
    

conf = util.parseOptions(sys.argv, options, quit=True)
main(conf)

