#!/usr/bin/env python

import os, sys

from rasmus import util


options = [
    ["r:", "rows=", "rows", "<number of row to align at a time>",
     {"single": True,
      "default": 100,
      "parser": int}],
    ["s:", "spacing=", "spacing", "<column spacing>",
     {"single": True,
      "default": 2,
      "parser": int}],
    ["d:", "delim=", "delim", "<delimiter>",
     {"single": True,
      "default": "\t"}],
]

conf = util.parseOptions(sys.argv, options, quit=True)


if len(conf["REST"]) > 0 and conf["REST"][0] != "-":
    infile = file(conf["REST"][0])
else:
    infile = sys.stdin


mat = []
for line in infile:
    tokens = line.rstrip().split(conf["delim"])
    mat.append(tokens)
    
    if len(mat) >= conf["rows"]:
        util.printcols(mat, spacing=conf["spacing"])
        mat = []

util.printcols(mat, spacing=conf["spacing"])
    
