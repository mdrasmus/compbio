#!/usr/bin/env python

import os, sys

from rasmus import util, tablelib


options = [
    ["r:", "rows=", "rows", "<number of row to align at a time>",
     {"single": True,
      "default": 100,
      "parser": int,
      "help": "use -1 to align entire file"}],
    ["s:", "spacing=", "spacing", "<column spacing>",
     {"single": True,
      "default": 2,
      "parser": int}],
    ["w:", "maxwidth=", "maxwidth", "<# chars>",
     {"single": True,
      "default": util.INF,
      "parser": int}],
    ["d:", "delim=", "delim", "<delimiter>",
     {"single": True,
      "default": "\t"}],
    ["H:", "header=", "header", "",
     {"single": True,
      "parser": int,
      "default": 0}],
    ["t", "table", "table", "",
     {"single": True}]
]

conf = util.parseOptions(sys.argv, options, quit=True)


if len(conf["REST"]) > 0 and conf["REST"][0] != "-":
    infile = file(conf["REST"][0])
else:
    infile = sys.stdin


if conf["table"]:
    tab = tablelib.readTable(infile)
    tab.writePretty()
else:
    mat = []
    headers = []
    
    for i, line in enumerate(infile):
        tokens = line.rstrip().split(conf["delim"])
        mat.append(tokens)

        if i < conf["header"]:
            headers.append(tokens)

        if conf["rows"] > 0 and len(mat) >= conf["rows"]:
            if len(headers) > 0 and i > len(mat):
                mat = headers + mat
            
            util.printcols(mat, spacing=conf["spacing"], colwidth=conf["maxwidth"])
            mat = []

    util.printcols(mat, spacing=conf["spacing"], colwidth=conf["maxwidth"])

