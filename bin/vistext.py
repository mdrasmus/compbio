#!/usr/bin/python -i

import sys

import summon
from summon.core import *
from summon import shapes

from rasmus import util


options = [
    ]


conf = util.parseOptions(sys.argv, options, quit=True)

filename = conf["REST"][0]

stream = file(filename)

win = summon.Window(filename)

win.set_bgcolor(1,1,1)
win.add_group(group(color(0,0,0)))

y = 0
height = 1
lines = stream.readlines()

width = max(map(len, lines))

vis = []
for line in lines:
    vis.append(text_scale(line, 0, y-height*.1, width, y-height*.9, "left", "middle"))
    y -= height
    
win.add_group(group(*vis))    
win.add_group(group(shapes.boxStroke(0,0, width, y)))
win.home()

