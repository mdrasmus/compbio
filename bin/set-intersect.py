#!/usr/bin/env python

from rasmus import util
import sys

for filename in sys.argv[1:2]:
    set = util.makeset(util.readStrings(filename))

for filename in sys.argv[2:]:
    set2 = util.makeset(util.readStrings(filename))
    set = util.intersect(set, set2)

keys = set.keys()
keys.sort()

for key in keys:
    print key
