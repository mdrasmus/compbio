#!/usr/bin/env python

import os, sys
from rasmus import util



options = [
    ["p:", "part=", "part", "<part file>"],
]



# check args
try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


parts = util.readDelim(param["part"][0])
lookup = {}
for i in xrange(len(parts)):
    for word in parts[i]:
        lookup[word] = i

for fn in param["part"][1:]:
    for line in file(fn):
        words = line.rstrip().split()
        for i in xrange(len(words)):
            for j in xrange(i+1, len(words)):
                if words[i] in lookup and words[j] in lookup:
                    if lookup[words[i]] != lookup[words[j]]:
                        print words[i], words[j]


