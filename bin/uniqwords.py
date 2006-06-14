#!/usr/bin/env python

from rasmus import util
import sys

words = {}

for filename in sys.argv[1:]:
    data = util.readDelim(filename)
    for row in data:
        for word in row:
            words[word] = 1

for row in sys.stdin:
    tokens = row.rstrip().split()
    for token in tokens:
        words[token] = 1

keys = words.keys()
keys.sort()

for word in keys:
    print word
