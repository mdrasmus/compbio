#!/usr/bin/env python


import sys
import random


if len(sys.argv) < 3:
    print >>sys.stderr, "usage: sample.py <total> <sample>"
    sys.exit(1)


total, k = map(int, sys.argv[1:3])

sample = random.sample(range(total), k)
sample.sort(lambda a,b: cmp(b,a))

i = 0
for line in sys.stdin:
    if i == sample[-1]:
        print line,
        sample.pop()
        if len(sample) == 0:
            break
    i += 1
