#!/usr/bin/python

import sys
import random

random.seed()

nsample = int(sys.argv[1])

sample = []
for i in range(nsample):
    sample.append("")

lineno = 0
for line in sys.stdin:
    if (lineno < nsample):
        sample[lineno] = line
        lineno += 1
    else:
        break

for line in sys.stdin:
    lineno += 1
    if random.random() < (1.0 / lineno):
        sample[int(random.random() * nsample)] = line

for line in sample:
    print line

