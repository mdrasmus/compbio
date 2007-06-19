#!/usr/bin/env python

import sys
import os


for line in sys.stdin:
    line = line.rstrip()
    if os.path.exists(line):
        print line

