#!/usr/bin/env python


import os
import sys


path = sys.argv[1]
ext = sys.argv[2]


infile = os.popen("ls -U '%s'" % path)

for line in infile:
    line = line.rstrip()
    print os.path.join(path, line, line + ext)
