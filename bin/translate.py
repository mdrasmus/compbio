#!/usr/bin/env python

import sys
from rasmus import alignlib


for line in sys.stdin:
    print alignlib.translate(line.rstrip())

