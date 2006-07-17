#!/usr/bin/env python

import sys

from rasmus import env


env.addEnvPaths("DATAPATH")

for filename in sys.argv[1:]:
    print env.findFile(filename)
