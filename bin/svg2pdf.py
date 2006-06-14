#!/usr/bin/env python

import sys
from reportlab import svglib
from rasmus import util





if "-r" in sys.argv:
    remove = True
else:
    remove = False

for svgfile in sys.argv[1:]:
    if svgfile.startswith("-"):
        continue

    util.tic("converting '%s' to pdf" % svgfile)
    svglib.svg2pdf(svgfile)
    if remove:
        os.remove(svgfile)
    util.toc()

