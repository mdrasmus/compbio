#!/usr/bin/env python

import sys

url = sys.argv[1]
print url.replace("&", "&amp;")
