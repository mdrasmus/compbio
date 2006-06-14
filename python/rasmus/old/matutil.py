#!/usr/bin/python

import sys
import util
import matrix

filename = sys.argv[1]
filenameout = filename + ".mat"

#mat = matrix.AdjList(filename)

#mat.write_cluto(filenameout)

infile = file(filename)

(nrows, ncols, nnz) = map(int, infile.next().split())

prog = util.Progress(0, 0, nnz, .01)

rows = []
for i in range(nrows):
   rows.append([])

print "reading..."
for line in infile:
   prog.update()
   (r, c, v) = line.split()
   row = int(r)
   col = int(c)
   val = float(v)
   
   rows[row].append((col, val))

infile.close()

out = file(filenameout, "w")

out.write("%d %d %d\n" % (nrows, ncols, nnz))

prog = util.Progress(0, 0, nrows, .01)
for row in range(len(rows)):
   prog.update()
   for entry in rows[row]:
      out.write("%d %f " % (entry[0], entry[1]))
   out.write("\n")
