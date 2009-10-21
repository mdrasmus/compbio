#!/usr/bin/python
#


import sys, os

from rasmus import treelib, util
from rasmus.bio import alignlib, fasta, phylip, blast


options = [
    ["f:", "fasta=", "fasta", "<fasta sequences>",
        {"single": True}],
    ["a:", "align=", "align", "<fasta alignment>",
        {"single": True}],
    ["F:", "fastaext=", "fastaext", "<fasta extension>",
        {"single": True,
         "default": ".fasta"}],
    ["D:", "distext=", "distext", "<blast matrix extension>",
        {"single": True,
         "default": ".blast.dist"}],
    ["B:", "bitsiteext=", "bitsiteext", "<bits/site matrix extension>",
        {"single": True,
         "default": ".bitsites.mat"}],
    ["C:", "coverageext=", "coverageext", "<coverage matrix extension>",
        {"single": True,
         "default": ".coverage.mat"}],
    ]


conf = util.parseOptions(sys.argv, options)

      

# create blast db
os.system("formatdb -o -i %s" % conf["fasta"])



util.tic("perform blast")

# blast fasta with itself
hits = blast.iterBestHitPerTarget(
            blast.blastp(conf["fasta"], conf["fasta"], split=100))

seqs = fasta.read_fasta(conf["fasta"])

# store hits in a matrix ordered to match alignment
if "align" in conf:
    aln = fasta.read_fasta(conf["align"])
    #aln = alignlib.mapalign(aln, valfunc=lambda seq: filter(lambda x: x != '-', seq))
    keys = aln.keys()
else:
    keys = seqs.keys()

mat = util.make_matrix(len(keys), len(keys))
mat2 = util.make_matrix(len(keys), len(keys))
mat3 = util.make_matrix(len(keys), len(keys))
lookup = util.list2lookup(keys)
seqlengths = [len(x) for x in util.mget(seqs, keys)]


# save hits to file
blasthits = file(conf["fasta"] + ".blastp", "w")
for hit in hits:
    print >>blasthits, "\t".join(hit)
    
    i = lookup[hit[0]]
    j = lookup[hit[1]]
    val = blast.bitscore(hit)
    length = blast.alignLength(hit)
    mat[i][j] = val
    mat[j][i] = val
    mat2[i][j] = val / float(length)
    mat2[j][i] = val / float(length)
    mat3[i][j] = length / float(max(seqlengths[i], seqlengths[j]))
    mat3[j][i] = length / float(max(seqlengths[i], seqlengths[j]))
blasthits.close()

blastmat = util.replace_ext(conf["fasta"], conf["fastaext"], conf["distext"])
phylip.write_dist_matrix(mat, out=blastmat)

blastmat2 = util.replace_ext(conf["fasta"], conf["fastaext"], conf["bitsiteext"])
phylip.write_dist_matrix(mat2, out=blastmat2)

blastmat3 = util.replace_ext(conf["fasta"], conf["fastaext"], conf["coverageext"])
phylip.write_dist_matrix(mat3, out=blastmat3)

util.toc()
