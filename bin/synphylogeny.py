#!/usr/bin/python

import sys
from rasmus import synphyllib
from rasmus import fasta, util, genomeio
from rasmus import muscle, phylip, algorithms, clustalw

if len(sys.argv) < 4:
    print "usage: synphyl.py [-g <gene file>] [-d <desc>] [-p <part file>]"
    print "       [-s <syntenic components>] [-f <fasta>] [-o <output prefix>]"
    print "       [-x <part-part sim>] [-i <id>,<id>,...]"
    sys.exit(1)

param, rest = util.getopt(sys.argv[1:], "g:d:s:o:f:p:x:i:")

util.tic("read input")
parts = []
if "-g" in param:
    genes = []
    for f in param["-g"]:
        genes.extend(util.readStrings(f))
    parts.append(genes)

if "-p" in param:
    for f in param["-p"]:
        parts.extend(util.readDelim(f))


# read sequence
seqs = {}
for f in param["-f"]:
    seqs.update(fasta.readFasta(f, valuefunc=fasta.removestar))


# read synteny components
comps = []
for f in param["-s"]:
    comps.extend(util.readDelim(f))
complookup = {}
for i in xrange(len(comps)):
    for gene in comps[i]:
        complookup[gene] = i



# read descriptions
desc = {}
for f in param["-d"]:
    desc.update(genomeio.readGeneDesc(f))
util.toc()

# read part-part sim file
partsim = {}
if "-x" in param:
    for line in file(param["-x"][-1]):
        part1, part2, score = line.split()
        part1 = int(part1)
        part2 = int(part2)
        score = float(score)

        if part1 == part2:
            continue

        if part1 not in partsim or score > partsim[part1][1]:
            partsim[part1] = [part2, score]
        if part2 not in partsim or score > partsim[part2][1]:
            partsim[part2] = [part1, score]

if "-i" in param:
    indices = []
    for token in param["-i"]:
        indices.extend(map(int, token.split(",")))
else:
    indices = range(len(parts))    


outprefix_base = param["-o"][0]

def writeOrths(f, orths):
    for orth in orths:
        for gene in orth:
            print >>f, gene,
        print >>f

easyfile = file(outprefix_base + "easy.proper.orths", "w")
nosynfile = file(outprefix_base + "rest.nosyn.orths", "w")
properfile = file(outprefix_base + "all.proper.orths", "w")
extendedfile = file(outprefix_base + "all.extended.orths", "w")
ambiguousfile = file(outprefix_base + "all.ambiguous.orths", "w")


def processGenes(genes, outprefix):    
    # find subset of sequences
    seqs2 = util.getkeys(seqs, genes)    
                
    # determine how many different synteny components in sequence set
    # only build tree if synteny needs to be untangled or there are a lot of
    # genes
    nsyn = len(util.unique(filter(lambda x: x, 
                    map(lambda x: complookup.get(x, None), genes))))
    
    if True:  #nsyn > 1 or len(seqs2) > 6:
        # output input data
        util.writeVector(outprefix + ".genes", genes)
        fasta.writeFasta(file(outprefix + ".fasta", "w"), seqs2)
        
        # get outgroup
        #outgroup = parts[partsim[i][0]]
        #outgroup = filter(lambda x: x in seqs, outgroup)
        
        #util.writeVector(outprefix + ".outgroup", outgroup)
        
        # add sequence for outgroup
        #seqs2.update(util.getkeys(seqs, outgroup))
        
        
        # build initial tree
        util.tic("build tree")
        
        aln = muscle.muscle(seqs2)
        #aln = clustalw.clustalw(seqs2)
        
        fasta.writeFasta(file(outprefix + ".align", "w"), aln)
        
        tree = phylip.proml(aln, force=True)
        util.toc()
        
        if tree.root == None:
            util.log("could not create tree")
            file(outprefix + ".BADTREE", "w")
            return
        
        #tree.writeNewick(file(outprefix + ".full.tree", "w"))
        #tree, good = algorithms.removeOutgroup(tree, outgroup)
        #if not good:
        #    file(outprefix + ".UNROOTED", "w")
        #    util.log("WARNING: outgroup could not be separated from tree")
        

        # write stats/intermediate data on homology groups
        util.tic("write output")        
        synphyl.writeStats(file(outprefix + ".stats", "w"), tree, comps, desc)
        tree.writeNewick(file(outprefix + ".tree", "w"))
        util.toc()

        stats = synphyl.allGroups(tree, comps)
        writeOrths(file(outprefix + ".proper.orths", "w"), stats["proper"])
        writeOrths(file(outprefix + ".extended.orths", "w"), stats["extended"])
        writeOrths(file(outprefix + ".ambiguous.orths", "w"), stats["ambiguous"])
        
        writeOrths(properfile, stats["proper"])
        writeOrths(extendedfile, stats["extended"])
        writeOrths(ambiguousfile, stats["ambiguous"])
    else:
        if nsyn == 1:
            writeOrths(easyfile, [genes])
        else:
            writeOrths(nosynfile, [genes])
    


util.tic("SYNPHYL")

if "-p" in param:
    for i in indices:
        util.tic("process part %d of %d (size=%d)" % (i, len(parts), len(parts[i])))

        outprefix = outprefix_base + str(i)    
        genes = parts[i]
        processGenes(genes, outprefix)

        util.toc()
elif "-g" in param:
    processGenes(genes, outprefix_base)
util.toc()

easyfile.close()
nosynfile.close()
properfile.close()
extendedfile.close()
ambiguousfile.close()

