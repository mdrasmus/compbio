#!/usr/bin/env python

from rasmus import util, fasta, genomeio, genomeutil, ensembl, exonutil
import sys




options = [
    ["g:", "genomes=", "genomes", "AUTO<genome1>,<genome2>,..."],
    ["n", "nexus", "nexus", "AUTO"]
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)

genomeNames = param["genomes"][-1].split(",")




# read peptides
peps = {}
for genome in genomeNames:
    util.tic("read %s peptides" % genome)
    peps.update(genomeio.readPeptides(genomeio.pepfile(genome), "trans"))
    util.toc()
lookup = util.revdict(peps)

# read genes
genomes = {}
genes = {}
for genome in genomeNames:
    util.tic("read %s genes" % genome)
    genomes[genome] = genomeutil.Genome(genome)
    genomeio.readEnsmartGenes(genomes[genome])
    genes.update(genomes[genome].genes)
    util.toc()


if "nexus" in param:
    for fn in rest:
        util.log("processing '%s'" % fn)
        out = file(fn.replace(".align", ".nexus"), "w")
        names, seqs = fasta.readFastaOrdered(fn)
        exonutil.writeNexusIntrons(out, names, seqs, genes, lookup, slide=1)
        out.close()


else:
    # usual intron fasta vis


    def makeGene2trans(names, seqs, lookup):
        gene2trans = {}
        for name, seq in zip(names, seqs):
            try:
                trans = lookup[seq.replace("-", "")]
            except:
                trans = lookup[seq.replace("-", "")[:-1]]
            gene2trans[name] = trans
        return gene2trans


    def readExons(genes):
        genomes = util.Dict(1, [])
        for x in genes:
            genomes[x.chrom.genome.name].append(x.chrom.name)
        for x in genomes:
            genomes[x] = util.unique(genomes[x])

        genomes2 = {}
        for genome in genomes:
            util.tic("read %s exons" % genome)
            genomes2[genome] = genomeutil.Genome(genome)

            for chrom in genomes[genome]:
                genomeio.readEnsmartExons(genomes2[genome], 
                                          genomeio.structfile(genome, chrom))
            util.toc()
        return genomes2


    for fn in rest:
        util.log("processing '%s'" % fn)
        names, seqs = fasta.readFastaOrdered(fn)
        gene2trans = makeGene2trans(names, seqs, lookup)
        genes2 = util.subdict(genes, names).values()
        genomes2 = readExons(genes2)

        # get all transcripts
        transcripts = {}
        for genome in genomes2.values():
            for gene in genome.genes.values():
                transcripts.update(gene.trans)

        out = file(fn.replace(".align", ".exon.align"), "w")

        for name, seq in zip(names, seqs):
            trans = transcripts[gene2trans[name]]

            elens = []
            for exon in trans.exons:
                if exon.codeStart != None:
                    elens.append((exon.codeEnd - exon.codeStart + 1))

            util.log(name, elens)

            print >>out, ">%s" % name
            j = 0
            k = 0
            for i in xrange(len(seq)):
                if seq[i] != "-":
                    j += 3
                if k < len(elens) and j >= elens[k]:
                    if k < len(elens) - 1:
                        ilen = min(abs(trans.exons[k+1].start - trans.exons[k].end),
                                   abs(trans.exons[k+1].end - trans.exons[k].start))
                        if ilen < 10:
                            out.write(".")
                        else:
                            out.write("#")
                    else:
                        out.write(seq[i])
                    j -= elens[k]
                    k += 1
                else:
                    out.write(seq[i])
            print >>out

        out.close()
