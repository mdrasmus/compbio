import util, genomeutil, graph, genomeio, fasta, graph, ensembl, mrbayes


def findAlignedExons(aln, genomes):
    unionTrans = []
    exonData = {}
    exons = {}
    
    # find exons, their rows and aligned coords
    for i in xrange(len(aln)):
        row = aln[i]
        unionTrans.append(findExons(genomes[row.genome].chroms[row.chrom].genes,
                                    row.start, row.end))
        
        alignLookup = genomeutil.mkAlignLookup(row.seq)
        
        for name,exon in unionTrans[-1].iteritems():
            start, end = exonAlignCoords(exon, row, alignLookup)
            exonData[exon] = {"row":i, "start": start, "end": end}
            exons[exon] = 1
    
    # find connected components of overlaped exons
    mat = util.Dict(2)
    for exon1 in exons:
        for exon2 in exons:
            if util.overlap(exonData[exon1]["start"], 
                            exonData[exon1]["end"],
                            exonData[exon2]["start"], 
                            exonData[exon2]["end"]):
                mat[exon1][exon2] = 1
                mat[exon2][exon1] = 1
    
    def overlapping(vertex):
        return mat[vertex].keys()
    comps = graph.connectedComponents(exons.keys(), overlapping)
        
    return comps, exonData



    
     


def printStats(stats):    
    print
    print "total     %d" % stats["total"]
    print "conserved %d   [perfect %d (%f),  movement %d (%f)]" % \
        (stats["conserved"], stats["perfect"], stats["perfect"]/float(stats["conserved"]), 
         stats["movement"], stats["movement"]/float(stats["conserved"]))
    print "gain      %d" % stats["gain"]
    print "loss      %d" % stats["loss"]
    print 
    print "alternatively spliced %d" % stats["alt"]

    return stats




def isPerfectConserved(exons, exonData):
    starts = {}
    ends = {}
    
    for exon in exons:
        starts[exonData[exon]["start"]] = 1
        ends[exonData[exon]["end"]] = 1
    
    return len(starts) == 1 and len(ends) == 1

def isConserved(exons, exonData, nrows):
    rows = getRowComposition(exons, exonData)
    return len(rows) == nrows

def isWellAligned(exons, exonData, maxslip = 10):
    starts = [exonData[x]["start"] for x in exons]
    ends   = [exonData[x]["end"] for x in exons]
    
    return max(starts) - min(starts) < 10 and \
           max(ends) - min(ends) < 10
        

def getExonRows(exons, exonData):
    dct = util.Dict(1, [])
    for exon in exons:
        dct[exonData[exon]["row"]].append(exon)
    
    top = max(dct.keys()) + 1
    rows = []
    for i in xrange(top): rows.append([])
    for i in dct:
        rows[i] = dct[i]
    return rows

def getRowComposition(exons, exonData):
    rows = getExonRows(exons, exonData)
    return [len(x) for x in rows]


def exonAlignCoords(exon, row, alignLookup):
    a = genomeutil.global2align(exon.start, row.start, row.end, 
                                row.strand, alignLookup)
    b = genomeutil.global2align(exon.end, row.start, row.end, 
                                row.strand, alignLookup)
    
    if a < b:
        return (a, b)
    else:
        return (b, a)


def findExons(genes, start, end):
    exons = {}
    for gene in genes:
        for trans in gene.trans.values():
            for exon in trans.exons:
                if util.overlap(start, end, exon.start, exon.end):
                    exons[exon.name] = exon
    return exons





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

def readPepLookup(genomenames):
    peps = {}
    for genome in genomenames:
        util.tic("read %s peptides" % genome)
        peps.update(genomeio.readPeptides(genomeio.pepfile(genome), "trans"))
        util.toc()
    lookup = util.revdict(peps)
    return lookup


def getRelatedTranscripts(names, seqs, genes, lookup):
    gene2trans = makeGene2trans(names, seqs, lookup)
    genes2 = util.subdict(genes, names).values()
    genomes2 = readExons(genes2)
    
    # get all transcripts
    trans = genomeutil.getTranscripts(genomes2)
    
    trans2 = {}
    for gene in names:
        trans2[gene] = trans[gene2trans[gene]]
    
    return trans2


def readRelatedTranscripts(names, seqs, genomenames):
    m = genomeutil.Matching()
    genomeio.readGenomes(m, genomenames)
    genes = m.getGenes()
    lookup = readPepLookup(genomenames)
    trans = getRelatedTranscripts(names, seqs, genes, lookup)
    
    return genes, lookup, trans


def getIntronPositions(aln, transcripts):   
    ipos = {}
    
    # loop through old transcript
    for name, seq in aln.items():
        trans = transcripts[name]
        
        ipos[name] = []
        tot = 0
        for exon in trans.exons:
            if exon.codeStart != None:
                ipos[name].append(tot + exon.codeEnd - exon.codeStart + 1)
                tot = ipos[name][-1]
                
        # remove last boundary (it's not an intron)
        ipos[name].pop()
        
    return ipos


def insertIntrons(names, seqs, transcripts):
    aln = fasta.array2dict(names, seqs)
    ipos = getIntronPositions(aln, transcripts)

    names2 = []
    seqs2 = []
    
    # loop through old transcript
    for name, seq in zip(names, seqs):
        trans = transcripts[name]
        pos = ipos[name]
        
        names2.append(name)
        seq2 = ""
        j = 0
        k = 0
        for i in xrange(len(seq)):
            if seq[i] != "-":
                j += 3
            if k < len(pos) and j >= pos[k]:
                if k < len(pos):
                    ilen = min(abs(trans.exons[k+1].start - trans.exons[k].end),
                               abs(trans.exons[k+1].end - trans.exons[k].start))
                    if ilen < 10:
                        seq2 += "."
                    else:
                        seq2 += "#"
                else:
                    seq2 += seq[i]
                k += 1
            else:
                seq2 += seq[i]
        seqs2.append(seq2)

    return names2, seqs2        


def intronAlignment(aln, transcripts, slide=2):
    ipos = getIntronPositions(aln, transcripts)
    introns = []
    lookups = {}
    
    for name in ipos:
        for pos in ipos[name]:
            introns.append([name, pos])
        lookups[name] = genomeutil.mkAlignLookup(aln[name])
    
    # create graph
    mat = util.Dict(2, 0)
    for i in xrange(len(introns)):
        for j in xrange(i+1, len(introns)):
            name1, name2 = introns[i][0], introns[j][0]
        
            # cannot connect introns from same sequence
            if name1 == name2:
                continue
            
            p1 = lookups[name1][introns[i][1] / 3]
            p2 = lookups[name2][introns[j][1] / 3]
            
            if p1 > p2:
                p1, p2 = p2, p1
            seq1 = aln[name1][p1:p2].replace("-", "")
            seq2 = aln[name2][p1:p2].replace("-", "")
            
            if min(len(seq1), len(seq2)) <= slide:
                mat[i][j] = 1
                mat[j][i] = 1
    
    # find components
    def getNeighbor(v):
        return mat[v].keys()
    comps = graph.connectedComponents(range(len(introns)), getNeighbor)
    
    # sort components by their location in alignment (left to right)
    comps.sort(lambda a, b: cmp(introns[a[0]][1], introns[b[0]][1]))
    
    # now create intron alignment
    aln2 = {}
    for name in aln:
        aln2[name] = ""
    
    for comp in comps:
        names = map(lambda x: introns[x][0], comp)
        for name in aln:
            if name in names:
                aln2[name] += "1"
            else:
                aln2[name] += "0"
    
    return aln2 #, comps, introns, mat
    

def writeNexusIntrons(out, names, seqs, genes, lookup, slide=2, options={}):
    genomenames = util.unique(map(genomeutil.gene2species, names))
    trans = getRelatedTranscripts(names, seqs, genes, lookup)
    aln = fasta.array2dict(names, seqs)
    ialn = intronAlignment(aln, trans, slide)
    
    iseqs = util.sublist(ialn, names)
    
    if len(iseqs[0]) > 0:    
        mrbayes.writeNexusIntrons(out, names, seqs, iseqs, options=options)
    else:
        mrbayes.writeNexus(out, names, seqs, options=options)



