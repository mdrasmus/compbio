import util, algorithms, genomeutil, blast




def plenFilter(plens, minalign):
    def filterFunc(line):
        query = blast.query(line)
        subject = blast.subject(line)
    
        # alignment percentage requirement
        if query in plens and \
           blast.queryLength(line) / float(plens[query]) < minalign:
            return False

        if subject in plens and \
           blast.subjectLength(line) / float(plens[subject]) < minalign:
            return False

        return True

    return filterFunc
    


def outgroupCutoff(files, orders, blastFilter=lambda x: True):
    util.tic("outgroup cutoffs")
    inBest = util.Dict(1, [0, None])    
    outBest = util.Dict(1, [0, None])

    # read in blast files
    for f, order in zip(files, orders):
        util.tic(f)
        for line in blast.BlastReader(f):
            # determine in and out genes
            if order:
                inGene = blast.query(line)
                outGene = blast.subject(line)
            else:
                inGene = blast.subject(line)
                outGene = blast.query(line)
            score = blast.bitscore(line)

            if not blastFilter(line):
                continue
                
            # find best unidirectional hits
            if score > inBest[inGene][0]:
                inBest[inGene] = [score, outGene]
            
            if score > outBest[outGene][0]:
                outBest[outGene] = [score, inGene]
        util.toc()
    util.toc()
    
    return inBest.data, outBest.data


def clusterIngroup(inBest, parts, sortHitsFile,
                   blastFilter=lambda x: True):
    util.tic("cluster")
    
    gene2comp = {}
    comp2cutoff = {}
    this = util.Closure(ncomps = 0,
                        printLimit = 0)

    def gene2cutoff(gene):
        return comp2cutoff[gene2comp[gene].root()]

    def getComp(gene):
        if gene not in gene2comp:
            comp = algorithms.UnionFind([gene])
            gene2comp[gene] = comp
            comp2cutoff[comp] = 0
            this.ncomps += 1
            this.printLimit += 1
            return comp
        else:
            return gene2comp[gene].root()
        
    
    # Make a set for each gene
    for gene in inBest:
        gene2comp[gene] = algorithms.UnionFind([gene])
        comp2cutoff[gene2comp[gene]] = inBest[gene][0]
    this.ncomps = len(gene2comp)
    this.printLimit = this.ncomps

    # merge based on initial partition
    for part in parts:
        comp1 = getComp(part[0])
        for gene2 in part[1:]:
            comp2 = getComp(gene2)
            cutoff = max(comp2cutoff[comp1.root()], comp2cutoff[comp2.root()])
            comp1.union(comp2)
            comp2cutoff[comp1.root()] = cutoff
            this.ncomps -= 1
    
    maxcomp = max(map(len, parts) + [1])
    nmatches = 0

    # read hits
    for line in blast.BlastReader(sortHitsFile):
        nmatches += 1
        gene1 = blast.query(line)
        gene2 = blast.subject(line)        
        score = blast.bitscore(line)

        if not blastFilter(line):
            continue
        
        comp1 = getComp(gene1)
        comp2 = getComp(gene2)
        
        if comp1 == comp2:
            continue

        cutoff = max(comp2cutoff[comp1], comp2cutoff[comp2])

        if score > cutoff:
            comp1.union(comp2)
            comp2cutoff[comp1.root()] = cutoff
            this.ncomps -= 1
            maxcomp = max(maxcomp, len(comp1.root().items))
            if this.ncomps < this.printLimit:
                print this.ncomps, maxcomp, nmatches, score, cutoff, gene1, gene2
                this.printLimit -= 1

    # find unique sets
    sets = {}
    for set in gene2comp.values():
        sets[set.root()] = 1

    # convert to component array
    comps = []
    for set in sets:
        comps.append(set.members())

    # find final cutoff scores for each component
    scores = map(lambda set: comp2cutoff[set], sets.keys())
    
    util.toc()

    return comps, scores
    

def getProteinLengths(fastaFiles):
    plens = {}

    tic("determine protein lengths")
    for fast in fastaFiles:
        tic(fast)
        seqs = fasta.readFasta(fast)
        for key, value in seqs.iteritems():
            plens[key] = len(value)
        toc()
    toc()
    
    return plens


def clusterGenes(crossBlast, orders, sortedBlast, plens, plenCutoff=.6):
    inBest, outBest = outgroupCutoff(crossBlast, orders, plenFilter(plens, plenCutoff))
    comps, scores = clusterIngroup(inBest, outBest, parts,
                                   sortedBlast,
                                   plenFilter(plens, plenCutoff))
    return comps, scores
