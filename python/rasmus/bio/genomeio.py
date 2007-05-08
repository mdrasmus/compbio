"""


Tue Feb 20 23:16:41 EST 2007
Matt Rasmussen
This all old code that I should throw away some day



"""


# python libs
import os

# rasmus libs
from rasmus.bio.genomeutil import *
from rasmus import util
from rasmus import env




# ensembl field names
FIELD_GENE_NAME       = 'Ensembl Gene ID'
FIELD_GENE_START      = 'Start Position (bp)'
FIELD_GENE_END        = 'End Position (bp)'
FIELD_STRAND          = 'Strand'
FIELD_CHROM_NAME      = 'Chromosome Name'
FIELD_TRANS_NAME      = 'Ensembl Transcript ID'
FIELD_EXON_NAME       = 'Ensembl Exon ID (versioned)'
FIELD_EXON_START      = 'Exon Start (Chr bp)'
FIELD_EXON_END        = 'Exon End (Chr bp)'
FIELD_EXON_CODE_START = 'Coding Start (Chr bp)'
FIELD_EXON_CODE_END   = 'Coding End (Chr bp)'


# default filename conventions
def codingfile(genome):
    return "%s.nt.fasta" % genome

def pepfile(genome):
    return "%s.fasta" % genome

def coordfile(genome):
    return "%s.coord" % genome

def genefile(genome):
    return "%s.gff" % genome


# REMOVE ME
env.genomepath = ""




def readGenomes(matching, genomes, gene2species=gene2species, 
                paths = env.datapaths):
    util.tic("reading genomes")
    for genome in genomes:
        util.tic("reading %s" % genome)
        
        try:
            fn = env.findFile(genefile(genome), paths)
        except env.PathError:    
            fn = env.findFile(coordfile(genome), paths)
        
        matching.readGenomes(fn, gene2species)
        util.log("read %d genes" % len(matching.genomes[genome].genes))
        util.toc()
    matching.autoconf(genomes)
    util.toc()


def readGenomeExons(genomes, datapath = env.genomepath):
    for genome in genomes:
        util.tic("read exons for " + genome)
        readEnsmartExons(genomes[genome], structfile(genome, datapath=datapath))
        util.toc()


def readGenomeGo(genomes, datapath = env.genomepath):
    for genome in genomes:
        util.tic("read go for "+genome)
        readGo(genomes[genome], gofile(genome, datapath))
        util.toc()

def readGenomeSeq(genomes, datapath = env.genomepath, 
                           protein = False,
                           dna = True):
    for genome in genomes:
        genomes[genome].seqdb = SequenceDB()
        if dna:
            genomes[genome].seqdb.addPath(datapath + genome + "/")
        if protein:
            genomes[genome].seqdb.addProteinFasta(
                datapath + genome + "/" + genome +".fasta")


def readGenes(genomes, key="name", datapath = env.genomepath):
    tables = []
    headers = []
    genes = {}
    for genome in genomes:
        if os.path.isfile(genefile(genome, datapath)):
            table, header = util.readTable(genefile(genome, datapath))
            genes.update(util.lookupTable(table, key))
            tables.append(table)
            headers.append(header)
    return genes


def readBlockDimensions(matching, filename, matchAssign = False):
    infile = util.DelimReader(filename, "\t")
    
    blockstart = len(matching.blocks)
    
    for line in infile:
        if line[0] not in matching.genomes or \
           line[4] not in matching.genomes:
            continue
        
        block = SyntenyBlock()
        block.blockid = len(matching.blocks)
        block.genome1 = matching.genomes[line[0]]
        if line[1] not in block.genome1.chroms:
            continue
        block.chrom1  = block.genome1.chroms[line[1]]
        block.start1  = int(line[2])
        block.end1    = int(line[3])
        block.genome2 = matching.genomes[line[4]]
        if line[5] not in block.genome2.chroms:
            continue
        block.chrom2  = block.genome2.chroms[line[5]]
        block.start2  = int(line[6])
        block.end2    = int(line[7])
        block.direction = (line[8] == "1")
        matching.blocks.append(block)
    
    if matchAssign:
        def between(x, a, b):
            return a <= x and x <= b
        
        # hash blocks
        blockhash = util.Dict(2, [])
        
        for block in matching.blocks[blockstart:]:
            blockhash[block.chrom1][block.chrom2].append(block)
        
        #prog = util.ProgressBar(len(matching.matches))
        for match in matching.matches:
            if match.block != None:
                continue
            
            #prog.update()
            gene1, gene2 = match.genes
            
            for block in blockhash[gene1.chrom][gene2.chrom]:
                if (block.start1 <= gene1.start <= block.end1) and \
                   (block.start2 <= gene2.start <= block.end2):
                    block.add(match, init = False)
                    match.block = block
                    continue
            
            for block in blockhash[gene2.chrom][gene1.chrom]:
                if (block.start1 <= gene2.start <= block.end1) and \
                   (block.start2 <= gene1.start <= block.end2):
                    block.add(match, init = False)
                    match.block = block
                    continue

                    


def readGeneDesc(filename):
    infile = file(filename)
    desc = {}

    for line in infile:
        line = line.rstrip()
        tokens = line.split("\t")
        if len(tokens) > 1:
            desc[tokens[0]] = tokens[1]
    
    return desc



def readSomeMatches(matching, infilename, genome1, genome2, pred):
    infile = file(infilename)
    
    for line in infile:
        (gene1, gene2, score) = line.split()[:3]
        
        if pred(gene1, gene2, score):
            matching.addMatch(matching.genomes[genome1].genes[gene1],
                              matching.genomes[genome2].genes[gene2],
                              float(score))

def readSimpleMatches(matching, filename):
    genes = matching.getGenes()

    for line in file(filename):
        words = line.rstrip().split()
        for i in range(len(words)):
            for j in range(i+1, len(words)):
                if words[i] in genes and words[j] in genes:
                    matching.addMatch(genes[words[i]], genes[words[j]], 0)


def readEnsmartGenes(genome, filename  = None):
    if filename == None:
        filename = structfile(genome.name)
    infile = file(filename)
    delim = "\t"
    
    config = {"source": "ensembl"} #readConfig(genome.name)
    
    # determine keys of ensmart
    keys = infile.next().rstrip().split(delim)
    lookup = util.list2lookup(keys)
    fields = {}
    
    # check which fields are in use
    assert (FIELD_GENE_NAME in keys and \
            FIELD_GENE_START in keys and \
            FIELD_GENE_END in keys and \
            FIELD_STRAND in keys and \
            FIELD_CHROM_NAME in keys), "required fields missing"

    
    # field ids
    geneId = lookup[FIELD_GENE_NAME]
    geneStartId = lookup[FIELD_GENE_START]
    geneEndId = lookup[FIELD_GENE_END]
    strandId = lookup[FIELD_STRAND]
    chromId = lookup[FIELD_CHROM_NAME]
    
    
    # read lines
    for line in infile:
        tokens = line.rstrip().split(delim)

        if config["source"] == "ensembl":
            geneName = tokens[geneId].split(".")[0]
        else:
            geneName = tokens[geneId]
            
        try:
            gene = genome.genes[geneName]
        except KeyError:
            gene = Gene()
            gene.name = geneName
            genome.genes[geneName] = gene
            
            gene.start = int(tokens[geneStartId])
            gene.end = int(tokens[geneEndId])
            gene.direction = int(tokens[strandId])

            chromName = tokens[chromId]
            try:
                chrom = genome.chroms[chromName]
            except KeyError:
                chrom = Chromosome(chromName)
                chrom.genome = genome
                genome.chroms[chromName] = chrom

            gene.chrom = chrom

def genekey(key):
    return key.split("|")[0]

def transkey(key):
    return key.split("|")[1]

def pepkey(key):
    return key.split("|")[2]

def readPeptides(filename, key = "gene"):
    keys = {
        "gene": genekey,
        "trans": transkey,
        "pep": pepkey
    }
    
    if key not in keys:
        raise "unknown key type"
    else:
        return fasta.readFasta(filename, keyfunc = keys[key])

    

def readEnsmartExons(genome, filename = None, genes=None):
    if filename == None:
        filename = structfile(genome.name)
    infile = file(filename)
    delim = "\t"
    
    # determine keys of ensmart
    keys = infile.next().rstrip().split(delim)
    lookup = util.list2lookup(keys)
    fields = {}
    
    config = readConfig(genome.name)
    
    # check which fields are in use
    assert (FIELD_GENE_NAME in keys and \
            FIELD_GENE_START in keys and \
            FIELD_GENE_END in keys and \
            FIELD_STRAND in keys and \
            FIELD_CHROM_NAME in keys and \
            FIELD_EXON_NAME in keys and \
            FIELD_EXON_START in keys and \
            FIELD_EXON_END in keys and \
            FIELD_EXON_CODE_START in keys and \
            FIELD_EXON_CODE_END in keys), "required fields missing"
    
    # field ids
    geneId = lookup[FIELD_GENE_NAME]
    geneStartId = lookup[FIELD_GENE_START]
    geneEndId = lookup[FIELD_GENE_END]
    strandId = lookup[FIELD_STRAND]
    chromId = lookup[FIELD_CHROM_NAME]
    transId = lookup[FIELD_TRANS_NAME]
    exonId = lookup[FIELD_EXON_NAME]
    exonStartId = lookup[FIELD_EXON_START]
    exonEndId = lookup[FIELD_EXON_END]
    exonCodeStartId = lookup[FIELD_EXON_CODE_START]
    exonCodeEndId = lookup[FIELD_EXON_CODE_END]
    
    
    # read lines
    for line in infile:
        tokens = line.rstrip().split(delim)
        
        if config["source"] == "ensembl":
            geneName = tokens[geneId].split(".")[0]
        else:
            geneName = tokens[geneId]
        
        
        if genes != None:
            if geneName not in genes:
                continue
        
        try:
            gene = genome.genes[geneName]
        except KeyError:
            gene = Gene()
            gene.name = geneName
            genome.genes[geneName] = gene

            gene.start = int(tokens[geneStartId])
            gene.end = int(tokens[geneEndId])
            gene.direction = int(tokens[strandId])

            chromName = tokens[chromId]
            try:
                chrom = genome.chroms[chromName]
            except KeyError:
                chrom = Chromosome(chromName)
                chrom.genome = genome
                genome.chroms[chromName] = chrom

            gene.chrom = chrom

        # read transcript
        transName = tokens[transId]
        try:
            trans = gene.trans[transName]
        except KeyError:
            trans = Transcript(gene)
            trans.name = transName
            gene.trans[transName] = trans
        
        # read exon
        if exonCodeStartId < len(tokens):
            codeStart = int(tokens[exonCodeStartId])
            codeEnd   = int(tokens[exonCodeEndId])
        else:
            codeStart = None
            codeEnd   = None
        
        trans.exons.append(Exon(trans, 
                                int(tokens[exonEndId]),
                                int(tokens[exonStartId]),
                                int(tokens[strandId]),
                                codeStart,
                                codeEnd))

        trans.exons[-1].name = tokens[exonId]

    sortExons(genome)




def sortExons(genome):
    for gene in genome.genes.values():
        for trans in gene.trans.values():
            if gene.direction == 1:
                trans.exons.sort(lambda a,b: a.start - b.start)
            else:
                trans.exons.sort(lambda a,b: b.start - a.start)


"""
def writeGenome(filename, genome):
    outfile = file(filename, "w")
    
    for gene in genome.genes.values():
        d = ""
        if gene.direction == 1:
            d = "+"
        else:
            d = "-"
        arr = [gene.name, gene.chrom.name, str(gene.start), str(gene.end), d]
        print >>outfile, "\t".join(arr)
"""



def writeMatches(filename, matches, scores=False):
    out = file(filename, "w")
    if scores:
        for match in matches:
            print >>out, match.genes[0].name, match.genes[1].name, match.score
    else:    
        for match in matches:
            print >>out, match.genes[0].name, match.genes[1].name
    out.close()

def writeSynteny(filename, blocks):
    out = file(filename, "w")
    for block in blocks:
        out.write("\t".join(map(str, [block.genome1.name, block.chrom1.name, \
              block.start1, block.end1, \
              block.genome2.name, block.chrom2.name, \
              block.start2, block.end2, \
              block.getDirection()])))
        
        #for match in block.matches:
        #    out.write("\t")
        #    out.write(match.genes[0].name)
        #    out.write("\t")
        #    out.write(match.genes[1].name)
        out.write("\n")
    out.close()





# OLD FUNCTIONS
    
"""
def readDesc(genome, filename):
    infile = file(filename)

    for line in infile:
        line = line.rstrip()
        tokens = line.split("\t")
        if genome.genes.has_key(tokens[0]) and len(tokens) > 1:
            genome.genes[tokens[0]].description = tokens[1]


def readExons(genome, filename):
    infile = file(filename)

    for line in infile:
        tokens = line.split()
        if not(tokens[0] in genome.genes):
            continue
        gene = genome.genes[tokens[0]]

        # skip gene if it already has exons (ignore alt. splicings)
        if len(gene.trans) == 0:    
            gene.trans["1"] = Transcript(gene, "1")
            trans = gene.trans["1"]
            for i in range(1, len(tokens), 3):
                trans.exons.append(Exon(trans, 
                    int(tokens[i]), int(tokens[i+1]), int(tokens[i+2])))


def readGo(genome, filename):
    infile = file(filename)

    for line in infile:
        tokens = line.split()
        if not(tokens[0] in genome.genes):
            continue
        gene = genome.genes[tokens[0]]
        gene.gos = []
        for i in range(1, len(tokens)):
            gene.gos.append(tokens[i])


# default filename conventions
def structfile(genome, chrom = None, datapath = env.genomepath):
    if chrom == None:
        return datapath + genome + "/" + genome + ".struct"
    else:
        return datapath + genome + "/" + genome + "_" + chrom + ".struct"

def genefile(genome, datapath = env.genomepath):
    return datapath + genome + "/" + genome + "_genes.tab"

def pepfile(genome, datapath = env.genomepath):
    return datapath + genome + "/" + genome + ".pep.fasta"

def simplepepfile(genome, datapath = env.genomepath):
    return datapath + genome + "/" + genome + ".fasta"

def genomefile(genome, datapath = env.genomepath):
    return datapath + genome + "/" + genome + ".genome"

def configfile(genome, datapath = env.genomepath):
    return datapath + genome + "/config"

def gofile(genome, datapath = env.genomepath):
    return datapath + genome + "/" + genome + ".go"

def descfile(genome, datapath = env.genomepath):
    return datapath + genome + "/" + genome + ".desc"

    

def readConfig(genome):
    config = {
        "source": "standard"
    }
    
    if os.path.isfile(configfile(genome)):
        for line in file(configfile(genome)):
            tokens = line.split("\t")
            config[tokens[0]] = tokens[1]
        
    return config



# automatic reading
def readGenomes2(matching, genomes, datapath = env.genomepath):
    order = []
    
    util.tic("reading genomes")
    for genome in genomes:
        util.tic("reading %s" % genome)
        fn = structfile(genome, datapath=datapath)
        matching.genomes[genome] = Genome(genome)
        readEnsmartGenes(matching.genomes[genome], fn)
        order.append(matching.genomes[genome])
        util.toc()
    matching.autoconf(order)
    util.toc()
"""

# END OLF FUNCTIONS
