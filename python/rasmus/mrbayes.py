import fasta, util, treelib, exonutil, phylip
import os, StringIO



def mrbayes(aln, nexfilename = "", format="pep", options=None, 
            verbose=True, saveOutput=""):
    util.tic("mrbayes on %d of length %d" % (len(aln), len(aln.values()[0])))
    
    
    if nexfilename == "":
        cwd = phylip.createTempDir()
    else:
        cwd = None
    
    # setup options
    if nexfilename == "":
        nexfilename = "infile.nex"    
    if not options: 
        options = {}
    setDefaultOptions(options)
    
    burnin = .25 * options["ngen"] / options["samplefreq"]
    options["extra"] += "sumt;"
    
    # write input file
    writeNexus(file(nexfilename, "w"), aln.keys(), aln.values(), format, options)
    
    # exec mrbayes
    if verbose:
        os.system("echo exe %s | mb" % nexfilename)
    else:
        os.system("echo exe %s | mb >/dev/null 2>&1" % nexfilename)
    
    # read tree
    tree = readNexusConTree(file(nexfilename + ".con"))
    
    # clean up
    if cwd != None:
        if saveOutput != "":
            phylip.saveTempDir(cwd, saveOutput)
        else:
            phylip.cleanupTempDir(cwd)
    
    util.toc()
    
    return tree
    
    

def writeNexus(out, names, seqs, format="pep", options=None):
    # setup options
    if options is None:
        options = {}
    
    # determine sequence format
    if format == "pep":
        format = "protein"
    elif format == "dna":
        format = "dna"
    else:
        raise Exception("unknown sequence format")
    
    # write Nexus file header
    print >>out, """\
#NEXUS
begin data;
dimensions ntax=%d nchar=%d;
format datatype=%s interleave=yes gap=-;
matrix
""" % (len(seqs), len(seqs[0]), format)
    
    # write sequences
    for name, seq in zip(names, seqs):
        print >>out, "%s\t%s" % (name, seq)
    
    print >>out, """\
;
end;
"""
    # write options
    writeMrbayesOptions(out, options)




# todo: make slide an argument
def mrbayesIntrons(aln, genes, lookup, slide=1, filename = "", format="protein", options=None):
    if not options: options = {}
    if filename == "":
        filename = util.tempfile(".", "mrbayes-in", ".nex")
    util.tic("mrbayes on %d of length %d" % (len(aln), len(aln.values()[0])))
    
    setDefaultOptions(options)
    
    burnin = .25 * options["ngen"] / options["samplefreq"]
    options["extra"] += "sumt;"
    
    
    exonutil.writeNexusIntrons(file(filename, "w"), aln.keys(), aln.values(), 
        genes, lookup, slide=slide, options=options)
    os.system("echo exe %s | mb" % filename)
    
    tree = readNexusConTree(file(filename + ".con"))
    util.toc()
    
    return tree    

def writeNexusIntrons(out, names, seqs, intronseqs, format="protein", options=None):
    if not options: options = {}
    lseq = len(seqs[0])
    lintrons = len(intronseqs[0])

    print >>out, """\
#NEXUS
begin data;
dimensions ntax=%d nchar=%d;
format datatype=mixed(%s:1-%d,restriction:%d-%d) interleave=yes gap=-;
matrix
""" % (len(seqs), lseq+lintrons, format, lseq, lseq+1, lseq+lintrons)

    for name, seq in zip(names, seqs):
        print >>out, "%s\t%s" % (name, seq)
    
    for name, seq in zip(names, intronseqs):
        print >>out, "%s\t%s" % (name, seq)
    
    print >>out, """\
;
end;
"""

    writeMrbayesOptions(out, options)


def setDefaultOptions(options):
    defaults = {
        "aamodelpr": "fixed(blosum)",
        "ngen": 10000,
        "samplefreq": 10,
        "extra": "",
        "stoprule" : "yes"
    }
    
    # set defaults
    for key, val in defaults.items():    
        options.setdefault(key, val)  
    return options
    
    
def writeMrbayesOptions(out, options=None):
    if not options: options = {}
    setDefaultOptions(options)

    print >>out, """\
begin mrbayes;
    [ batch mode ]
    set autoclose=yes nowarn=yes;
    
    [ amino acid mutation prior ]
    prset aamodelpr = %(aamodelpr)s;
    
    [ begin MCMC run ]
    mcmc ngen=%(ngen)d samplefreq=%(samplefreq)d stoprule=%(stoprule)s;
    %(extra)s
end;
""" % options


def readMrbayesProb(infile):
    runs = []
    delim = "\t"
    data = None
    
    for line in infile:
        # hack to handel comments
        if line.startswith("["):
            data = []
            runs.append(data)
            
            # read header
            headers = infile.next().rstrip().split(delim)
        else:
            vals = line.rstrip().split(delim)
            
            dct = {}
            for i in xrange(len(vals)):
                dct[headers[i]] = float(vals[i])
            data.append(dct)
    
    return runs
        

def readNexusConTree(infile):
    count = 0
    for line in infile:
        if line.startswith("   tree con_50_majrule ="):
            count += 1
            
            # only read the second tree
            if count == 2:
                line = line.replace("   tree con_50_majrule =", "")
                tree = treelib.Tree()
                tree.readNewick(StringIO.StringIO(line))
                
                return tree
        
