import os, StringIO
from itertools import chain

from rasmus import util
from rasmus import treelib
from . import fasta
from . import phylip




def mrbayes(aln, nexfilename = "", seqtype="pep", options=None, 
            usertree=None, bootiter=0,
            verbose=True, saveOutput=""):
    util.tic("mrbayes on %d of length %d" % (len(aln), len(aln.values()[0])))
    
    
    if nexfilename == "":
        cwd = phylip.create_temp_dir()
    else:
        cwd = None
    
    # setup options
    if nexfilename == "":
        nexfilename = "infile.nex"    
    if not options: 
        options = {}
    setDefaultOptions(options)
    
    options["burninfrac"] = .25
    options["relburnin"] = "yes"
    
    # force best binary tree (if possible)
    options["extra"] += "sumt contype=allcompat;"
    
    # get gene names
    names = []
    namemap = {}
    
    for key in aln.keys():
        if "+" in key:
            key2 = key.replace("+", "_")
            names.append(key2)
            namemap[key2] = key
        else:
            names.append(key)
    
    
    # write input file
    out = file(nexfilename, "w")
    writeNexus(out, names, aln.values(), seqtype, options)

    # write options
    writeMrbayesOptions(out, options, seqtype=seqtype)
    out.close()
    
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
            phylip.save_temp_dir(cwd, saveOutput)
        else:
            phylip.cleanup_temp_dir(cwd)
    
    util.toc()
    
    
    for tmpname, origname in namemap.iteritems():
        tree.rename(tmpname, origname)
    
    return tree
    
    

def write_nexus(out, names, seqs, format="pep", options=None, seqwidth=1000):
    # setup options
    if options is None:
        options = {}
    
    # determine sequence format
    if format == "pep":
        format2 = "protein"
    elif format == "dna":
        format2 = "dna"
    else:
        raise Exception("unknown sequence format")
    
    # write Nexus file header
    print >>out, """\
#NEXUS
begin data;
dimensions ntax=%d nchar=%d;
format datatype=%s interleave=yes gap=-;
matrix
""" % (len(seqs), len(seqs[0]), format2)
    
    totalsize = len(seqs[0])
    size = 0
    
    # write sequences
    while size < totalsize:
        for name, seq in zip(names, seqs):
            print >>out, "%s  %s" % (name, seq[size:size+seqwidth])
        size += seqwidth
    
    print >>out, """\
;
end;
"""
writeNexus = write_nexus    



'''
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
'''

def setDefaultOptions(options):
    defaults = {
        "aamodelpr": "fixed(blosum)",
        "ngen": 10000,
        "samplefreq": 10,
        "extra": "",
        "stoprule" : "yes",
        "burninfrac": .25,
        "relburnin": "yes"
    }
    
    # set defaults
    for key, val in defaults.items():    
        options.setdefault(key, val)  
    return options
    
    
def writeMrbayesOptions(out, options=None, seqtype="pep"):
    if not options: options = {}
    setDefaultOptions(options)

    print >>out, """\
begin mrbayes;
    [ batch mode ]
    set autoclose=yes nowarn=yes;
    """
    
    if seqtype == "pep":
        print >>out, """\
        [ amino acid mutation prior ]
        prset aamodelpr = %(aamodelpr)s;
        """ % options
    elif seqtype == "dna":
        print >>out, """\
        lset Nst = 2;
        """
    else:
        raise Exception("unknown seqtype '%s'" % seqtype)
    
    print >>out, """\
    [ begin MCMC run ]
    mcmc ngen=%(ngen)d samplefreq=%(samplefreq)d stoprule=%(stoprule)s
         relburnin=%(relburnin)s burninfrac=%(burninfrac)f;
    %(extra)s
end;
""" % options


def readMrbayesProb(infile):
    runs = []
    delim = "\t"
    data = None
    
    for line in infile:
        # hack to handle comments
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
        if line.startswith("   tree con_all_compat ="):
            count += 1
            
            # only read the second tree
            if count == 1:
                line = line.replace("   tree con_all_compat =", "")
                tree = treelib.Tree()
                tree.read_newick(StringIO.StringIO(line))
                
                return tree
    raise Exception("No tree found in output file")

 

def read_nexus_trees(infile):
    """Iterate over trees from trees run file"""
    
    infile = iter(infile)

    # skip until translate
    for line in infile:
        if "translate" in line:
            break

    # read translate
    names = {}
    for line in infile:
        if "tree rep" in line:
            break
        num, name = line.strip().replace(",", "").replace(";", "").split()
        names[num] = name

    # read trees
    for line in chain([line], infile):
        if "tree rep" not in line:
            break

        tree = treelib.parse_newick(line.split("=")[1])
        for oldname in tree.leaf_names():
            tree.rename(oldname, names[oldname])
        yield tree
