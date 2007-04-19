#!/usr/bin/env python


import sys, os


# get data environments
dataenv = os.getenv("DATAENV")
if dataenv == None:
    dataenv = []
else:
    dataenv = dataenv.split(":")


quiet = False

# process args
for arg in sys.argv[1:]:
    if arg == "-q":
        quiet = True
    elif arg == "-c":
        dataenv = []
    elif arg.startswith("-"):
        name = arg[1:]
        if name in dataenv:
            dataenv.remove(name)
    elif arg.startswith("+"):
        name = arg[1:]
        if name not in dataenv:
            dataenv.append(name)
    else:
        name = arg
        if name not in dataenv:
            dataenv.append(name)


filepath = os.getenv("FILEPATH")


# build datapaths
datapath = [filepath+"/data"]

if "fungi" in dataenv:
    datapath.append("%s/data/fungi" % filepath)
    datapath.append("%s/data/fungi/gene_coord" % filepath)
    datapath.append("%s/data/fungi/pep" % filepath)
    datapath.append("%s/data/fungi/nt" % filepath)
    datapath.append("%s/data/fungi/blast_pep" % filepath)

if "fungi-2" in dataenv:
    datapath.append("%s/data/fungi-2" % filepath)
    datapath.append("%s/data/fungi-2/gene_coord" % filepath)
    datapath.append("%s/data/fungi-2/gff" % filepath)    
    datapath.append("%s/data/fungi-2/pep" % filepath)
    datapath.append("%s/data/fungi-2/nt" % filepath)
    datapath.append("%s/data/fungi-2/blast_pep" % filepath)


if "fly_flybase" in dataenv:
    datapath.append("%s/data/flies/flybase" % filepath)
    datapath.append("%s/data/flies/flybase/gene_coord" % filepath)
    datapath.append("%s/data/flies/flybase/pep" % filepath)
    datapath.append("%s/data/flies/flybase/nt" % filepath)
    datapath.append("%s/data/flies/flybase/blast_pep" % filepath)


if "fly_flybase_single" in dataenv:
    datapath.append("%s/data/flies/flybase" % filepath)
    datapath.append("%s/data/flies/flybase/gene_coord_single" % filepath)
    datapath.append("%s/data/flies/flybase/pep" % filepath)
    datapath.append("%s/data/flies/flybase/nt" % filepath)
    datapath.append("%s/data/flies/flybase/blast_pep" % filepath)


if "fly_aaa" in dataenv:
    datapath.append("%s/data/flies/aaa" % filepath)
    datapath.append("%s/data/flies/aaa/gene_coord" % filepath)
    datapath.append("%s/data/flies/aaa/pep" % filepath)
    datapath.append("%s/data/flies/aaa/nt" % filepath)
    datapath.append("%s/data/flies/aaa/blast_pep" % filepath)



# export dataenv
print "export DATAENV=%s" % ":".join(dataenv)
print

# export datapath
print "export DATAPATH=%s" % ":".join(datapath)


if not quiet:
    print >>sys.stderr
    print >>sys.stderr, "Data Environments"
    print >>sys.stderr, "  " + " ".join(dataenv)
    print >>sys.stderr 
    
    print >>sys.stderr, "Data Paths"
    for path in datapath:
        print >>sys.stderr, "  " + path
    
