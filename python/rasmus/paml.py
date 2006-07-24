
# python libs
import os

# rasmus libs
import phylip
import util



def dndsMatrix(seqs, saveOutput="", verbose=False):
    
    phylip.validateSeq(seqs)
    cwd = phylip.createTempDir()

    util.tic("yn00 on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = phylip.fasta2phylip(file("seqfile.phylip", "w"), seqs)
    util.writeVector(file("labels", "w"), labels)    
    
    # create control file
    out = file("yn00.ctl", "w")
    print >>out, "seqfile = seqfile.phylip"
    print >>out, "outfile = outfile"
    out.close()
    
    # run phylip
    if verbose:
        os.system("yn00 yn00.ctl")
    else:
        os.system("yn00 yn00.ctl > /dev/null")
    
    dnmat = phylip.readDistMatrix("2YN.dN")
    dsmat = phylip.readDistMatrix("2YN.dS")
    
    
    if saveOutput != "":
        phylip.saveTempDir(cwd, saveOutput)
    else:
        phylip.cleanupTempDir(cwd)
    
    util.toc()
    
    return dnmat, dsmat

