
from rasmus import util
from rasmus import treelib
import os

def writeDenseMatrix(filename, mat):
    out = util.openStream(filename, "w")
    
    print >>out, len(mat), len(mat[0])
    
    for row in mat:
        for val in row[:-1]:
            out.write(str(val) + "\t")
        out.write(str(row[-1]) + "\n")


def writeSquareMatrix(filename, mat):
    out = util.openStream(filename, "w")
    
    print >>out, len(mat)
    
    for row in mat:
        for val in row[:-1]:
            out.write(str(val) + "\t")
        out.write(str(row[-1]) + "\n")            
    

def cluster(mat, nclusters=10, prog="vcluster", options="", verbose=False):
    
    tmpfile = util.tempfile("./", "tmpcluto_", ".mat")
    partfile = tmpfile + ".clustering.%d" % nclusters
    
    if prog == "scluster":
        writeSquareMatrix(tmpfile, mat)
    else:
        writeDenseMatrix(tmpfile, mat)
    
    if verbose:
        os.system("%s %s %d %s" % (prog, tmpfile, nclusters, options))
    else:
        os.system("%s %s %d %s > /dev/null" % (prog, tmpfile, nclusters, options))
    
    partids = util.readInts(partfile)
    
    os.remove(tmpfile)
    os.remove(partfile)
    
    return partids


def clusterTree(mat, nclusters=10, prog="vcluster", options="", verbose=False):
    
    # determine files
    tmpfile = util.tempfile("./", "tmpcluto_", ".mat")
    partfile = tmpfile + ".clustering.%d" % nclusters
    treefile = tmpfile + ".tree"
    
    # write input matrix
    if prog == "scluster":
        writeSquareMatrix(tmpfile, mat)
    else:
        writeDenseMatrix(tmpfile, mat)
    
    if verbose:
        os.system("%s %s %d %s -fulltree" % (prog, tmpfile, nclusters, options))
    else:
        os.system("%s %s %d %s -fulltree > /dev/null" % (prog, tmpfile, nclusters, options))
    
    tree = treelib.Tree()
    tree.readParentTree(treefile)
    
    os.remove(tmpfile)
    os.remove(partfile)
    os.remove(treefile)
    
    return tree


def reorderTree(tree, mat, prog="vcluster"):
    matfile = util.tempfile("./", "tmpcluto_", ".mat")
    treefile = util.tempfile("./", "tmpcluto_", ".tree")
    
    if prog == "scluster":
        writeSquareMatrix(matfile, mat)
    else:
        writeDenseMatrix(matfile, mat)
    tree.writeParentTree(treefile, map(str, range(len(tree.leaves()))))

    permfile = os.popen("cluto_reorder_tree %s %s 2>/dev/null" % 
                        (matfile, treefile))
    perm = util.readInts(permfile)
    
    os.remove(matfile)
    os.remove(treefile)
    
    return perm


def reorderPartids(partids, mat, prog="vcluster"):
    matfile = util.tempfile("./", "tmpcluto_", ".mat")
    partidsfile = util.tempfile("./", "tmpcluto_", ".partids")
    
    if prog == "scluster":
        writeSquareMatrix(matfile, mat)
    else:
        writeDenseMatrix(matfile, mat)
    util.writeVector(partidsfile, partids)

    permfile = os.popen("cluto_reorder %s %s full 2>/dev/null" % 
                        (matfile, partidsfile))
    perm = util.readInts(permfile)
    
    os.remove(matfile)
    os.remove(partidsfile)
    
    return perm

def partids2perm(partids):
    perm = range(len(partids))
    perm.sort(key=lambda x: partids[x])
    
    return perm



