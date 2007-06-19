"""

 python wrapper for CLUTO clustering toolkit

"""


# python libs
import os

# rasmus libs
from rasmus import matrix
from rasmus import util
from rasmus import treelib


def writeDenseMatrix(filename, mat):
    """Write a CLUTO formatted dense matrix file"""
    
    out = util.openStream(filename, "w")
    
    print >>out, len(mat), len(mat[0])
    
    for row in mat:
        for val in row[:-1]:
            out.write(str(val) + "\t")
        out.write(str(row[-1]) + "\n")


def writeSquareMatrix(filename, mat):
    """Write a CLUTO formatted square matrix file"""
    
    out = util.openStream(filename, "w")
    
    print >>out, len(mat)
    
    for row in mat:
        for val in row[:-1]:
            out.write(str(val) + "\t")
        out.write(str(row[-1]) + "\n")            


def readDenseMatrix(filename):
    """Read a CLUTO formatted dense matrix file"""
    
    infile = util.openStream(filename)
    
    header = map(int, infile.next().split())
    
    if len(header) == 1:
        nrows = ncols = header[0]
    elif len(header) == 2:
        nrows, ncols = header
    else:
        raise Exception("Not a CLUTO dense matrix file")
    
    mat = []
    for line in infile:
        mat.append(map(float, line.split()))
        assert len(mat[-1]) == ncols, Exception("File has wrong number of columns")
    
    assert len(mat) == nrows, Exception("File has wrong number of rows")
    
    return mat

    

# TODO: add saveOutput
def cluster(mat, nclusters=10, prog="vcluster", options="", verbose=False):
    """Cluster a matrix into 'nclusters'
        
       mat  - a matrix or a matrix filename
       prog - the program to use:
              "vcluster" will treat the matrix as a series of row vectors
              "scluster" will treat the matrix as a square matrix of
                         similarities between objects

       returns a partition id list (partid) which indicates the partition id
       of each item (row) in the input matrix (mat)
    """
    
    
    if isinstance(mat, str):
        matfile = mat
        tmpfile = None
    else:
        tmpfile = util.tempfile("./", "tmpcluto_", ".mat")
        matfile = tmpfile
        if prog == "scluster":
            writeSquareMatrix(tmpfile, mat)
        else:
            writeDenseMatrix(tmpfile, mat)
    
    partfile = tmpfile + ".clustering.%d" % nclusters
    
    
    if verbose:
        os.system("%s %s %d %s" % (prog, matfile, nclusters, options))
    else:
        os.system("%s %s %d %s > /dev/null" % (prog, matfile, nclusters, options))
    
    partids = util.readInts(partfile)
    
    if tmpfile != None:
        os.remove(tmpfile)
    os.remove(partfile)
    
    return partids



# TODO: add saveOutput
def clusterTree(mat, nclusters=10, prog="vcluster", options="", verbose=False):
    """Cluster a matrix into 'nclusters'
        
       mat  - a matrix or a matrix filename
       prog - the program to use:
              "vcluster" will treat the matrix as a series of row vectors
              "scluster" will treat the matrix as a square matrix of
                         similarities between objects

       returns a hierarchical tree representing the clustering
    """
    
    if isinstance(mat, str):
        matfile = mat
        tmpfile = None
    else:
        tmpfile = util.tempfile("./", "tmpcluto_", ".mat")
        matfile = tmpfile
        if prog == "scluster":
            writeSquareMatrix(tmpfile, mat)
        else:
            writeDenseMatrix(tmpfile, mat)
        
    # determine files
    partfile = matfile + ".clustering.%d" % nclusters
    treefile = matfile + ".tree"

    
    if verbose:
        os.system("%s %s %d %s -fulltree" % (prog, matfile, nclusters, options))
    else:
        os.system("%s %s %d %s -fulltree > /dev/null" % (prog, matfile, nclusters, options))
    
    tree = treelib.Tree()
    tree.readParentTree(treefile)
    
    if tmpfile != None:
        os.remove(tmpfile)
    os.remove(partfile)
    os.remove(treefile)
    
    return tree


def reorderTree(tree, mat, prog="vcluster"):
    """
    Reorders a tree such that leaves appears in a visually pleasing order

    returns a permutation list
    """
    
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

#
# TODO: add a function that flips a tree based on a perm
#


def reorderPartids(partids, mat, prog="vcluster"):
    """
    Reorders the rows of a matrix such that rows appears in a
    visually pleasing order

    returns a permutation list
    """
    
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
    """Trivially produce a permutation list from a partid list"""
    
    perm = range(len(partids))
    perm.sort(key=lambda x: partids[x])
    
    return perm





##############################
# testing code
#
if __name__ == "__main__":
    from rasmus.common import *
    import random
    import summon
    from summon import matrix as summatrix
    
    mat = []
    nclusters = 3
    nrows = 100
    ncols = 30
    clustersize = ncols // nclusters

    # make data for cluster
    for i in xrange(nrows):
        k = i % nclusters
        
        cluster1 = []
        for j in xrange(ncols):
            if j // clustersize == k:
                cluster1.append(2)
            else:
                cluster1.append(0)

        mat.append([])

        for j in xrange(ncols):
            mat[-1].append(random.normalvariate(cluster1[j], 1))


    # randomize matrix
    random.shuffle(mat)
    
    # display shuffled matrix
    #heatmap(mat)
    viewer1 = summatrix.DenseMatrixViewer(mat)
    viewer1.show()

    #########################
    # cluster by vcluster
    partids = cluster(mat, 2)
    perm = reorderPartids(partids, mat)

    # apply permutation to matrix
    mat2 = mget(mat, perm)

    # display clustered matrix
    viewer2 = summatrix.DenseMatrixViewer(mat2)
    viewer2.show()

if 0:    
    #########################
    # cluster by scluster
    import stats, matrix
    simmat = stats.corrmatrix(mat)
    partids = cluster(simmat, 2, prog="scluster")
    perm = reorderPartids(partids, mat)

    # apply permutation to matrix
    mat2 = mget(mat, perm)
    simmat2 = matrix.submatrix(simmat, perm, perm)

    # display similarity matrix
    viewer3 = summatrix.DenseMatrixViewer(simmat, cuttoff=.4)
    viewer3.show()

    # display clustered similarity matrix
    viewer4 = summatrix.DenseMatrixViewer(simmat2, cuttoff=.4)
    viewer4.show()



    


