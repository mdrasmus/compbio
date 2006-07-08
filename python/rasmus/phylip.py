import os
import shutil
import sys

import algorithms
import fasta
import matrix
import util



def validateSeq(seqs):
    # ensure sequences are all same size
    sizes = map(len, seqs.values())
    assert max(sizes) == min(sizes), "sequences are not same length"




def checkTempFiles(force):    
    if force:
        os.system("rm -f infile outfile outtree")
    elif os.path.isfile("infile") or \
         os.path.isfile("outfile") or \
         os.path.isfile("outtree"):
        raise Exception("Can't run phylip, 'infile'/'outfile'/'outtree' is in current dir!")


def execPhylip(cmd, args, verbose):
    if verbose:
        util.logger("exec: %s" % cmd)
        util.logger("args: %s" % args)
        assert os.system("""cat <<EOF | %s
%s""" % (cmd, args)) == 0
    else:
        assert os.system("""cat <<EOF | %s >/dev/null 2>&1
%s""" % (cmd, args)) == 0



def cleanup(files=["infile", "outfile", "outtree"]):
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass


def createTempDir():
    directory = os.path.split(util.tempfile(".", "tmpphylip_", ""))[1]
    os.mkdir(directory)
    os.chdir(directory)
    return directory

def cleanupTempDir(directory):
    os.chdir("..")
    assert "/" not in directory
    assert os.path.isdir(directory)
    util.deldir(directory)

def saveTempDir(directory, newname):
    os.chdir("..")
    assert "/" not in directory
    assert os.path.isdir(directory)
    
    if os.path.exists(newname):
        util.deldir(newname)
    
    os.rename(directory, newname)
    

#
# common input/output
#    



def fasta2phylip(out, seqs):
    validateSeq(seqs)
    
    i = 0
    print >>out, len(seqs), len(seqs.values()[0])
    for name in seqs.keys():
        #print >>out, "%9d %s" % (i, seqs[name])
        print >>out, "%8s  %s" % (phylipPadding(str(i), 8), seqs[name])
        i += 1

    return seqs.keys()



def readLogl(filename):
    # parse logl
    logl = None
    for line in file(filename):
        if line.startswith("Ln Likelihood ="):
            logl = float(line.replace("Ln Likelihood =", ""))
    assert logl != None, "could not find logl in outfile"
    
    return logl


def readOutTree(filename, labels, iters=1):
    if iters == 1:
        # parse output
        tree = algorithms.Tree()
        tree.readNewick(filename)
        renameTreeWithNames(tree, labels)
        return tree
    else:
        trees = []
        infile = file("outtree")
        for i in xrange(iters):
            tree = algorithms.Tree()
            tree.readNewick(infile)
            renameTreeWithNames(tree, labels)
            trees.append(tree)
        infile.close()        
        return trees


def writeInTree(filename, tree, labels):
    tree2 = tree.copy()
    renameTreeWithIds(tree2, labels)
    for node in tree2.nodes.values():
        node.dist = 0
    tree2.writeNewick(filename)


def writeBootTrees(filename, trees):
    out = util.openStream(filename, "w")
    for tree in trees:
        tree.writeNewick(out)



def readDistMatrix(filename):
    infile = util.openStream(filename)

    size = int(util.readWord(infile))
    mat = matrix.makeMatrix(size, size)
    names = []
    
    """
     I must be able to read all of these matrices
    
    
      11
_______0    0.00000  0.60810  0.46709  0.57693  0.67485  0.62632  0.64763
            0.67709  0.70192  0.70949  0.68634
_______1    0.60810  0.00000  0.45522  0.49033  0.47842  0.47278  0.47224
            0.47160  0.52655  0.50293  0.49679
_______2    0.46709  0.45522  0.00000  0.57586  0.57433  0.57300  0.56020
            0.57763  0.54225  0.58722  0.58559
_______3    0.57693  0.49033  0.57586  0.00000  0.20713  0.20357  0.21252
            0.46120  0.49081  0.50956  0.49340
_______4    0.67485  0.47842  0.57433  0.20713  0.00000  0.11210  0.13503
            0.45915  0.46692  0.48844  0.47421
_______5    0.62632  0.47278  0.57300  0.20357  0.11210  0.00000  0.10037
            0.45525  0.50959  0.48943  0.49588
_______6    0.64763  0.47224  0.56020  0.21252  0.13503  0.10037  0.00000
            0.46078  0.49727  0.53117  0.51126
_______7    0.67709  0.47160  0.57763  0.46120  0.45915  0.45525  0.46078
            0.00000  0.20980  0.21216  0.20121
_______8    0.70192  0.52655  0.54225  0.49081  0.46692  0.50959  0.49727
            0.20980  0.00000  0.18209  0.13265
_______9    0.70949  0.50293  0.58722  0.50956  0.48844  0.48943  0.53117
            0.21216  0.18209  0.00000  0.08389
______10    0.68634  0.49679  0.58559  0.49340  0.47421  0.49588  0.51126
            0.20121  0.13265  0.08389  0.00000
    
    As well as
    
    
          11
_______0    
_______1    0.60810  
_______2    0.46709  0.45522  
_______3    0.57693  0.49033  0.57586  
_______4    0.67485  0.47842  0.57433  0.20713  
_______5    0.62632  0.47278  0.57300  0.20357  0.11210  
_______6    0.64763  0.47224  0.56020  0.21252  0.13503  0.10037  
_______7    0.67709  0.47160  0.57763  0.46120  0.45915  0.45525  0.46078
_______8    0.70192  0.52655  0.54225  0.49081  0.46692  0.50959  0.49727
            0.20980  
_______9    0.70949  0.50293  0.58722  0.50956  0.48844  0.48943  0.53117
            0.21216  0.18209  
______10    0.68634  0.49679  0.58559  0.49340  0.47421  0.49588  0.51126
            0.20121  0.13265  0.08389  
            
    """
    
    def isName(token):
        try:
            float(token)
            return False
        except:
            return True
    
    
    i = -1
    j = 0
    for line in infile:
        row = line.split()
        if isName(row[0]):
            names.append(row[0])
            row = row[1:]
            i += 1
            j = 0
        
        assert i != -1
        for val in row:
            mat[i][j] = float(val)
            mat[j][i] = float(val)
            j += 1
    
    """
    for i in xrange(size):
        names.append(util.readWord(infile))
        for j in xrange(size):
            mat[i][j] = float(util.readWord(infile))
    """
    
    return names, mat


def writeDistMatrix(mat, labels=None, out=sys.stdout):
    out = util.openStream(out, "w")
    
    out.write("%d\n" % len(mat))
    
    for i in range(len(mat)):
        if labels == None:
            out.write("%8s  " % phylipPadding(str(i)))
        else:
            out.write("%8s  " % labels[i])
        
        for val in mat[i]:
            out.write("%10f" % val)
        out.write("\n")




#
# common conversions
#

def phylipPadding(name, width=8):
    return "_" * (width - len(name)) + name


def renameTreeWithNames(tree, labels):
    names = tree.nodes.keys()
    for name in names:
        if tree.nodes[name].isLeaf():
            num = int(name.replace("_", ""))
            tree.rename(name, labels[num])


def renameTreeWithIds(tree, labels):
    lookup = util.list2lookup(labels)
    
    names = tree.nodes.keys()
    for name in names:
        if tree.nodes[name].isLeaf():
            tree.rename(name, phylipPadding(str(lookup[name])))




#
# tree building programs
#

def align2tree(prog, seqs, verbose=True, force = False, args=None, 
               usertree=None, saveOutput="", 
               bootiter=1,
               seed=1,
               jumble=1):
    validateSeq(seqs)
    cwd = createTempDir()

    util.tic("%s on %d of length %d" % (prog, len(seqs), len(seqs.values()[0])))

    # create input
    labels = fasta2phylip(file("infile", "w"), seqs)

    # initialize default arguments
    if args == None:
        args = "y"
    
    # create user tree if given
    if usertree != None:
        writeInTree("intree", usertree, labels)
        args = "u\n" + args # add user tree option
    
    
    # bootstrap alignment if needed
    if bootiter > 1:
    	execPhylip("seqboot", "r\n%d\ny\n%d" % (bootiter, seed), verbose)
        os.rename("outfile", "infile")    
    	
    	# add bootstrap arguments
    	args = "m\nD\n%d\n%d\n%d\n%s" % (bootiter, seed, jumble, args)
        
    
    # run phylip
    execPhylip(prog, args, verbose)
    
    # parse tree
    if bootiter == 1:
        tree = readOutTree("outtree", labels, bootiter)
        
        # parse likelihood
        tree.data["logl"] = readLogl("outfile")
        
    else:
        trees = readOutTree("outtree", labels, bootiter)
    
    
    if saveOutput != "":
        saveTempDir(cwd, saveOutput)
    else:
        cleanupTempDir(cwd)
    
    util.toc()
    
    
    if bootiter == 1:
        return tree
    else:
        return trees


def bootNeighbor(seqs, iters=100, seed=1, output=None, 
                 verbose=True, force=False):
    validateSeq(seqs)
    cwd = createTempDir()
    util.tic("bootNeighbor on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = fasta2phylip(file("infile", "w"), seqs)
    
    execPhylip("seqboot", "r\n%d\ny\n%d" % (iters, seed), verbose)
    
    os.rename("outfile", "infile")
    execPhylip("protdist", "m\nd\n%d\ny" % iters, verbose)
    
    os.rename("outfile", "infile")
    execPhylip("neighbor", "m\n%d\n%d\ny" % (iters, seed), verbose)

    util.toc()        
    
    # read tree samples
    if output != None:
        os.rename("outtree", "../" + output)
        cleanupTempDir(cwd)
        return labels
    else:
        trees = []
        infile = file("outtree")
        for i in xrange(iters):
            tree = algorithms.Tree()
            tree.readNewick(infile)
            renameTreeWithNames(tree, labels)
            trees.append(tree)
        infile.close()
        cleanupTempDir(cwd)
        return trees


def protpars(seqs, verbose=True, force = False, args="y", 
          usertree=None, saveOutput="", bootiter=1):
    return align2tree("protpars", seqs, 
                      verbose=verbose, 
                      force=force,
                      args=args, 
                      usertree=usertree,
                      saveOutput=saveOutput,
                      bootiter=bootiter)


def proml(seqs, verbose=True, force = False, args="y", 
          usertree=None, saveOutput="", bootiter=1):
    return align2tree("proml", seqs, 
                      verbose=verbose, 
                      force=force,
                      args=args, 
                      usertree=usertree,
                      saveOutput=saveOutput,
                      bootiter=bootiter)


def dnaml(seqs, verbose=True, force = False, args="y", 
          usertree=None, saveOutput="", bootiter=1):
    return align2tree("dnaml", seqs, 
                      verbose=verbose, 
                      force=force,
                      args=args, 
                      usertree=usertree,
                      saveOutput=saveOutput,
                      bootiter=bootiter)


def dnapars(seqs, verbose=True, force = False, args="y", 
          usertree=None, saveOutput="", bootiter=1):
    return align2tree("dnapars", seqs, 
                      verbose=verbose, 
                      force=force,
                      args=args, 
                      usertree=usertree,
                      saveOutput=saveOutput,
                      bootiter=bootiter)



def promlTreelk(aln, tree, verbose=True, force = False, args="u\ny"):
    validateSeq(aln)
    cwd = createTempDir()

    util.tic("proml on %d of length %d" % (len(aln), len(aln.values()[0])))

    # create input
    labels = fasta2phylip(file("infile", "w"), aln)
    writeInTree("intree", tree, labels)
    
    # run phylip
    execPhylip("proml", args, verbose)
    
    # parse logl
    logl = readLogl("outfile")
    
    # parse tree
    tree = readOutTree("outtree", labels)
    
    cleanupTempDir(cwd)
    util.toc()
    
    return logl, tree
    

#
# distance estimation programs
#



    


def protdist(seqs, output=None, verbose=True, force = False, args=None):
    if args == None:
        args = "y"

    validateSeq(seqs)
    cwd = createTempDir()
    util.tic("protdist on %d of length %d" % (len(seqs), len(seqs.values()[0])))
    
    # create input
    labels = fasta2phylip(file("infile", "w"), seqs)
    
    # run phylip
    execPhylip("protdist", args, verbose)
    
    util.toc()    
    
    # parse output
    if output != None:
        os.rename("outfile", "../" + output)
        cleanupTempDir(cwd)
        return labels
    else:
        name, mat = readDistMatrix("outfile")    
        cleanupTempDir(cwd)
        return labels, mat


def dnadist(seqs, output=None, verbose=True, force = False, args=None):
    if args == None:
        args = "y"
    
    validateSeq(seqs)
    cwd = createTempDir()
    util.tic("dnadist on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = fasta2phylip(file("infile", "w"), seqs)
    
    # run phylip
    execPhylip("dnadist", args, verbose)
    
    util.toc()    
    
    # parse output
    if output != None:
        os.rename("outfile", "../" + output)
        cleanupTempDir(cwd)
        return labels
    else:
        name, mat = readDistMatrix("outfile")    
        cleanupTempDir(cwd)
        return labels, mat


def bootNeighbor(seqs, iters=100, seed=1, output=None, 
                 verbose=True, force=False):
    validateSeq(seqs)
    cwd = createTempDir()
    util.tic("bootNeighbor on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = fasta2phylip(file("infile", "w"), seqs)
    
    execPhylip("seqboot", "y\n%d" % seed, verbose)
    
    os.rename("outfile", "infile")
    execPhylip("protdist", "m\nd\n%d\ny" % iters, verbose)
    
    os.rename("outfile", "infile")
    execPhylip("neighbor", "m\n%d\n%d\ny" % (iters, seed), verbose)

    util.toc()        
    
    # read tree samples
    if output != None:
        os.rename("outtree", "../" + output)
        cleanupTempDir(cwd)
        return labels
    else:
        trees = []
        infile = file("outtree")
        for i in xrange(iters):
            tree = algorithms.Tree()
            tree.readNewick(infile)
            renameTreeWithNames(tree, labels)
            trees.append(tree)
        infile.close()
        cleanupTempDir(cwd)
        return trees


def bootProml(seqs, iters = 100, seed = 1, jumble=5, output=None, 
                 verbose=True, force = False):
    validateSeq(seqs)
    cwd = createTempDir()
    util.tic("bootProml on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = fasta2phylip(file("infile", "w"), seqs)
    
    execPhylip("seqboot", "y\n%d" % seed, verbose)
    
    os.rename("outfile", "infile")
    execPhylip("proml", "m\nD\n%d\n%d\n%d\ny" % (iters, seed, jumble), verbose)
    
    util.toc()        
    
    # read tree samples
    if output != None:
        os.rename("outtree", "../" + output)
        cleanupTempDir(cwd)
        return labels
    else:
        trees = []
        infile = file("outtree")
        for i in xrange(iters):
            tree = algorithms.Tree()
            tree.readNewick(infile)
            renameTreeWithNames(tree, labels)
            trees.append(tree)
        infile.close()
        cleanupTempDir(cwd)
        return trees


def consenseFromFile(intrees, verbose=True):
    validateSeq(seqs)
    cwd = createTempDir()
    
    shutil.copy(os.path.join("..", intrees), "intree")
    
    execPhylip("consense", "y", verbose)
    
    tree = algorithms.Tree()
    tree.readNewick("outtree")
    
    cleanupTempDir(cwd)
    return tree

def consense(trees, verbose=True):
    cwd = createTempDir()
    
    writeBootTrees("intree", trees)
    
    execPhylip("consense", "y", verbose)
    
    tree = algorithms.Tree()
    tree.readNewick("outtree")
    
    cleanupTempDir(cwd)
    return tree





"""
old code


def refine(tree, seqs, maxsize, phyfunc, alnfunc):
    trees = algorithms.smallSubtrees(tree, maxsize)
    prog = util.Progress(0, len(tree.leaves()), .01)
    
    for subtree in trees:
        names = tree.leaveNames(subtree.root)
        util.log("subtree size %d" % len(names))
        for i in range(len(names)):
            prog.update()
    
        refineNode(tree, subtree.root, seqs, phyfunc, alnfunc)



def refineNode(tree, node, seqs, phyfunc, alnfunc):
    # get names of sequences to refine
    names = tree.leaveNames(node)
    
    # get outgroup sequence
    if node.parent:
        for child in node.parent.children:
            if child != node:
                while len(child.children) > 0:
                    child = child.children[0]
                outgroup = child.name
                names.append(outgroup)
                break
    else:
        outgroup = None
    
    print "outgroup:", outgroup
    
    if len(names) > 2:
        seqs2 = util.getkeys(seqs, names)
        aln   = alnfunc(seqs2)    
        tree2 = phyfunc(aln, outgroup=outgroup)
        
        if outgroup:
            tree2.remove(tree2.nodes[outgroup])
            if len(tree2.root.children) == 1:
                tree2.root = tree2.root.children[0]

        # only replace tree if a new one is created
        if len(tree2.nodes) != 0:
            parent = node.parent
            index = parent.children.index(node)
            tree.removeTree(node)
            tree.addTree(parent, tree2)
            
            # ensure new subtree appears in same spot as old subtree
            parent.children = parent.children[:index] + [parent.children[-1]] + \
                              parent.children[index:-1]
"""





# testing
if __name__ == "__main__":
    seqs = fasta.readFasta("test/dna-align.fa")
    del seqs["target"]

    reload(algorithms)
    tree = protpars(seqs, force=True, verbose=False)
    tree.write()
