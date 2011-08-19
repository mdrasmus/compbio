"""

    PHYLIP wrapper for python
    
    author: Matt Rasmussen
    date:   2/4/2007


The following programs should be in your PATH

    dnaml       -- Maximum likelihood (nucleotide)
    proml       -- Maximum likelihood (peptide)
    dnapars     -- Parsimony (nucleotide)
    protpars    -- Parsimony (peptide)
    neighbor    -- Neighbor Joining
    dnadist     -- Distance estimation (nucleotide)
    protdist    -- Distance estimation (peptide)
    seqboot     -- Sequence bootstrapping
    consense    -- Consensus tree building

"""


# python imports
import os
import shutil
import sys

# rasmus imports
from rasmus import util
from rasmus import treelib


# compbio imports
from . import fasta


#=============================================================================
# managing input, output, and execution of PHYLIP-like programs
#

def validate_seqs(seqs):
    """Ensures sequences are all same size"""
    
    sizes = map(len, seqs.values())
    assert util.equal(* sizes), "sequences are not same length"


def check_temp_files(force=False):
    """Ensure PHYLIP tempfiles do not already exist in current directory"""
    
    if force:
        os.system("rm -f infile outfile outtree")
    elif os.path.isfile("infile") or \
         os.path.isfile("outfile") or \
         os.path.isfile("outtree"):
        raise Exception("Can't run phylip, 'infile'/'outfile'/'outtree' is in current dir!")


def exec_phylip(cmd, args, verbose=False):
    """Execute a phylip-like program that expects arguments from stdin"""

    util.logger("exec: %s" % cmd)
    util.logger("args: %s" % args)


    if verbose:
        util.logger("exec: %s" % cmd)
        util.logger("args: %s" % args)
        assert os.system("""cat <<EOF | %s
%s""" % (cmd, args)) == 0
    else:
        assert os.system("""cat <<EOF | %s >/dev/null 2>&1
%s""" % (cmd, args)) == 0




def create_temp_dir(prefix="tmpphylip_"):
    """Create a temporary directory for executing PHYLIP"""
    
    directory = os.path.split(util.tempfile(".", prefix, ""))[1]
    os.mkdir(directory)
    os.chdir(directory)
    return directory


def cleanup_temp_dir(directory):
    """Exit and delete a temporary directory for executing PHYLIP"""

    os.chdir("..")
    assert "/" not in directory
    assert os.path.isdir(directory)
    util.deldir(directory)


def save_temp_dir(directory, newname):
    """Exit and save a temporary directory for executing PHYLIP"""

    os.chdir("..")
    assert "/" not in directory
    assert os.path.isdir(directory)
    
    if os.path.exists(newname):
        util.deldir(newname)
    
    os.rename(directory, newname)



#=============================================================================
# common input/output
#    


def read_phylip_align(filename):
    """
    Read a PHYLIP alignment.  Can be interleaved or not.
    
    returns a FastaDict object.
    """
    
    infile = util.open_stream(filename)
    
    seqs = fasta.FastaDict()
    
    # read sequences and length
    nseq, seqlen = infile.next().split()
    nseq = int(nseq)
    
    i = 0
    first = True
    names = []
    
    # parse remaining lines
    for line in infile:
        line = line.rstrip()
    
        if len(line) > 0:
            if first:
                name = line[:10].strip()
                seq = line[10:].strip().replace(" ", "")
                names.append(name)
            else:
                seq = line.strip().replace(" ", "")
                name = names[i]
            i += 1                
        
            if not name in seqs:
                seqs[name] = seq
            else:
                seqs[name] += seq
        else:
            i = 0
            first = False
    return seqs


def write_phylip_align(out, seqs, strip_names=True):
    """
    Write a PHYLIP alignment
    """
    
    validate_seqs(seqs)
    
    if strip_names:
        print >>out, len(seqs), len(seqs.values()[0])
        for i, name in enumerate(seqs.keys()):
            print >>out, "%8s  %s" % (phylip_padding(str(i), 8), seqs[name])
    else:
        print >>out, len(seqs), len(seqs.values()[0])
        for i, name in enumerate(seqs.keys()):
            if name <= 8:
                print >>out, "%8s  %s" % (name, seqs[name])
            else:
                print >>out, "%s  %s" % (name, seqs[name])

    return seqs.keys()


def read_logl(filename):
    # parse logl
    logl = None
    for line in file(filename):
        if line.startswith("Ln Likelihood ="):
            logl = float(line.replace("Ln Likelihood =", ""))
    assert logl != None, "could not find logl in outfile"
    
    return logl


def read_out_tree(filename, labels, iters=1):
    infile = file(filename)
    
    # skip any numbers that may appear on the first line
    line = infile.readline()
    if not line[0].isdigit():
        # reopen file
        infile = file(filename)

    
    if iters == 1:
        # parse output
        tree = treelib.Tree()
        tree.read_newick(infile)
        rename_tree_with_name(tree, labels)
        return tree
    else:
        trees = []
        for i in xrange(iters):
            tree = treelib.Tree()
            tree.read_newick(infile)
            rename_tree_with_name(tree, labels)
            trees.append(tree)
        infile.close()        
        return trees


def write_in_tree(filename, tree, labels):
    tree2 = tree.copy()
    rename_tree_with_ids(tree2, labels)
    for node in tree2.nodes.values():
        node.dist = 0
    tree2.write(filename)


def write_boot_trees(filename, trees, counts=None):
    out = util.open_stream(filename, "w")
    
    if counts == None:
        counts = [1] * len(trees)
    
    for tree, count in zip(trees, counts):
        for i in range(count):
            out.write(tree.get_one_line_newick() + "\n")



def read_dist_matrix(filename):
    infile = util.open_stream(filename)

    size = int(util.read_word(infile))
    mat = util.make_matrix(size, size)
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
        if len(row) == 0:
            continue
        
        if isName(row[0]):
            names.append(row[0])
            row = row[1:]
            i += 1
            j = 0
        
        assert i != -1
        for val in row:
            if val == "nan" or val == "inf":
                val = None
            else:
                val = float(val)
            
            
            mat[i][j] = val
            mat[j][i] = val
            j += 1
    
    # remove nasty infinities
    top = util.max2(mat)
    for i in range(size):
        for j in range(size):
            if mat[i][j] == None:
                mat[i][j] = 10 * top
    
    """
    for i in xrange(size):
        names.append(util.read_word(infile))
        for j in xrange(size):
            mat[i][j] = float(util.read_word(infile))
    """
    
    return names, mat


def write_dist_matrix(mat, labels=None, out=sys.stdout):
    out = util.open_stream(out, "w")
    
    out.write("%d\n" % len(mat))
    
    for i in range(len(mat)):
        if labels == None:
            out.write("%8s  " % phylip_padding(str(i)))
        else:
            out.write("%8s  " % labels[i])
        
        for val in mat[i]:
            out.write("%10f " % val)
        out.write("\n")





#=============================================================================
# common conversions
#

def phylip_padding(name, width=8):
    return "_" * (width - len(name)) + name


def rename_tree_with_names(tree, labels):
    names = tree.nodes.keys()
    for name in names:
        if tree.nodes[name].is_leaf():
            num = int(name.replace("_", ""))
            tree.rename(name, labels[num])


def rename_tree_with_ids(tree, labels):
    lookup = util.list2lookup(labels)
    
    names = tree.nodes.keys()
    for name in names:
        if tree.nodes[name].is_leaf():
            tree.rename(name, phylip_padding(str(lookup[name])))



#=============================================================================
# tree building programs
#

def align2tree(prog, seqs, verbose=True, force = False, args=None, 
               usertree=None, saveOutput="", 
               bootiter=1,
               seed=1,
               jumble=1):
    validate_seqs(seqs)
    cwd = create_temp_dir()

    util.tic("%s on %d of length %d" % (prog, len(seqs), len(seqs.values()[0])))

    # create input
    labels = write_phylip_align(file("infile", "w"), seqs)
    util.write_list(file("labels", "w"), labels)

    # initialize default arguments
    if args == None:
        args = "y"
    
    # create user tree if given
    if usertree != None:
        write_in_tree("intree", usertree, labels)
        args = "u\n" + args # add user tree option
    
    
    # bootstrap alignment if needed
    if bootiter > 1:
    	exec_phylip("seqboot", "r\n%d\ny\n%d" % (bootiter, seed), verbose)
        os.rename("outfile", "infile")    
    	
    	# add bootstrap arguments
    	args = "m\nD\n%d\n%d\n%d\n%s" % (bootiter, seed, jumble, args)
        
    
    # run phylip
    exec_phylip(prog, args, verbose)
    
    # check for PHYLIP GIVE UP
    if is_phylip_give_up("outfile"):
        tree = treelib.Tree()
        tree.make_root()
        
        # make star tree
        for key in seqs:
            tree.add_child(tree.root, treelib.TreeNode(key))
        
    else:
        # parse tree
        if bootiter == 1:
            tree = read_out_tree("outtree", labels, bootiter)

            # parse likelihood
            if prog in ["dnaml", "proml"]:
                tree.data["logl"] = read_logl("outfile")

        else:
            trees = read_out_tree("outtree", labels, bootiter)
    
    
    if saveOutput != "":
        save_temp_dir(cwd, saveOutput)
    else:
        cleanup_temp_dir(cwd)
    
    util.toc()
    
    
    if bootiter == 1:
        return tree
    else:
        return trees


def is_phylip_give_up(filename):
    for line in file(filename):
        if "0 trees in all found" in line:
            return True
    return False
    




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



def proml_treelk(aln, tree, verbose=True, force = False, args="u\ny"):
    validate_seqs(aln)
    cwd = create_temp_dir()

    util.tic("proml on %d of length %d" % (len(aln), len(aln.values()[0])))

    # create input
    labels = write_phylip_align(file("infile", "w"), aln)
    write_in_tree("intree", tree, labels)
    
    # run phylip
    exec_phylip("proml", args, verbose)
    
    # parse logl
    logl = read_logl("outfile")
    
    # parse tree
    tree = read_out_tree("outtree", labels)
    
    cleanup_temp_dir(cwd)
    util.toc()
    
    return logl, tree


def draw_tree(tree, plotfile, verbose=False, args=None, saveOutput = ""):
    cwd = create_temp_dir()
    
    fontfile = os.popen("which font4", "r").read().rstrip()
    
    # create input
    tree.write("intree")
    
    # initialize default arguments
    if args == None:
        args = "%s\nv\nn\ny" % fontfile
    
    # run phylip
    exec_phylip("drawgram", args, verbose)
    
    os.rename("plotfile", "../" + plotfile)
    
    if saveOutput != "":
        save_temp_dir(cwd, saveOutput)
    else:
        cleanup_temp_dir(cwd)
    
    

#=============================================================================
# distance estimation programs
#



    


def protdist(seqs, output=None, verbose=True, force = False, args=None):
    if args == None:
        args = "y"

    validate_seqs(seqs)
    cwd = create_temp_dir()
    util.tic("protdist on %d of length %d" % (len(seqs), len(seqs.values()[0])))
    
    # create input
    labels = write_phylip_align(file("infile", "w"), seqs)
    
    # run phylip
    exec_phylip("protdist", args, verbose)
    
    util.toc()    
    
    # parse output
    if output != None:
        os.rename("outfile", "../" + output)
        cleanup_temp_dir(cwd)
        return labels
    else:
        name, mat = read_dist_matrix("outfile")    
        cleanup_temp_dir(cwd)
        return labels, mat


def dnadist(seqs, output=None, verbose=True, force = False, args=None):
    if args == None:
        args = "y"
    
    validate_seqs(seqs)
    cwd = create_temp_dir()
    util.tic("dnadist on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = write_phylip_align(file("infile", "w"), seqs)
    
    # run phylip
    exec_phylip("dnadist", args, verbose)
    
    util.toc()    
    
    # parse output
    if output != None:
        os.rename("outfile", "../" + output)
        cleanup_temp_dir(cwd)
        return labels
    else:
        name, mat = read_dist_matrix("outfile")    
        cleanup_temp_dir(cwd)
        return labels, mat


def correct_dist_matrix(distmat, maxdist=40, fardist=None):
    """remove -1 and extremely large distances (>maxdist), replace them with 
       fatdist (defaults to maximum distance in matrix)"""

    if fardist == None:
        fardist = 0

        for row in distmat:
            for x in row:
                if x < maxdist:
                    fardist = max(fardist, x)
    
    distmat2 = []
    for row in distmat:
        distmat2.append([])
        for x in row:
            if x == -1 or x > maxdist:
                distmat2[-1].append(fardist)
            else:
                distmat2[-1].append(x)
    
    return distmat2
    

def boot_neighbor(seqs, iters=100, seed=None, output=None, 
                 verbose=True, force=False):
    
    if seed == None:
        seed = random.randInt(0, 1000) * 2 + 1
    
    validate_seqs(seqs)
    cwd = create_temp_dir()
    util.tic("boot_neighbor on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = write_phylip_align(file("infile", "w"), seqs)
    
    exec_phylip("seqboot", "r\n%d\ny\n%d" % (iters, seed), verbose)
    
    os.rename("outfile", "infile")
    exec_phylip("protdist", "m\nd\n%d\ny" % iters, verbose)
    
    os.rename("outfile", "infile")
    exec_phylip("neighbor", "m\n%d\n%d\ny" % (iters, seed), verbose)

    util.toc()        
    
    # read tree samples
    if output != None:
        os.rename("outtree", "../" + output)
        cleanup_temp_dir(cwd)
        return labels
    else:
        trees = []
        infile = file("outtree")
        for i in xrange(iters):
            tree = treelib.Tree()
            tree.read_newick(infile)
            rename_tree_with_name(tree, labels)
            trees.append(tree)
        infile.close()
        cleanup_temp_dir(cwd)
        return trees



def boot_proml(seqs, iters = 100, seed = 1, jumble=5, output=None, 
               verbose=True, force = False):
    validate_seqs(seqs)
    cwd = create_temp_dir()
    util.tic("bootProml on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = write_phylip_align(file("infile", "w"), seqs)
    
    exec_phylip("seqboot", "y\n%d" % seed, verbose)
    
    os.rename("outfile", "infile")
    exec_phylip("proml", "m\nD\n%d\n%d\n%d\ny" % (iters, seed, jumble), verbose)
    
    util.toc()        
    
    # read tree samples
    if output != None:
        os.rename("outtree", "../" + output)
        cleanup_temp_dir(cwd)
        return labels
    else:
        trees = []
        infile = file("outtree")
        for i in xrange(iters):
            tree = treelib.Tree()
            tree.read_newick(infile)
            rename_tree_with_names(tree, labels)
            trees.append(tree)
        infile.close()
        cleanup_temp_dir(cwd)
        return trees


def consense_from_file(intrees, verbose=True, args="y"):

    # read all trees
    trees = util.open_stream(intrees).readlines()
    ntrees = len(trees)

    cwd = create_temp_dir()
    out = open("intree", "w")
    for tree in trees:
        out.write(tree)
    out.close()
    
    exec_phylip("consense", args, verbose)
    
    tree = treelib.read_tree("outtree")
    
    cleanup_temp_dir(cwd)
    return tree, ntrees


def consense(trees, counts=None, verbose=True, args="y"):
    cwd = create_temp_dir()
    
    write_boot_trees("intree", trees, counts=counts)
    
    exec_phylip("consense", args, verbose)
    
    tree = treelib.Tree()
    tree.read_newick("outtree")
    
    cleanup_temp_dir(cwd)
    return tree




# testing
if __name__ == "__main__":
    seqs = fasta.read_fasta("test/dna-align.fa")
    del seqs["target"]

    tree = protpars(seqs, force=True, verbose=False)
    tree.write()
