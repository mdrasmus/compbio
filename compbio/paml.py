
# python libs
import os, re

# rasmus libs
from rasmus import util, tablelib, treelib

# compbio imports
from . import alignlib, seqlib
from . import fasta, phylip, nexus


def removeStopCodons(seq, stops=None):
    assert len(seq) % 3 == 0
    
    # ensure stop codons are defined
    if stops is None:
        stops = set(["TAA", "TAG", "TGA"])          
    
    seq2 = []
    
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        if codon in stops:
            seq2.append("NNN")
        else:
            seq2.append(codon)
    return "".join(seq2)


def dndsMatrix(seqs, saveOutput="", verbose=False, safe=True):
    
    if safe:
        seqs = alignlib.mapalign(seqs, valfunc=removeStopCodons)
    
    phylip.validate_seqs(seqs)
    cwd = phylip.create_temp_dir()

    util.tic("yn00 on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = phylip.write_phylip_align(file("seqfile.phylip", "w"), seqs)
    util.write_list(file("labels", "w"), labels)    
    
    # create control file
    out = file("yn00.ctl", "w")
    print >>out, "seqfile = seqfile.phylip"
    print >>out, "outfile = outfile"
    out.close()
    
    # run yn00
    if verbose:
        os.system("yn00 yn00.ctl")
    else:
        os.system("yn00 yn00.ctl > /dev/null")
    
    try:
        dnmat = phylip.read_dist_matrix("2YN.dN")
        dsmat = phylip.read_dist_matrix("2YN.dS")
    except:
        # could not make distance matrix
        if safe:
            # make dummy matrix
            dnmat = labels, [[0] * len(labels)] * len(labels)
            dsmat = labels, [[0] * len(labels)] * len(labels)
        else:
            raise Exception("could not read dn or ds matrices")
    
    if saveOutput != "":
        phylip.save_temp_dir(cwd, saveOutput)
    else:
        phylip.cleanup_temp_dir(cwd)
    
    util.toc()
    
    return dnmat, dsmat


def pamp(seqs, tree, seqtype="dna", saveOutput="", verbose=False, safe=True):
    
    if safe and seqtype == "dna":
        seqs = alignlib.mapalign(seqs, valfunc=removeStopCodons)
    
    phylip.validate_seqs(seqs)
    cwd = phylip.create_temp_dir()

    util.tic("pamp on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    nex = nexus.Nexus("align", "w")
    nex.write_matrix(seqs.keys(), seqs.values(), seqtype, seqs.alignlen())
    nex.close()
    treefile = open("tree", "w")
    treefile.write("%d 1\n" % len(tree.leaves()))
    tree.write(treefile, writeData=lambda x: "")
    treefile.close()
    
    # create control file
    out = file("pamp.ctl", "w")
    print >>out, "seqfile = align"
    print >>out, "treefile = tree"
    print >>out, "outfile = out"    

    if seqtype == "dna":
        print >>out, "seqtype = 0"
    elif seqtype == "pep":
        print >>out, "seqtype = 2"
    else:
        raise Exception("unknown seqtype '%s'" % seqtype)
    print >>out, "ncatG = 8"
    print >>out, "nhomo = 0"
    out.close()
    
    # run pamp
    if verbose:
        os.system("pamp paml.ctl")
    else:
        os.system("pamp paml.ctl > /dev/null")

    res = PamlResults("out")
    aln = res.getPampReconstruction()
    aln.write("recon.mfa")
    tree2 = res.getBranchNames()
    renameTreeAlign(tree2, aln)
    
    if saveOutput != "":
        phylip.save_temp_dir(cwd, saveOutput)
    else:
        phylip.cleanup_temp_dir(cwd)
    
    util.toc()

    return tree2, aln
    

def renameTreeAlign(tree, align):
    """Rename the ancestral nodes of a tree"""

    # rename leaves
    for leaf in tree.leaves():
        tree.rename(leaf.name, align.keys()[int(leaf.name)-1])

    names = align.keys()
    name_lookup = util.list2lookup(names)
    lookup = dict(align)
    align.clear()
    order = {}

    # rename ancestral nodes by preorder
    next = [1]
    def walk(node):
        oldname = node.name
        if not node.is_leaf():
            tree.rename(node.name, next[0])
            next[0] += 1
            align[str(node.name)] = lookup[str(oldname)]
            order[str(node.name)] = name_lookup[str(oldname)]
        else:
            align[node.name] = lookup[node.name]
            order[node.name] = name_lookup[node.name]
        for child in node.children:
            walk(child)
    walk(tree.root)

    align.names.sort(key=lambda x: order[x])


def countTreeSub(tree, align):

    for node in tree:
        if node == tree.root:
            node.dist = 0
        else:
            seq1 = align[str(node.name)]
            seq2 = align[str(node.parent.name)]

            node.dist = util.count(lambda i: seq1[i] != seq2[i],
                                   xrange(align.alignlen()))

        

def align2tree(seq, seqtype="dna", 
               saveOutput="", usertree=None, verbose=False, safe=True):
    if safe:
        seqs = alignlib.mapalign(seqs, valfunc=removeStopCodons)
    
    # create input
    labels = phylip.write_phylip_align(file("seqfile.phylip", "w"), seqs)
    
    # set verbose level
    if verbose:
        verbose_num = 1
    else:
        verbose_num = 0
        
    if usertree:
        runmode = 0
    else:
        runmode = 2
    
    if seqtype == "dna":
        seqtype = 1
    else:
        seqtype = 2
    
    # create control file
    out = file("yn00.ctl", "w")
    print >>out, util.evalstr("""
    seqfile = seqfile.phylip
    treefile = out.tree
    outfile = out.stats
    
    noisy = 9
    verbose = ${verbose_num}
    runmode = ${runmode}
    
    seqtype = ${seqtype}
    CodonFreq = 2
    
    """)


def read_tree_data(node, data):
    """Note: PAML does not use bootstraps"""

    tokens = data.split(":", 1)

    if len(tokens) == 1:
        other = ""
        dist = tokens[0]
    else:
        other, dist = tokens
        
    # set branch length
    node.dist = float(dist)
            
    if len(other) > 0:

        # read branch labels: integers starting from 0
        if "#" in other:
            other, label = other.split("#")
            node.data["label"] = int(label)


def write_tree_data(node):

    text = ""

    if "label" in node.data:
        text += " # %d " % node.data["label"]

    text += ":%f" % node.dist
    return text
    

def read_tree(treefile):
    """Read a tree for PAML, includes branch labels"""

    tree = treelib.Tree()
    tree.read_newick(treefile, readData=read_tree_data)
    return tree


def write_tree(treefile, tree):
    tree.write(treefile, writeData=write_tree_data)


class PamlResults (object):
    def __init__(self, resultsFile):
        self._lines = util.open_stream(resultsFile).readlines()

    def setupLines(self, lines):
        if lines is None:
            return iter(self._lines)
        else:
            return iter(lines)

        
    def iterBebLines(self, lines=None):
        """Iterate over the BEB section of PamlResults"""

        lines = self.setupLines(lines)
        
        for line in lines:
            if line.startswith("Bayes Empirical Bayes"):
                # skip next four lines
                stream.next()
                stream.next()
                stream.next()
                stream.next()
                break
        
        for line in lines:
            if line == "\n":
                break
            yield line
        
    def getBeb(self, lines=None):
        """Iterate over a parsed matrix of BEB data"""
        
        if lines is None:
            lines = self.iterBebLines()
        
        for line in lines:
            pos, aa, prob, omega, plusmin, omega_stdev = line.rstrip().split()
            yield [int(pos), aa, float(prob.replace("*", "")),
                        float(omega), float(omega_stdev)]
    
    def getBebTable(self, lines=None):
        """Returns a table of BEB data"""
        
        return tablelib.Table(list(self.getBeb(lines)), 
                              headers=["pos", "aa", "prob", "omega",
                                       "omega_stdev"],
                              types={"pos": int,
                                     "aa": str,
                                     "prob": float,
                                     "omega": float,
                                     "omega_stdev": float})
    
    
    def getM0_dnds(self, lines=None, header=True):
        """Return the omega value from the Model 0 section"""
        
        omega = None
        dn = None
        ds = None
        dt = None
        kappa = None

        lines = self.setupLines(lines)        
        
        if header:
            for line in lines:
                if line.startswith("Model 0"):
                    break
            else:
                # paml file does not contain model 0
                return None, None, None, None, None
        
        for line in lines:
            if line.startswith("omega (dN/dS) ="):
                omega = float(line.rstrip().split()[-1])
            
            elif line.startswith("kappa (ts/tv) ="):
                kappa = float(line.rstrip().split()[-1])
            
            elif line.startswith("tree length for dN:"):
                dn = float(line.split(":")[1])
            
            elif line.startswith("tree length for dS:"):
                ds = float(line.split(":")[1])
                
            elif line.startswith("tree length ="):
                dt = float(line.split("=")[1])
            
            # quit if entering next model
            elif line.startswith("Model"):
                break
                
        assert None not in (omega, dn, ds, dt, kappa)
        return omega, dn, ds, dt, kappa
    
    
    def get_ng_dnds(self, lines=None):
        """Return the Nei & Gojobori 1986 dN, dS  matrices"""
        
        genes1 = []
        genes2 = []
        omegas = []
        dn = []
        ds = []

        lines = self.setupLines(lines)
        
        for line in lines:
            if line.startswith("Nei & Gojobori 1986"):
                break
        else:
            # paml file does not contain Nei Gojobori matrix
            return None, None, None, None, None
        
        # skip to first blank line
        for line in lines:
            if len(line) <= 1:
                break
        
        genes = []
        pattern = re.compile("[ \t\(\)]+")
        for line in lines:
            # quit on first blank line
            if len(line) <= 1:
                break
            
            tokens = re.split(pattern, line.rstrip("\n\t )"))
            if len(tokens) != 1 + 3 * len(genes):
                # early truncation
                return None, None, None, None, None
            
            genes.append(tokens[0])           
            for i in xrange(0, 3*(len(genes) - 1), 3):
                genes1.append(tokens[0])
                genes2.append(genes[i/3])
                omegas.append(float(tokens[1+i]))
                dn.append(float(tokens[1+i+1]))
                ds.append(float(tokens[1+i+2]))
        
        return genes1, genes2, omegas, dn, ds


    def getM1_lnl(self, lines=None):

        # nearly neutral
        lines = self.setupLines(lines)
        for line in lines:
            if line.startswith("Model 1:"):
                break

        for line in lines:
            if line.startswith("lnL"):
                return float(line.split(":")[3].strip().split()[0])

        return None

    def getM2_lnl(self, lines=None):

        # positive selection
        lines = self.setupLines(lines)        
        for line in lines:
            if line.startswith("Model 2:"):
                break

        for line in lines:
            if line.startswith("lnL"):
                return float(line.split(":")[3].strip().split()[0])

        return None


    def get_lnl(self, lines=None):
        """Get the tree likelihood"""

        lines = self.setupLines(lines)
        for line in lines:
            if line.startswith("lnL"):
                return float(line.split(":")[3].strip().split()[0])

        return None


    def getBranchNames(self, lines=None):
        """Get numbering of ancestral nodes"""

        lines = self.setupLines(lines)

        for line in lines:
            if line.startswith("(1) Branch lengths and substitution pattern"):
                break
        else:
            raise Exception("no branch names found")

        branches = lines.next().strip().split()
        names = [branch.split("..") for branch in branches]
        dists = map(float, lines.next().strip().split())

        # add nodes to tree
        tree = treelib.Tree()        
        for name in set(util.flatten(names)):                    
            tree.add(treelib.TreeNode(name))

        # link up nodes
        for (top, bot), dist in zip(names, dists):
            tree.add_child(tree.nodes[top], tree.nodes[bot])
            tree.nodes[bot].dist = dist

        # find root
        for node in tree:
            if node.parent is None:
                tree.root = node
                break

        return tree


    def getPampReconstruction(self, lines=None):
        """Get alignment of ancestral sequences"""
        lines = self.setupLines(lines)

        # find start of reconstruction
        for line in lines:
            if line.startswith("list of extant and reconstructed sequences"):
                lines.next()
                break
        else:
            raise Exception("no reconstruction found")
                      

        aln = fasta.FastaDict()

        # parse out sequences
        for line in lines:

            if line.startswith("node #"):
                # parse ancestral seq
                _, num, rest = line.split(" ", 2)
                name = num[1:]
            else:
                # parse extant seq
                name, rest = line.split(" ", 1)

            seq = rest.strip().replace(" ", "")
            aln[name] = seq
            
            # end of block
            if len(line) <= 1:
                break

        return aln


    def get_branch_site_params(self, lines=None):
        """Return the parameters found in the branch site model"""
        lines = self.setupLines(lines)
        
        # find start of params
        for line in lines:
            if line.startswith("dN/dS for site classes (K=4)"):
                lines.next() # skip blank line
                line = lines.next() 
                break
        else:
            raise Exception("no reconstruction found")

        tokens = line.rstrip().split()
        site_names = tokens[2:6]
        line = lines.next()

        tokens = line.rstrip().split()
        proportions = dict(zip(site_names,
                               map(float, tokens[1:5])))
        lines.next()

        # skip line
        line = lines.next()

        tokens = line.rstrip().split()
        omegas = dict(zip(site_names,
                          map(float, tokens[2:6])))
        
        return proportions, omegas
                      
            


def lrt_M1M2(m1lnl, m2lnl):

    """
    l0 = Model 1 = nearly neutral
    l1 = Model 2 = positively selected
    df = 2
    2dl = 2(l1 - l0)
    """
    import rpy

    dl = 2 * (m2lnl - m1lnl)
    pvalue = 1.0 - rpy.r.pchisq(dl, df=2)
    return pvalue


def lrt_branch_site(lnl_null, lnl_alt):
    """
    Do LRT for ModelA_null-ModelA

    The null distribution should be a 50:50 mixture of point mass 0
    and chi_1^2 , so that the critical values at the 5% and 1% levels
    are 2.71 and 5.41, respectively. We recommend use of chi_1^2 (with
    critical values 3.84 and 5.99) to guide against violations of
    model assumptions.
    """

    #assert lnl_alt >= lnl_null, (lnl_alt, lnl_null)
    
    import rpy

    dl = 2 * (lnl_alt - lnl_null)
    pvalue = rpy.r.pchisq(dl, df=1, lower_tail=False) / 2.0
    return pvalue

        

# Example codeml.ctl
"""
      seqfile = stewart.aa * sequence data filename
     treefile = stewart.trees      * tree structure file name
      outfile = mlc           * main result file name

        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
*        ndata = 10
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 2
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
        method = 0   * 0: simultaneous; 1: one branch at a time


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
"""
