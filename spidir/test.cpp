/*=============================================================================

    Test SPIDIR functions

=============================================================================*/

#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>


#include "branchlen.h"
#include "common.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "Tree.h"
#include "ConfigParam.h"


using namespace std;





int test_reconstruct(int argc, char **argv)
{    
    // parameters
    string alignfile;    
    string smapfile;
    string streefile;
    string paramsfile;
    int niter = 0;
    bool help = false;
    
    
    // parse arguments
    ConfigParser config;
    config.add(new ConfigParam<string>(
        "-a", "--align", "<alignment fasta>", &alignfile, 
        "sequence alignment in fasta format"));
    config.add(new ConfigParam<string>(
        "-S", "--smap", "<species map>", &smapfile, 
        "gene to species map"));
    config.add(new ConfigParam<string>(
        "-s", "--stree", "<species tree>", &streefile, 
        "species tree file in newick format"));
    config.add(new ConfigParam<string>(
        "-p", "--param", "<spidir params file>", &paramsfile, 
        "SPIDIR branch length parameters file"));
    
    
    config.add(new ConfigParamComment("Miscellaneous options"));
    config.add(new ConfigParam<int>(
        "-i", "--niter", "<# iterations>", &niter, 100, 
        "number of iterations"));
    config.add(new ConfigSwitch(
        "-h", "--help", &help, "display help information"));

    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }
    
    if (help) {
        config.printHelp();
        return 0;
    }

    
    //============================================================
    // read species tree
    SpeciesTree stree;
    stree.readNewick(streefile.c_str());
    stree.setDepths();
    
    // read gene2species map
    Gene2species g;
    g.read(smapfile.c_str());
    

    // read sequences
    Sequences *aln;
    
    if ((aln = readAlignFasta(alignfile.c_str())) == NULL ||
        !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
        return 1;
    }

    // read SPIDIR parameters
    SpidirParams *params;
    if ((params = readSpidirParams(paramsfile.c_str())) == NULL)
    {
        printError("bad parameters file");
        return 1;
    }
    
    if (!params->order(&stree)) {
        printError("parameters do not correspond to the given species tree");
        return 1;
    }
    
    
    int nnodes = aln->nseqs * 2 - 1;

    // produce mapping array
    ExtendArray<string> genes(0, nnodes);
    ExtendArray<string> species(stree.nnodes);
    ExtendArray<int> gene2species(nnodes);
    genes.extend(aln->names, aln->nseqs);
    for (int i=aln->nseqs; i<nnodes; i++)
        genes.append("");
    stree.getLeafNames(species);   
    
    g.getMap(genes, nnodes, species, stree.nnodes, gene2species);
    
    // search
    Tree *toptree = searchMCMC(NULL, &stree,
                               params, gene2species,
                               aln->nseqs, aln->seqlen, aln->seqs,
                               niter);
    
    toptree->setLeafNames(genes);
    toptree->writeNewick();
    
    delete toptree;
    
    return 0;
}


int test_mledist(int argc, char **argv)
{
    // parameters
    string alignfile;
    
    // parse arguments
    ConfigParser config;
    config.add(new ConfigParam<string>(
        "-a", "--align", "<alignment fasta>", &alignfile, 
        "sequence alignment in fasta format"));
    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }
    
    
    // read sequences
    Sequences *aln;
    
    if ((aln = readAlignFasta(alignfile.c_str())) == NULL ||
        !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
        return 1;
    }
    
    int nnodes = aln->nseqs * 2 - 1;
    int nseqs = aln->nseqs;
    int seqlen = aln->seqlen;
    char **seqs = aln->seqs;
    
    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    Matrix<float> distmat(nseqs, nseqs);

    calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
    neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setLeafNames(aln->names);
    
    parsimony(&tree, nseqs, seqs);
    tree.writeNewick("before.tree");
    
    float bgfreq[4] = {.25, .25, .25, .25};
    float ratio = .5;
    int maxiter = 20;
    
    findMLBranchLengthsHky(&tree, nseqs, seqs, bgfreq, ratio, maxiter);
    tree.writeNewick("after.tree");
    
    return 0;
}


int test_gene2species(int argc, char **argv)
{
    Gene2species g;
    Tree tree, stree;
    
    tree.readNewick(argv[1]);
    stree.readNewick(argv[2]);
    g.read(argv[3]);
    
    ExtendArray<string> genes(tree.nnodes);
    tree.getLeafNames(genes);

    ExtendArray<string> species(stree.nnodes);
    stree.getLeafNames(species);

    
    for (int i=0; i<tree.nnodes; i++) {
        printf("'%s' -> '%s'\n", genes[i].c_str(), 
               g.getSpecies(genes[i]).c_str());
    }
    
    ExtendArray<int> map(tree.nnodes);
    g.getMap(genes, tree.nnodes, species, stree.nnodes, map);
    
    printIntArray(map, tree.nnodes);
    
    return 0;
}




int test_reroot(int argc, char **argv)
{
    Tree tree;
    
    tree.readNewick(argv[1]);
    
    char filename[100];
    
    for (int i=0; i<tree.nnodes; i++) {
        tree.reroot(tree.nodes[i]);
        snprintf(filename, 100, "%d.reroot.tree", i);
        tree.writeNewick(filename);
    }
}


int main(int argc, char **argv)
{

    if (argc < 2) {
        printf("choose a test:\n"
               "  reconstruct\n"
               "  gene2species\n"
               "  mledist\n"
               "  reroot\n");
        return 1;
    }

    string testname = argv[1];
    
    if (testname == "reconstruct") {
        test_reconstruct(argc-1, &argv[1]);
        
    } else if (testname == "gene2species") {
        test_gene2species(argc-1, &argv[1]);
        
    } else if (testname == "mledist") {
        test_mledist(argc-1, &argv[1]);
    } else if (testname == "reroot") {
        test_reroot(argc-1, &argv[1]);
    }
    
    return 0;
}
