/*=============================================================================

    Test SPIDIR functions

=============================================================================*/

#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>


#include "common.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "Tree.h"
#include "ConfigParam.h"


using namespace std;




int main(int argc, char **argv)
{
    // TODO: add default values
    
    // parameters
    string alignfile;    
    string smapfile;
    string streefile;
    string paramsfile;
    int niter = 0;
    
    
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
    
    config.add(new ConfigParamComment("Miscellaneous"));
    config.add(new ConfigParam<int>("-i", "--niter", "<# iterations>", 
                                    &niter, 100, "number of iterations"));
    
    
    if (argc < 2)
        config.printHelp();
    if (!config.parse(argc, (const char**) argv)) {
        return 1;
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
    
    
    // calc distmatrix for neighbor joining
    Matrix<float> distmat(aln->nseqs, aln->nseqs);
    calcDistMatrix(aln->nseqs, aln->seqlen, aln->seqs, distmat);
    
    
    // do neighbor joining and parsimony distance
    int nnodes = aln->nseqs * 2 - 1;
    /*
    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    
    neighborjoin(aln->nseqs, distmat, ptree, dists);
    Tree tree(nnodes);    
    ptree2tree(nnodes, ptree, &tree);
    parsimony(&tree, aln->nseqs, aln->seqs);
    */

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
