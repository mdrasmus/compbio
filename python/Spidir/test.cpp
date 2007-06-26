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
    
    ConfigParser config;
    
    string smap;
    int niter = 0;
    
    config.add(new ConfigParam<string>("-S", "--smap", "<species map>", 
                                       &smap, "gene to species map"));
    
    config.add(new ConfigParamComment("Miscellaneous"));
    config.add(new ConfigParam<int>("-i", "--niter", "<# iterations>", 
                                       &niter, "number of iterations"));
    
    
    
    if (argc < 2)
        config.printHelp();
    if (!config.parse(argc, (const char**) argv)) {
        return 1;
    }
    

/*
    // read sequences
    Sequences *aln = readAlignFasta(argv[1]);
    writeFasta("out.fa", aln);
    assert(checkSequences(aln->nseqs, aln->seqlen, aln->seqs));
    
    // calc distmatrix
    Matrix<float> distmat(aln->nseqs, aln->nseqs);
    calcDistMatrix(aln->nseqs, aln->seqlen, aln->seqs, 
                   distmat);
    
    // write dist matrix
    if (argc > 2)
        writeDistMatrix(argv[2], aln->nseqs, distmat, aln->names);
    
    
    //SpidirParams params = SpidirParams(nsnodes, mu, sigma, alpha, beta);
    
    // do neighbor joining
    int nnodes = aln->nseqs * 2 - 1;
    Tree tree(nnodes);    
    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    
    neighborjoin(aln->nseqs, distmat, ptree, dists);
    ptree2tree(nnodes, ptree, &tree);
    parsimony(&tree, aln->nseqs, aln->seqs);
    
    
    for (int i=0; i<100; i++) {
        Node *node = tree.root;
        while (node->parent == NULL || node->name < aln->nseqs)
            node = tree.nodes[int(rand() / float(RAND_MAX) * tree.nnodes)];
        proposeNni(&tree, node, node->parent, int(rand() / float(RAND_MAX) * 2));
        
        parsimony(&tree, aln->nseqs, aln->seqs);
        //printTree(&tree);
    }
*/
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
