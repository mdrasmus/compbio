/*=============================================================================

    Test SPIDIR functions

=============================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "common.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "Tree.h"


using namespace std;

float maxCubicRoot(float a, float b, float c);

int main(int argc, char **argv)
{
    Tree tree;
    SpeciesTree stree;
    
    tree.readNewick(argv[1]);
    stree.readNewick(argv[2]);
    stree.setDepths();
    
    stree.writeNewick();
    
    ExtendArray<string> names(tree.nnodes);
    tree.getLeafNames(names);
    
    ExtendArray<string> snames(stree.nnodes);
    stree.getLeafNames(snames);
    
    /*
    for (int i=0; i<tree.nnodes; i++) {
        printf("name %d = '%s'\n", i, names[i].c_str());
    }
    
    for (int i=0; i<stree.nnodes; i++) {
        printf("sname %d = '%s'\n", i, snames[i].c_str());
    }*/
    
    int gene2species[17] = {5,   0,  11, 10,  9,  7,  6, 3, 4, 
                            -1, -1, -1, -1, -1, -1, -1, -1};
    
    
    
    ExtendArray<int> recon(tree.nnodes);
    ExtendArray<int> events(tree.nnodes);
    reconcile(&tree, &stree, gene2species, recon);
    labelEvents(&tree, recon, events);
    
    printTree(&tree);
    printTree(&stree);
    
    printIntArray(recon, recon.size());
    printIntArray(events, events.size());
    
    return 0;

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
