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


int main(int argc, char **argv)
{

    // read sequences
    Sequences *seqs = readFasta(argv[1]);
    seqs->setAlignLength();
    assert(checkSequences(seqs->nseqs, seqs->seqlen, seqs->seqs.getArray()));
    
    // calc distmatrix
    Matrix<float> distmat(seqs->nseqs, seqs->nseqs);
    calcDistMatrix(seqs->nseqs, seqs->seqlen, seqs->seqs.getArray(), 
                   distmat.getMatrix());
    
    // write dist matrix
    if (argc > 2)
        writeDistMatrix(argv[2], seqs->nseqs, distmat.getMatrix(), 
                        seqs->names.getArray());
    
    
    // do neighbor joining
    int nnodes = seqs->nseqs * 2 - 1;
    int *ptree = new int [nnodes];
    float *dists = new float [nnodes];
    
    neighborjoin(seqs->nseqs, distmat.getMatrix(), ptree, dists);
    
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    parsimony(&tree, seqs->nseqs, seqs->seqs.getArray());
    
    writeNewick(&tree, seqs->names.getArray());
    
    proposeNni(&tree, tree.nodes[0].parent, tree.nodes[0].parent->parent, 0);

    writeNewick(&tree, seqs->names.getArray());
}
