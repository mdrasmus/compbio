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
    
    
    // do neighbor joining
    int nnodes = aln->nseqs * 2 - 1;
    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    
    neighborjoin(aln->nseqs, distmat, ptree, dists);
    
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    parsimony(&tree, aln->nseqs, aln->seqs);
    
    writeNewick(&tree, aln->names.get());
    
    proposeNni(&tree, tree.nodes[0].parent, tree.nodes[0].parent->parent, 0);
    
    writeNewick(&tree, aln->names);
}
