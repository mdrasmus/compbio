#ifndef SPIDIR_BRANCHLEN_H
#define SPIDIR_BRANCHLEN_H


#include "Tree.h"

float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                            float *bgfreq, float ratio, int maxiter=100);

template <class Model>
float getTotalLikelihood(ExtendArray<float*> &lktable, Tree *tree, 
                         int nseqs, int seqlen, char **seqs, Model &model,
                         float *bgfreq);


#endif // SPIDIR_BRANCHLEN_H
