#ifndef SPIDIR_BRANCHLEN_H
#define SPIDIR_BRANCHLEN_H


#include "Tree.h"

namespace spidir {

float findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                             float *dists, float *bgfreq, float ratio, 
                             int maxiter, bool parsinit=false);

float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                            float *bgfreq, float ratio, int maxiter=100);

template <class Model>
float getTotalLikelihood(ExtendArray<float*> &lktable, Tree *tree, 
                         int nseqs, int seqlen, char **seqs, Model &model,
                         float *bgfreq);

void makeHkyMatrix(float *bgfreq, float ratio, float t, float *matrix);

} // namespace spidir

#endif // SPIDIR_BRANCHLEN_H
