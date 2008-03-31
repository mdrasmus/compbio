#ifndef SPIDIR_BRANCHLEN_H
#define SPIDIR_BRANCHLEN_H


#include "Tree.h"

namespace spidir {

float calcHkySeqProb(Tree *tree, int nseqs, char **seqs, 
                     const float *bgfreq, float ratio);


float findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                             float *dists, const float *bgfreq, float ratio, 
                             int maxiter, bool parsinit=false);

float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                            const float *bgfreq, float ratio, 
                            int maxiter=100, int samples=0);

template <class Model>
float getTotalLikelihood(ExtendArray<float*> &lktable, Tree *tree, 
                         int nseqs, int seqlen, char **seqs, Model &model,
                         const float *bgfreq);

void makeHkyMatrix(const float *bgfreq, float ratio, float t, float *matrix);

} // namespace spidir

#endif // SPIDIR_BRANCHLEN_H
