#ifndef SPIDIR_PARSIMONY_H
#define SPIDIR_PARSIMONY_H

#include "Tree.h"

namespace spidir {

void parsimony(int nnodes, int *ptree, int nseqs, char **seqs, float *dists,
               bool buildAncestral=false, char **ancetralSeqs=NULL);

void parsimony(Tree *tree, int nseqs, char **seqs,
               bool buildAncestral=false, char **ancetralSeqs=NULL);

} // namespace spidir

#endif
