#ifndef SPIDIR_LIKELIHOOD_H
#define SPIDIR_LIKELIHOOD_H

#include "Tree.h"
#include "spidir.h"

float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate, float disterror,
             float predupprob=1.0, float dupprob=1.0, float errorlogl=0,
             float alpha=0, float beta=0);

float treelk(Tree *tree,
             SpeciesTree *stree,
             int *recon, int *events, SpidirParams *params,
             float generate, float disterror,
             float predupprob, float dupprob, float errorlogl);



#endif // SPIDIR_LIKELIHOOD_H
