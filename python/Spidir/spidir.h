#ifndef SPIDIR_H
#define SPIDIR_H

/*=============================================================================

    SPIDIR - SPecies Informed DIstance-base Reconstruction
    Matt Rasmussen
    Copyright 2007
    
    Public interface to spidir C++ library.
    
=============================================================================*/

#include <stdlib.h>

// make functions linkable with C
extern "C" {


float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate, float disterror,
             float predupprob=1.0, float dupprob=1.0, float errorlogl=0,
             float alpha=0, float beta=0);

void parsimony(int nnodes, int *ptree, int nseqs, char **seqs, float *dists,
               bool buildAncestral=false, char **ancetralSeqs=NULL);

void findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                       float *dists, float *bgfreq, float ratio, int maxiter);

} // extern C


#endif // SPIDIR_H
