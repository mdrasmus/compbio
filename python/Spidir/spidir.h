#ifndef SPIDIR_H
#define SPIDIR_H

// publicly visible interface




// function prototypes

float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate, float disterror,
             float predupprob=1.0, float dupprob=1.0, float errorlogl=0,
             float alpha=0, float beta=0);

void parsimony(int nnodes, int *ptree, int nseqs, char **seqs, float *dists);

#endif // SPIDIR_H
