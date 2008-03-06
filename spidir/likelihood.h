#ifndef SPIDIR_LIKELIHOOD_H
#define SPIDIR_LIKELIHOOD_H

#include "Tree.h"
#include "spidir.h"

namespace spidir {

float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate, 
             float predupprob=1.0, float dupprob=1.0, 
             float alpha=0, float beta=0, bool onlyduploss=false);

float treelk(Tree *tree,
             SpeciesTree *stree,
             int *recon, int *events, SpidirParams *params,
             float generate, 
             float predupprob, float dupprob, bool onlyduploss=false);

float rareEventsLikelihood(Tree *tree, SpeciesTree *stree, int *recon, 
                           int *events,
                           float predupprob, float dupprob);


float maxPosteriorGeneRate(Tree *tree, SpeciesTree *stree,
                           int *recon, int *events, SpidirParams *params);

float maxPosteriorGeneRate(int nnodes, int *ptree, float *dists,
                           int nsnodes, int *pstree, 
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta);

void generateBranchLengths(int nnodes, int *ptree, 
                           int nsnodes, int *pstree,
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta,
                           float *dists);

void generateBranchLengths(Tree *tree,
                           SpeciesTree *stree,
                           int *recon, int *events,
                           SpidirParams *params);


float birthDeathDensity(float *times, int ntimes, float maxtime, 
                        float birthRate, float deathRate);

float birthDeathTreePrior(Tree *tree, SpeciesTree *stree, int *recon, 
                          int *events, float birthRate, float deathRate);

} // namespace spidir

#endif // SPIDIR_LIKELIHOOD_H
