#ifndef SPIDIR_LIKELIHOOD_H
#define SPIDIR_LIKELIHOOD_H

#include "Tree.h"
#include "spidir.h"
#include "birthdeath.h"

namespace spidir {

typedef void (*geneRateCallback) (float generate, Tree *tree, void *userdata);

float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate, 
             float predupprob=1.0, float dupprob=1.0, float lossprob=1.0,
             float alpha=0, float beta=0, bool onlyduploss=false);

float treelk(Tree *tree,
             SpeciesTree *stree,
             int *recon, int *events, SpidirParams *params,
             float generate, 
             float predupprob, float dupprob, float lossprob,
             bool onlyduploss=false, bool oldduploss=false,
             bool duploss=true);


float rareEventsLikelihood(Tree *tree, SpeciesTree *stree, int *recon, 
                           int *events,
                           float predupprob, float dupprob, float lossprob);

float rareEventsLikelihood_old(Tree *tree, SpeciesTree *stree, int *recon, 
                               int *events,
                               float predupprob, float dupprob);

float maxPosteriorGeneRate(Tree *tree, SpeciesTree *stree,
                           int *recon, int *events, SpidirParams *params);

float maxPosteriorGeneRate(int nnodes, int *ptree, float *dists,
                           int nsnodes, int *pstree, 
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta);


void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs,
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *gene2species,
                             SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata);


void samplePosteriorGeneRate(Tree *tree,
                             int nseqs, char **seqs, 
                             const float *bgfreq, float ratio,
                             SpeciesTree *stree,
                             int *recon, int *events, SpidirParams *params,
                             int nsamples,
                             geneRateCallback callback,
                             void *userdata);


void generateBranchLengths(int nnodes, int *ptree, 
                           int nsnodes, int *pstree,
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta,
                           float *dists);

void generateBranchLengths(Tree *tree,
                           SpeciesTree *stree,
                           int *recon, int *events,
                           SpidirParams *params,
                           float generate=-1.0, 
                           int subnode=-1, int subchild=-1);


} // namespace spidir

#endif // SPIDIR_LIKELIHOOD_H
