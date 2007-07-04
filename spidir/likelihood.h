#ifndef SPIDIR_LIKELIHOOD_H
#define SPIDIR_LIKELIHOOD_H



float treelk(Tree *tree,
             SpeciesTree *stree,
             int *recon, int *events, SpidirParams *params,
             float generate, float disterror,
             float predupprob, float dupprob, float errorlogl);



#endif // SPIDIR_LIKELIHOOD_H
