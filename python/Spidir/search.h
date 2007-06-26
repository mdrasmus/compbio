#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

void proposeNni(Tree *tree, Node *node1, Node *node2, int change=0);

Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                int nseqs, int seqlen, char **seqs,
                int niter=500);

#endif // SPIDIR_SEARCH_H
