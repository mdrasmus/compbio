#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

#include "spidir.h"



class TopologyProposer
{
public:
    TopologyProposer() {}
    virtual ~TopologyProposer() {}
    virtual void propose(Tree *tree) {}
    virtual void revert(Tree *tree) {}
};


class NniProposer: public TopologyProposer
{
public:
    NniProposer();

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    
    Node *node1;
    Node *node2;
    Node *node3;
    Node *node4;
    Node *oldroot;
    int change1;
    int change2;
};

class BranchLengthFitter
{
public:
    BranchLengthFitter() {}
    virtual ~BranchLengthFitter() {}
    virtual float findLengths(Tree *tree) {return 0.0;}
};


class ParsimonyFitter : public BranchLengthFitter
{
public:
    ParsimonyFitter(int nseqs, int seqlen, char **seqs);
    virtual float findLengths(Tree *tree);
    
    int nseqs;
    int seqlen;
    char **seqs;
};


class HkyFitter : public BranchLengthFitter
{
public:
    HkyFitter(int nseqs, int seqlen, char **seqs, 
              float *bgfreq, float tsvratio, int maxiter);
    virtual float findLengths(Tree *tree);

    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio;
    int maxiter;
    float logl;
};


extern NniProposer nniProposer;


void proposeNni(Tree *tree, Node *node1, Node *node2, int change=0);

Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                int nseqs, int seqlen, char **seqs,
                int niter=500, 
                TopologyProposer *proposer=&nniProposer,
                BranchLengthFitter *fitter=NULL
                );

#endif // SPIDIR_SEARCH_H
