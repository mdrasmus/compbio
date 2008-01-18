#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H

#include "spidir.h"


namespace spidir {

class TopologyProposer
{
public:
    TopologyProposer() {}
    virtual ~TopologyProposer() {}
    virtual void propose(Tree *tree) {}
    virtual void revert(Tree *tree) {}
    virtual bool more() { return false; }
};


class NniProposer: public TopologyProposer
{
public:
    NniProposer(SpeciesTree *stree=NULL, int *gene2species=NULL, int niter=500);

    virtual void propose(Tree *tree);
    virtual void revert(Tree *tree);
    virtual bool more();

protected:    
    Node *nodea;
    Node *nodeb;
    Node *nodec;
    Node *noded;
    Node *oldroot1;
    Node *oldroot2;
    SpeciesTree *stree;
    int *gene2species;
    int niter;
    int iter;
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
              float *bgfreq, float tsvratio, int maxiter, bool useLogl=true);
    virtual float findLengths(Tree *tree);

    int nseqs;
    int seqlen;
    char **seqs;    
    float *bgfreq;
    float tsvratio;
    int maxiter;
    bool useLogl;
};


class BranchLikelihoodFunc
{
public:
    BranchLikelihoodFunc() {}
    virtual ~BranchLikelihoodFunc() {}
    
    virtual float likelihood(Tree *tree) { return 0.0; }
};


class SpidirBranchLikelihoodFunc : public BranchLikelihoodFunc
{
public:
    SpidirBranchLikelihoodFunc(int nnodes, SpeciesTree *stree, 
                               SpidirParams *params, 
                               int *gene2species,
                               float predupprob, float dupprob,
                               bool estGenerate);
    virtual float likelihood(Tree *tree);
    
protected:
    int nnodes;
    SpeciesTree *stree;
    SpidirParams *params;
    int *gene2species;
    ExtendArray<int> recon;
    ExtendArray<int> events;
    float predupprob;
    float dupprob;  
    bool estGenerate;
};



Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species);

Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                string *genes, int nseqs, int seqlen, char **seqs,
                int niter=500, float predupprob=0.01, float dupprob=1.0);


Tree *searchMCMC(Tree *initTree, 
                 string *genes, int nseqs, int seqlen, char **seqs,
                 BranchLikelihoodFunc *lkfunc,
                 TopologyProposer *proposer,
                 BranchLengthFitter *fitter);


} // namespace spidir

#endif // SPIDIR_SEARCH_H
