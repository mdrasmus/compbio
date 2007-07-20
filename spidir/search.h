#ifndef SPIDIR_SEARCH_H
#define SPIDIR_SEARCH_H


class TopologyProposer
{
public:
    TopologyProposer()
    {}
    
    virtual ~TopologyProposer()
    {}

    virtual void propose(Tree *tree)
    {}
    
    virtual void revert(Tree *tree)
    {}
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

extern NniProposer nniProsposer;


void proposeNni(Tree *tree, Node *node1, Node *node2, int change=0);

Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                int nseqs, int seqlen, char **seqs,
                int niter=500, 
                TopologyProposer *proposer=&nniProsposer);

#endif // SPIDIR_SEARCH_H
