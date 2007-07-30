#ifndef SPIDIR_PHYLOGENY_H
#define SPIDIR_PHYLOGENY_H

#include "Tree.h"


namespace spidir {

// events
enum {
    EVENT_GENE = 0,
    EVENT_SPEC = 1,
    EVENT_DUP = 2
};

class SpeciesTree : public Tree
{
public:
    SpeciesTree(int nnodes=0) :
        Tree(nnodes),
        depths(NULL)
    {
    }
    
    virtual ~SpeciesTree()
    {
        delete [] depths;
    }
    
    
    void setDepths(Node *node=NULL, int depth=0)
    {
        if (node == NULL)
            node = root;
        
        if (depths == NULL)
            depths = new int [nnodes];
        
        depths[node->name] = depth;
        
        for (int i=0; i<node->nchildren; i++)
            setDepths(node->children[i], depth+1);
    }
    
    int *depths;
};




class Gene2speciesRule
{
public:
    Gene2speciesRule(int rule=PREFIX, string expr="", string species="") :
        rule(rule),
        expr(expr),
        species(species)
    {
    }
    
    enum {
        PREFIX,
        SUFFIX,
        EXACT
    };
    
    int rule;
    string expr;
    string species;
};


class Gene2species
{
public:
    Gene2species() :
        m_rules(0, 20)
    {}
    
    const static string NULL_SPECIES;
    
    bool read(const char *filename);
    string getSpecies(string gene);
    bool getMap(string *genes, int ngenes, string *species, int nspecies, 
                int *map);
    
protected:
    ExtendArray<Gene2speciesRule> m_rules;
    
};


//=============================================================================
// phylogenetic reconstruction

void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches);


//=============================================================================
// reconciliation functions

void reconRoot(Tree *tree, SpeciesTree *stree, int *gene2species);
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon);
void labelEvents(Tree *tree, int *recon, int *events);
Node *treeLca(SpeciesTree *stree, Node *node1, Node *node2);

inline int countDuplications(int nnodes, int *events)
{
    int dups = 0;
    for (int i=0; i<nnodes; i++)
        if (events[i] == EVENT_DUP)
            dups++;
    return dups;
}


inline int reconcileNode(Node *node, SpeciesTree *stree, int *recon) {
    Node *snode1 = stree->nodes[recon[node->children[0]->name]];
    Node *snode2 = stree->nodes[recon[node->children[1]->name]];
    Node *snode = treeLca(stree, snode1, snode2);
    return snode->name;
}


inline int labelEventsNode(Node *node, int *recon)
{
    if (!node->isLeaf()) {
        for (int i=0; i<node->nchildren; i++)
            if (recon[node->name] == recon[node->children[i]->name])
                return EVENT_DUP;
        return EVENT_SPEC;
    } else {
        return EVENT_GENE;
    }
}



} // namespace spidir

#endif // SPIDIR_PHYLOGENY_H
