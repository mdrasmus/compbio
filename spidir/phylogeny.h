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
// reconciliation functions
void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches);
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon);
void labelEvents(Tree *tree, int *recon, int *events);


} // namespace spidir

#endif // SPIDIR_PHYLOGENY_H
