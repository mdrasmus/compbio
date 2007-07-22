//=============================================================================
// SPIDIR Tree datastructure

#ifndef SPIDIR_TREE_H
#define SPIDIR_TREE_H

#include <stdlib.h>
#include <string>

#include "common.h"

using namespace std;


namespace spidir {

// events
enum {
    EVENT_GENE = 0,
    EVENT_SPEC = 1,
    EVENT_DUP = 2
};



// A node in the phylogenetic tree
class Node
{
public:
    Node(int nchildren=0) :
        name(-1),
        parent(NULL),
        children(NULL),
        nchildren(nchildren)
    {
        if (nchildren != 0)
            setChildren(nchildren);
    }
    
    ~Node()
    {
        if (children)
            delete [] children;
    }
    
    void setChildren(int _nchildren)
    {
        children = resize(children, nchildren, _nchildren);
        nchildren = _nchildren;
    }
    
    void allocChildren(int _nchildren)
    {
        children = new Node* [_nchildren];
    }
    
    bool isLeaf() {
        return nchildren == 0;
    }
    
    void addChild(Node *node)
    {
        setChildren(nchildren + 1);
        children[nchildren - 1] = node;
        node->parent = this;
    }

    int name;
    Node *parent;
    Node **children;
    int nchildren;
    float dist;
    string leafname;
};


// A phylogenetic tree
class Tree
{
public:
    Tree(int nnodes=0) :
        nnodes(nnodes),
        root(NULL),
        nodes(nnodes, 100)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i] = new Node();
    }
    
    virtual ~Tree()
    {
        for (int i=0; i<nnodes; i++)
            delete nodes[i];
    }
    
    void setDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i]->dist = dists[i];
    }
    
    void getDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            dists[i] = nodes[i]->dist;
    }
    
    void setLeafNames(string *names)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i]->leafname = names[i];
    }
    
    void getLeafNames(string *names)
    {
        for (int i=0; i<nnodes; i++)
            names[i] = nodes[i]->leafname;
    }
    
    void getNames(string *names)
    {
        for (int i=0; i<nnodes; i++) {
            if (nodes[i]->isLeaf())
                names[i] = nodes[i]->leafname;
            else {
                char numstr[21];
                snprintf(numstr, 20, "%d", nodes[i]->name);
                names[i] = numstr;
            }
        }
    }
    
    bool isRooted()
    {
        return (root != NULL && root->nchildren == 2);
    }
    
    Node *getNode(int name)
    {
        return nodes[name];
    }
    
    Node *addNode(Node *node)
    {
        nodes.append(node);
        node->name = nodes.size() - 1;
        nnodes = nodes.size();
        return node;
    }
    
    void reroot(Node *newroot, bool onBranch=true);
    Tree *copy();
    Node *readNode(FILE *infile, Node *parent, int &depth);
    bool readNewick(FILE *infile);
    bool readNewick(const char *filename);
    void writeNewick(FILE *out=stdout, Node *node=NULL, int depth=0);
    bool writeNewick(const char *filename);
    bool assertTree();
    
    int nnodes;
    Node *root;
    ExtendArray<Node*> nodes;
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



// right now works for pre-order traversal
class TreeWalker
{
public:
    TreeWalker(int nnodes, int **ftree, int start) :
        nnodes(nnodes),
        node(start),        
        ftree(ftree)
    {
        stack = new int [nnodes];
        stacki = 0;
        stack[0] = start;
    }
    
    ~TreeWalker()
    {
        delete [] stack;
    }
    
    bool recurse(int node)
    {
        // descend tree
        if (ftree[node][0] != -1) {
            stack[++stacki] = ftree[node][1];
            stack[++stacki] = ftree[node][0];
            return true;
        } else
            return false;
    }
    
    int next()
    {
        // done walking
        if (stacki < 0)
            return -1;
        
        // get next node
        node = stack[stacki];
        
        // pop node of stack
        stacki--;
        
        return node;
    }
    
    
    int nnodes;
    int node;
    int **ftree;
    
    int *stack;
    int stacki;
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



void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);
void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);
float readDist(FILE *infile, int &depth);


//=============================================================================
// visualization

void displayTree(Tree *tree, FILE *outfile=stdout, 
                 float xscale=20.0, int yscale=2);

//=============================================================================
// conversion functions

// creates a forward tree from a parent tree
void makeFtree(int nnodes, int *ptree, int ***ftree);
void freeFtree(int nnodes, int **ftree);
// create a tree object from a parent tree array
void ptree2tree(int nnodes, int *ptree, Tree *tree);
void tree2ptree(Tree *tree, int *ptree);

//=============================================================================
// reconciliation functions
void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches);
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon);
void labelEvents(Tree *tree, int *recon, int *events);

//=============================================================================
// Input/output
void printFtree(int nnodes, int **ftree);
void printTree(Tree *tree, Node *node=NULL, int depth=0);


} // namespace spidir

#endif // SPDIR_TREE_H
