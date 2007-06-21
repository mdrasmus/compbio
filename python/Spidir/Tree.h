//=============================================================================
// SPIDIR Tree datastructure

#ifndef SPIDIR_TREE_H
#define SPIDIR_TREE_H

#include <stdlib.h>
#include <string>

#include "common.h"

using namespace std;


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

    int name;
    Node *parent;
    Node **children;
    int nchildren;
    
    // thinking about whether to keep this here...
    float dist;
};


// A phylogenetic tree
class Tree
{
public:
    Tree(int nnodes) :
        nnodes(nnodes),
        root(NULL),
        nodes(NULL)
    {
        nodes = new Node* [nnodes];
        for (int i=0; i<nnodes; i++)
            nodes[i] = new Node();
    }
    
    virtual ~Tree()
    {
        if (nodes) {
            for (int i=0; i<nnodes; i++)
                delete nodes[i];
            delete [] nodes;
        }
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
    
    bool isRooted()
    {
        return (root != NULL && root->nchildren == 2);
    }
    
    Node *get(int name)
    {
        return nodes[name];
    }
    
    void reroot(Node *newroot, bool onBranch=true);
    
    int nnodes;
    Node *root;
    Node **nodes;
};


class SpeciesTree : public Tree
{
public:
    SpeciesTree(int nnodes) :
        Tree(nnodes)
    {
        depths = new int [nnodes];
    }
    
    
    void setDepths(Node *node=NULL, int depth=0)
    {
        if (node == NULL)
            node = root;
        
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
void writeNewick(Tree *tree, string *names, Node *node=NULL, int depth=0);


#endif // SPDIR_TREE_H
