#ifndef SPIDIR_COMMON_H
#define SPIDIR_COMMON_H

/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "spidir.h"



#define PI 3.1415926

// events
enum {
    EVENT_GENE = 0,
    EVENT_SPEC = 1,
    EVENT_DUP = 2
};

// fractional branches
enum {
    FRAC_NONE,
    FRAC_DIFF,
    FRAC_PARENT,
    FRAC_NODE
};



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
        nchildren = _nchildren;
        children = new Node* [nchildren];
    }
    
    void allocChildren(int _nchildren)
    {
        children = new Node* [_nchildren];
    }

    int name;
    Node *parent;
    Node **children;
    int nchildren;
    
    // thinking about whether to keep this here...
    float dist;
};


class Tree
{
public:
    Tree(int nnodes) :
        nnodes(nnodes),
        root(NULL),
        nodes(NULL)
    {
        nodes = new Node [nnodes];
    }
    
    virtual ~Tree()
    {
        if (nodes)
            delete [] nodes;
    }
    
    void setDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i].dist = dists[i];
    }
    
    int nnodes;
    Node *root;
    Node *nodes;
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




void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon);
void labelEvents(Tree *tree, int *recon, int *events);
void parsimony(int nnodes, int *ptree, int nseqs, char **seqs, float *dists);


//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s);
double gammln(double xx);
float gammalog(float x, float a, float b);

//=============================================================================
// data structure manipulation

// creates a forward tree from a parent tree
void makeFtree(int nnodes, int *ptree, int ***ftree);
void freeFtree(int nnodes, int **ftree);
// create a tree object from a parent tree array
void ptree2tree(int nnodes, int *ptree, Tree *tree);
void tree2ptree(Tree *tree, int *ptree);

//=============================================================================
// debug

void printIntArray(int *array, int size);
void printFloatArray(float *array, int size);
void printFtree(int nnodes, int **ftree);
void printTree(Tree *tree, Node *node=NULL, int depth=0);
void writeNewick(Tree *tree, Node *node=NULL, int depth=0);

#endif // SPIDIR_COMMON_H
