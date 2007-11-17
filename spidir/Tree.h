//=============================================================================
// SPIDIR Tree datastructure

#ifndef SPIDIR_TREE_H
#define SPIDIR_TREE_H

#include <stdlib.h>
#include <string>

#include "common.h"

using namespace std;


namespace spidir {



// A node in the phylogenetic tree
class Node
{
public:
    Node(int nchildren=0) :
        name(-1),
        parent(NULL),
        children(NULL),
        nchildren(nchildren),
        dist(0.0)
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

    void hashkey(int *key);
    
    
    void reroot(Node *newroot, bool onBranch=true);
    void reroot(Node *node1, Node *node2);
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



struct HashTopology {
    static unsigned int hash(const ExtendArray<int> &key)
    {
        unsigned int h = 0, g;
        
        for (int i=0; i<key.size(); i++) {
            h = (h << 4) + key[i];
            if ((g = h & 0xF0000000))
                h ^= g >> 24;
            h &= ~g;
        }
        
        return h;
    }    
};


void getTreeSortedPostOrder(Tree *tree, ExtendArray<Node*> *nodes, 
                      int *ordering, Node *node=NULL);
void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);
void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);


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
// Input/output
void printFtree(int nnodes, int **ftree);
void printTree(Tree *tree, Node *node=NULL, int depth=0);
float readDist(FILE *infile, int &depth);
char readChar(FILE *stream, int &depth);
char readUntil(FILE *stream, string &token, const char *stops, int &depth);


} // namespace spidir



/*
DEPRECATED: Not needed anymore

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

*/




#endif // SPDIR_TREE_H

