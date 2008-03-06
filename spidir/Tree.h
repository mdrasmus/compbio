//=============================================================================
// SPIDIR Tree datastructure

#ifndef SPIDIR_TREE_H
#define SPIDIR_TREE_H

#include <stdlib.h>
#include <string>

#include "common.h"

using namespace std;


/*

    Parent Array Format (ptree)
    

           4
          / \
         3   \
        / \   \
       /   \   \
       0   1    2

     Then the parent tree array representation is
       ptree = [3, 3, 4, 4, -1]
     such that ptree[node's id] = node's parent's id

     In addition, the following must be true
     1. tree must be binary: n leaves, n-1 internal nodes
     2. leaves must be numbered 0 to n-1
     3. internal nodes are numbered n to 2n-2
     4. root must be numbered 2n-2
     5. the parent of root is -1
     6. the length of ptree is 2n-1

*/



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
            allocChildren(nchildren);
    }
    
    ~Node()
    {
        if (children)
            delete [] children;
    }
    
    // Sets and allocates the number of children '_nchildren'
    void setChildren(int _nchildren)
    {
        children = resize(children, nchildren, _nchildren);
        nchildren = _nchildren;
    }
    
    // Allocates the number of children '_nchildren'
    void allocChildren(int _nchildren)
    {
        children = new Node* [_nchildren];
    }
    
    // Returns whether the node is a leaf
    bool isLeaf() {
        return nchildren == 0;
    }
    
    // Adds a node 'node' to be a child
    void addChild(Node *node)
    {
        setChildren(nchildren + 1);
        children[nchildren - 1] = node;
        node->parent = this;
    }
    
    int name;           // node name id (matches index in tree.nodes)
    Node *parent;       // parent pointer
    Node **children;    // array of child pointers (size = nchildren)
    int nchildren;      // number of children
    float dist;         // branch length above node
    string longname;    // node name (used mainly for leaves only)
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
    
    // Sets the branch lengths of the tree
    //  Arguments:
    //      dists: array of lengths (size = nnodes)
    void setDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            nodes[i]->dist = dists[i];
    }
    
    
    // Gets the branch lengths of the tree
    //  Arguments:
    //      dists: output array (size = nnodes) for storing branch lengths
    void getDists(float *dists)
    {
        for (int i=0; i<nnodes; i++)
            dists[i] = nodes[i]->dist;
    }
    
    // Sets the leaf names of the tree
    //  Arguments:
    //      names:      array of names (size > # leaves)
    //      leavesOnly: whether to only set names for leaves, or all nodes
    void setLeafNames(string *names, bool leavesOnly=true)
    {
        for (int i=0; i<nnodes; i++) {
            if (leavesOnly && !nodes[i]->isLeaf())
                nodes[i]->longname = "";
            else
                nodes[i]->longname = names[i];
        }
    }
    
    void reorderLeaves(string *names);
    
    // Gets leaf names of the nodes of a tree
    // Internal nodes are often named "" (empty string)
    //  Arguments:
    //      names:      output array for storing node names (size > # leaves)
    //      leavesOnly: whether to only get names for leaves, or all nodes
    void getLeafNames(string *names, bool leavesOnly=true)
    {
        for (int i=0; i<nnodes; i++)
            if (!leavesOnly || nodes[i]->isLeaf())
                names[i] = nodes[i]->longname;
    }
    
    
    // Gets names of the nodes of a tree
    // This differs from getLeafNames in that internal nodes will be named 
    // after their name id (int) converted to a string.
    //  Arguments:
    //      names: output array (size = nnodes) for node names
    void getNames(string *names)
    {
        for (int i=0; i<nnodes; i++) {
            if (nodes[i]->isLeaf())
                names[i] = nodes[i]->longname;
            else {
                char numstr[21];
                snprintf(numstr, 20, "%d", nodes[i]->name);
                names[i] = numstr;
            }
        }
    }
    
    // Returns whether tree is rooted
    bool isRooted()
    {
        return (root != NULL && root->nchildren == 2);
    }
    
    // Returns the pointer to a node given is name id 'name'
    Node *getNode(int name)
    {
        return nodes[name];
    }
    
    // Adds a node 'node' to the tree
    // This will also set the node's name id
    Node *addNode(Node *node)
    {
        nodes.append(node);
        node->name = nodes.size() - 1;
        nnodes = nodes.size();
        return node;
    }
    
    // Compute a topology hash of the tree
    //  Arguments:
    //      key: output array (size = nnodes) containing a unique sequence of
    //           integers for the tree
    void hashkey(int *key);
    
    bool sameTopology(Tree *other);
    
    // Roots the tree on branch 'newroot'
    void reroot(Node *newroot, bool onBranch=true);
    
    // Roots the tree on branch connecting 'node1' and 'node2'
    void reroot(Node *node1, Node *node2);
    
    // Returns a new copy of the tree
    Tree *copy();
    
    // Reads a tree structure from an input stream 'infile'
    // Returns true on success
    bool readNewick(FILE *infile);
    
    // Reads a tree structure from a file 'filename'
    // Returns true on success
    bool readNewick(const char *filename);
    
    // Writes a tree structure to an output stream 'out'
    void writeNewick(FILE *out=stdout, Node *node=NULL, int depth=0, bool oneline=false);
    
    // Writes a tree structure to a file 'filename'
    bool writeNewick(const char *filename, bool oneline=false);
    
    // Returns whether the tree is self consistent
    bool assertTree();
    
protected:
    // Reads a single node from an open file
    Node *readNode(FILE *infile, Node *parent, int &depth);    
    
public:    
    int nnodes;                 // number of nodes in tree
    Node *root;                 // root of the tree (NULL if no nodes)
    ExtendArray<Node*> nodes;   // array of nodes (size = nnodes)
};


// A hash function for a topology key to an integer
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


//=============================================================================
// Tree traversals

void getTreeSortedPostOrder(Tree *tree, ExtendArray<Node*> *nodes, 
                            int *ordering, Node *node=NULL);
void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);
void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL);


//=============================================================================
// visualization

void displayTree(Tree *tree, FILE *outfile=stdout, 
                 float xscale=20.0, int yscale=2);
void displayTreeMatrix(Tree *tree, float xscale, int yscale, 
                       char ***matrix, int *nrows, int *ncols);


//=============================================================================
// conversion functions

// Creates a 'forward tree' from a 'parent tree'
void makeFtree(int nnodes, int *ptree, int ***ftree);

// Deallocates a 'forward tree'
void freeFtree(int nnodes, int **ftree);

// Creates a tree object from a 'parent tree' array
void ptree2tree(int nnodes, int *ptree, Tree *tree);

// Creates a 'parent tree' from a tree object
void tree2ptree(Tree *tree, int *ptree);


//=============================================================================
// Input/output

void printFtree(int nnodes, int **ftree);
void printTree(Tree *tree, Node *node=NULL, int depth=0);


} // namespace spidir


#endif // SPDIR_TREE_H

