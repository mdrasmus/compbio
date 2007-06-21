//=============================================================================
// SPIDIR Tree datastructure

#include <assert.h>
#include <stdio.h>


#include "Tree.h"
#include "Matrix.h"

#define MAX_FLOAT 1e10


//=============================================================================
// phylogeny functions


// Neighbor-joining algorithm
void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches)
{
    Matrix<float> dists(ngenes*2-1, ngenes*2-1);
    float *restdists = new float [ngenes*2-1];
    int *leaves = new int [ngenes];
    int nleaves = ngenes;
    int newnode = ngenes;
    
    // initialize distances
    for (int i=0; i<ngenes; i++) {
        float r = 0.0;
        for (int j=0; j<ngenes; j++) {
            dists[i][j] = distmat[i][j];
            r += distmat[i][j];
        }
        restdists[i] = r / (ngenes - 2);
    }
    
    // initialize leaves
    for (int i=0; i<ngenes; i++)
        leaves[i] = i;
    
    
    // join loop
    while (nleaves > 2) {
        // search for closest genes
        float low = MAX_FLOAT;
        int lowi = -1, lowj = -1;
        
        for (int i=0; i<nleaves; i++) {
            for (int j=i+1; j<nleaves; j++) {
                int gene1 = leaves[i];
                int gene2 = leaves[j];
                float dist = dists[gene1][gene2] - restdists[gene1] 
                                                 - restdists[gene2];
                if (dist < low) {
                    low = dist;
                    lowi = i;
                    lowj = j;
                }
            }
        }
        
        // join gene1 and gene2
        int lowgene1 = leaves[lowi];
        int lowgene2 = leaves[lowj];
        int parent = newnode++;
        ptree[lowgene1] = parent;
        ptree[lowgene2] = parent;
        
        // set distances
        branches[lowgene1] = (dists[lowgene1][lowgene2] + 
                              restdists[lowgene1] - 
                              restdists[lowgene2]) / 2.0;
        branches[lowgene2] = dists[lowgene1][lowgene2] - branches[lowgene1];
        
        // gene1 and gene2 are no longer leaves, remove them from leaf set
        leaves[lowi] = parent;
        leaves[lowj] = leaves[nleaves-1];
        nleaves--;
        
        float r = 0;
        for (int i=0; i<nleaves; i++) {
            int gene = leaves[i];
            if (gene != parent) {
                float v = (dists[lowgene1][gene] + 
                           dists[lowgene2][gene] -
                           dists[lowgene1][lowgene2]) / 2.0;
                dists[parent][gene] = v;
                dists[gene][parent] = v;
                r += v;
            }
        }
        
        if (nleaves > 2)
            restdists[parent] = r / (nleaves - 2);
    }
    
    // join the last two genes, split the remaining dist evenly
    int gene1 = leaves[0];
    int gene2 = leaves[1];
    int parent = newnode++;
    
    ptree[gene1] = parent;
    ptree[gene2] = parent;
    ptree[parent] = -1;
    branches[gene1] = dists[gene1][gene2] / 2.0;
    branches[gene2] = dists[gene1][gene2] / 2.0;
    branches[parent] = 0.0;
    
    assert(parent == ngenes*2-2);
    
    delete [] restdists;
    delete [] leaves;
}



void reconRoot(Tree *tree, SpeciesTree *stree, int *gene2species)
{
    
}


/*
def reconRoot(gtree, stree, gene2species = gene2species, 
               rootby = "duploss", newCopy=True):
    # make a consistent unrooted copy of gene tree
    if newCopy:
        gtree = gtree.copy()
    treelib.unroot(gtree, newCopy=False)
    treelib.reroot(gtree, 
                   gtree.nodes[util.sort(gtree.leafNames())[0]].parent.name, 
                   onBranch=False, newCopy=False)
    
    
    # make recon root consistent for rerooting tree of the same names
    # TODO: there is the possibility of ties, they are currently broken
    # arbitrarily.  In order to make comparison of reconRooted trees with 
    # same gene names accurate, hashOrdering must be done, for now.
    hashOrderTree(gtree, gene2species)
    
    # get list of edges to root on
    edges = []
    def walk(node):
        edges.append((node, node.parent))
        if not node.isLeaf():
            node.recurse(walk)
            edges.append((node, node.parent))
    for child in gtree.root.children:
        walk(child)
    
    
    # try initial root and recon    
    treelib.reroot(gtree, edges[0][0].name, newCopy=False)
    recon = reconcile(gtree, stree, gene2species)
    events = labelEvents(gtree, recon)     
    
    # find reconciliation that minimizes loss
    minroot = edges[0]
    rootedge = sorted(edges[0])
    if rootby == "dup": 
        cost = countDup(gtree, events)
    elif rootby == "loss":
        cost = len(findLoss(gtree, stree, recon))
    elif rootby == "duploss":
        cost = countDupLoss(gtree, stree, recon, events)
    else:
        raise "unknown rootby value '%s'"  % rootby
    mincost = cost
    
    
    # try rooting on everything
    for edge in edges[1:-1]:
        if sorted(edge) == rootedge:
            continue
        rootedge = sorted(edge)
        
        node1, node2 = edge
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2, "%s %s" % (node1.name, node2.name)
        
        # uncount cost
        if rootby in ["dup", "duploss"]:
            if events[gtree.root] == "dup":
                cost -= 1
            if events[node2] == "dup":
                cost -= 1
        if rootby in ["loss", "duploss"]:
            cost -= len(findLossNode(gtree.root, recon))
            cost -= len(findLossNode(node2, recon))
        
        # new root and recon
        treelib.reroot(gtree, node1.name, newCopy=False)        
        
        recon[node2] = reconcileNode(node2, stree, recon)
        recon[gtree.root] = reconcileNode(gtree.root, stree, recon)
        events[node2] = labelEventsNode(node2, recon)
        events[gtree.root] = labelEventsNode(gtree.root, recon)
        
        if rootby in ["dup", "duploss"]:
            if events[node2] ==  "dup":
                cost += 1
            if events[gtree.root] ==  "dup":
                cost += 1
        if rootby in ["loss", "duploss"]:
            cost += len(findLossNode(gtree.root, recon))
            cost += len(findLossNode(node2, recon))
        
        #print edge[0].name, edge[1].name, cost
        
        # keep track of min cost
        if cost < mincost:
            mincost = cost
            minroot = edge
    
    # root tree by minroot
    if edge != minroot:
        node1, node2 = minroot
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2
        treelib.reroot(gtree, node1.name, newCopy=False)
    
    return gtree
*/



void Tree::reroot(Node *newroot, bool onBranch)
{
    // handle trivial case, newroot is root
    if (root == newroot ||
        (onBranch &&
         root->nchildren == 2 &&
         (root->children[0] == newroot ||
          root->children[1] == newroot)))
        return;
    
    // determine where to stop ascending
    Node *oldroot = root;
    Node *stop1=NULL, *stop2=NULL;
    
    if (isRooted()) {
        stop1 = root->children[0];
        stop2 = root->children[1];
    } else {
        stop1 = root;
    }

    Node *ptr = NULL, *ptr2 = NULL;
    
    if (onBranch) {
        if (isRooted()) {
            // just need to stick current root somewhere else
            oldroot->children[0] = newroot;
            oldroot->children[1] = newroot->parent;

            Node *other = newroot->parent;
            ptr = other;
            ptr2 = oldroot;
            
            if (ptr->children[0] == newroot)
                ptr->children[0] = ptr->parent;
            else
                ptr->children[1] = ptr->parent;
            
            ptr->parent = ptr2;
            
            
        } else {
            // need to add a new node to be root
            assert(0);
        }
    } else {
        if (isRooted()) {
            // need to remove the root node, nad make tribranch
            assert(0);
        } else {
            // just need to swap node positions
            assert(0);
        }
    }
    
    
    // reverse parent child relationships

    while (ptr != stop1 && ptr != stop2) {
        
    }
    
}


/*


def reroot(tree, newroot, onBranch=True):

    unroot(tree, newCopy=False)
    
    if onBranch:
        # add new root in middle of branch
        newNode = TreeNode(tree.newName())
        node1 = tree.nodes[newroot]
        rootdist = node1.dist
        node1.dist = rootdist / 2.0
        newNode.dist = rootdist / 2.0
        node2 = node1.parent
        node2.children.remove(node1)
        tree.addChild(newNode, node1)
        tree.addChild(node2, newNode)
        
        ptr = node2
        ptr2 = newNode
        newRoot = newNode
    else:
        # root directly on node
        ptr2 = tree.nodes[newroot]
        ptr = ptr2.parent
        newRoot = ptr2
    
    newRoot.parent = None
    
    # reverse parent child relationship of all nodes on path node1 to root
    oldroot = tree.root    
    nextDist = ptr2.dist
    ptr2.dist = 0
    while True:
        nextPtr = ptr.parent
        ptr.children.remove(ptr2)
        tree.addChild(ptr2, ptr)
        
        tmp = ptr.dist
        ptr.dist = nextDist
        nextDist = tmp
        
        ptr2 = ptr
        ptr = nextPtr
        
        if nextPtr is None:
            break
    tree.root = newRoot
    
    return tree

*/


// Find Last Common Ancestor
Node *treeLca(SpeciesTree *stree, Node *node1, Node *node2)
{
    int depth1 = stree->depths[node1->name];
    int depth2 = stree->depths[node2->name];
        
    // get nodes to same depth
    if (node1 != node2) {
        while (depth1 > depth2) {
            node1 = node1->parent;
            depth1 = stree->depths[node1->name];
        }
        
        while (depth2 > depth1) {
            node2 = node2->parent;
            depth2 = stree->depths[node2->name];
        }
    }
    
    // walk up both nodes until they meet
    while (node1 != node2) {
        node1 = node1->parent;
        node2 = node2->parent;
    }
    
    return node1;
}


// NOTE: assumes binary species tree
void reconcile_helper(Tree *tree, Node *node, SpeciesTree *stree, int *recon)
{
    // recurse
    for (int i=0; i<node->nchildren; i++)
        reconcile_helper(tree, node->children[i], stree, recon);
    
    if (node->nchildren > 0) {
        int sname1 = recon[node->children[0]->name];
        int sname2 = recon[node->children[1]->name];
    
        // this node's species is lca of children species
        recon[node->name] = treeLca(stree, 
                                    &(stree->nodes[sname1]), 
                                    &(stree->nodes[sname2]))->name;
    }
}


// reconcile a gene tree with a species tree
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon)
{  
    // label gene leaves with their species
    for (int i=0; i<tree->nnodes; i++)
        if (tree->nodes[i].nchildren == 0)
            recon[i] = gene2species[i];
    
    reconcile_helper(tree, tree->root, stree, recon);    
}


// label events for each node in tree
// NOTE: assumes binary gene tree
void labelEvents(Tree *tree, int *recon, int *events)
{
    Node *nodes = tree->nodes;

    for (int i=0; i<tree->nnodes; i++) {
        if (nodes[i].nchildren == 0)
            events[i] = EVENT_GENE;
        else 
        if (recon[i] == recon[nodes[i].children[0]->name] ||
            recon[i] == recon[nodes[i].children[1]->name])
            events[i] = EVENT_DUP;
        else
            events[i] = EVENT_SPEC;
    }
}



//=============================================================================
// conversion functions

// creates a forward tree from a parent tree
void makeFtree(int nnodes, int *ptree, int ***ftree)
{
    *ftree = new int* [nnodes];
    int **ftree2 = *ftree;
    
    // initialize
    for (int i=0; i<nnodes; i++) {
        ftree2[i] = new int [2];
        ftree2[i][0] = -1;
        ftree2[i][1] = -1;
    }
    
    // populate
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];
        
        if (parent != -1) {
            if (ftree2[parent][0] == -1)
                ftree2[parent][0] = i;
            else
                ftree2[parent][1] = i;
        }
    }
}


void freeFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++)
        delete [] ftree[i];
    delete [] ftree;
}


// create a tree object from a parent tree array
void ptree2tree(int nnodes, int *ptree, Tree *tree)
{
    Node *nodes = tree->nodes;
    
    // allocate children
    for (int i=0; i<nnodes; i++) {
        nodes[i].allocChildren(2);
        nodes[i].name = i;
        nodes[i].nchildren = 0;
    }
    
    // store parent and child pointers
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];
        
        if (parent != -1) {
            Node *parentnode = &nodes[parent];            
            parentnode->children[parentnode->nchildren++] = &nodes[i];
            nodes[i].parent = parentnode;
        } else {
            nodes[i].parent = NULL;
        }
    }
    
    // set root
    tree->root = &nodes[nnodes - 1];
}


// create a tree object from a parent tree array
void tree2ptree(Tree *tree, int *ptree)
{
    Node *nodes = tree->nodes;
    int nnodes = tree->nnodes;
    
    for (int i=0; i<nnodes; i++) {
        if (nodes[i].parent)
            ptree[i] = nodes[i].parent->name;
        else
            ptree[i] = -1;
    }
}


//=============================================================================
// Input/output


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}


// write out the newick notation of a tree
void printTree(Tree *tree, Node *node, int depth)
{
    if (node == NULL) {
        if (tree->root != NULL) {
            printTree(tree, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (node->nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d", node->name);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=(\n", node->name);
            
            for (int i=0; i<node->nchildren - 1; i++) {
                printTree(tree, node->children[i], depth+1);
                printf(",\n");
            }
            
            printTree(tree, node->children[node->nchildren-1], depth+1);
            printf("\n");
            
            for (int i=0; i<depth; i++) printf("  ");
            printf(")");
        }
    }
}

// write out the newick notation of a tree
void writeNewick(Tree *tree, string *names, Node *node, int depth)
{
    if (node == NULL) {
        if (tree->root != NULL) {
            writeNewick(tree, names, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (node->nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%s:%f", names[node->name].c_str(), node->dist);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("(\n");
            
            for (int i=0; i<node->nchildren - 1; i++) {
                writeNewick(tree, names, node->children[i], depth+1);
                printf(",\n");
            }
            
            writeNewick(tree, names, node->children[node->nchildren-1], depth+1);
            printf("\n");
            
            for (int i=0; i<depth; i++) printf("  ");
            printf(")");
            
            if (depth > 0)
                printf(":%f", node->dist);
        }
    }
}

