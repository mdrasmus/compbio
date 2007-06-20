//=============================================================================
// SPIDIR Tree datastructure

#include <assert.h>
#include <stdio.h>


#include "Tree.h"
#include "Matrix.h"

#define MIN_FLOAT -1e10


//=============================================================================
// phylogeny functions


// Neighbor-joining algorithm
void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches)
{
    Matrix<float> dists(ngenes*2-1, ngenes*2-1);
    float *restdists = new float [ngenes*2-1];
    bool *leaves = new bool [ngenes];
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
        float low = MIN_FLOAT;
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
    ptree[gene1] = parent;
    ptree[parent] = -1;
    branches[gene1] = dists[gene1][gene2] / 2.0;
    branches[gene2] = dists[gene1][gene2] / 2.0;
    branches[parent] = 0.0;
    
    assert(parent == ngenes*2);
    
    delete [] restdists;
    delete [] leaves;
}


/*
def neighborjoin(distmat, genes, usertree=None):
    """Neighbor joining algorithm"""
    
    tree = treelib.Tree()
    leaves = {}
    dists = util.Dict(2, None)
    restdists = {}
    
    
    # initialize distances
    for i in range(len(genes)):
        r = 0
        for j in range(len(genes)):
            dists[genes[i]][genes[j]] = distmat[i][j]
            r += distmat[i][j]
        restdists[genes[i]] = r / (len(genes) - 2)
        
    # initialize leaves
    for gene in genes:
        tree.add(treelib.TreeNode(gene))
        leaves[gene] = 1
    
    # if usertree is given, determine merging order
    merges = []
    newnames = {}
    if usertree != None:
        def walk(node):
            if not node.isLeaf():
                assert len(node.children) == 2, \
                    Exception("usertree is not binary")
            
                for child in node:
                    walk(child)
                merges.append(node)
                newnames[node] = len(merges)
            else:
                newnames[node] = node.name
        walk(usertree.root)
        merges.reverse()
    
    # join loop
    while len(leaves) > 2:
        # search for closest genes
        if not usertree:
            low = util.INF
            lowpair = (None, None)
            leaveslst = leaves.keys()

            for i in range(len(leaves)):
                for j in range(i+1, len(leaves)):
                    gene1, gene2 = leaveslst[i], leaveslst[j]
                    dist = dists[gene1][gene2] - restdists[gene1] \
                                               - restdists[gene2]
                    if dist < low:
                        low = dist
                        lowpair = (gene1, gene2)
        else:
            node = merges.pop()
            lowpair = (newnames[node.children[0]],
                       newnames[node.children[1]])
        
        # join gene1 and gene2
        gene1, gene2 = lowpair
        parent = treelib.TreeNode(tree.newName())
        tree.addChild(parent, tree.nodes[gene1])
        tree.addChild(parent, tree.nodes[gene2])
        
        # set distances
        tree.nodes[gene1].dist = (dists[gene1][gene2] + restdists[gene1] - 
                                  restdists[gene2]) / 2.0
        tree.nodes[gene2].dist = dists[gene1][gene2] - tree.nodes[gene1].dist
        
        # gene1 and gene2 are no longer leaves
        del leaves[gene1]
        del leaves[gene2]
        
        gene3 = parent.name
        r = 0
        for gene in leaves:
            dists[gene3][gene] = (dists[gene1][gene] + dists[gene2][gene] -
                                  dists[gene1][gene2]) / 2.0
            dists[gene][gene3] = dists[gene3][gene]
            r += dists[gene3][gene]
        leaves[gene3] = 1
        
        if len(leaves) > 2:
            restdists[gene3] = r / (len(leaves) - 2)
    
    # join the last two genes into a tribranch
    gene1, gene2 = leaves.keys()
    if type(gene1) != int:
        gene1, gene2 = gene2, gene1
    tree.addChild(tree.nodes[gene1], tree.nodes[gene2])
    tree.nodes[gene2].dist = dists[gene1][gene2]
    tree.root = tree.nodes[gene1]

    # root tree according to usertree    
    if usertree != None and treelib.isRooted(usertree):
        roots = set([newnames[usertree.root.children[0]],
                     newnames[usertree.root.children[1]]])
        newroot = None
        for child in tree.root.children:
            if child.name in roots:
                newroot = child
        
        assert newroot != None
        
        treelib.reroot(tree, newroot.name, newCopy=False)
    
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
void writeNewick(Tree *tree, Node *node, int depth)
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
                writeNewick(tree, node->children[i], depth+1);
                printf(",\n");
            }
            
            writeNewick(tree, node->children[node->nchildren-1], depth+1);
            printf("\n");
            
            for (int i=0; i<depth; i++) printf("  ");
            printf(")");
        }
    }
}

