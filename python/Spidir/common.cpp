/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "spidir.h"
#include "common.h"


#define matind(m, i, j) ((m)*(i) + (j))
#define DNA_A 0
#define DNA_C 1
#define DNA_G 2
#define DNA_T 3


#define MAX_COST 1000000000



int dna2int [256] = 
{
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 9
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 19
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 29
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 39
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 49
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 59
    -1, -1, -1, -1, -1,  0, -1,  1, -1, -1,   // 69
    -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,   // 79
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1,   // 89
    -1, -1, -1, -1, -1, -1, -1,  0, -1,  1,   // 99
    -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,   // 109
    -1, -1, -1, -1, -1, -1,  3, -1, -1, -1,   // 119
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 129
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 139
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 149
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 159
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 169
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 179
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 189
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 199
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 209
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 219
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 229
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 239
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 249
    -1, -1, -1, -1, -1, -1
};


//=============================================================================
// Phylogeny functions


// Find Last Common Ancestor
int treeLca(SpeciesTree *stree, int node1, int node2)
{
    int depth1 = stree->depths[node1];
    int depth2 = stree->depths[node2];
    Node *nodes = stree->nodes;
    
    // get nodes to same depth
    if (node1 != node2) {
        while (depth1 > depth2) {
            node1 = nodes[node1].parent;
            depth1 = stree->depths[node1];
        }
        
        while (depth2 > depth1) {
            node2 = nodes[node2].parent;
            depth2 = stree->depths[node2];
        }
    }
    
    // walk up both nodes until they meet
    while (node1 != node2) {
        node1 = nodes[node1].parent;
        node2 = nodes[node2].parent;
    }
    
    return node1;
}


// NOTE: assumes binary species tree
void reconcile_helper(Tree *tree, int node, SpeciesTree *stree, int *recon)
{
    Node *n = &tree->nodes[node];

    // recurse
    for (int i=0; i<n->nchildren; i++)
        reconcile_helper(tree, n->children[i], stree, recon);
    
    if (n->nchildren > 0) {
        // this node's species is lca of children species
        recon[node] = treeLca(stree, recon[n->children[0]], 
                                     recon[n->children[1]]);
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
        if (recon[i] == recon[nodes[i].children[0]] ||
            recon[i] == recon[nodes[i].children[1]])
            events[i] = EVENT_DUP;
        else
            events[i] = EVENT_SPEC;
    }
}


//=============================================================================
// Parsimony algorithm


float subcost[4][4] = {
    {0, 1, 1, 1},
    {1, 0, 1, 1},
    {1, 1, 0, 1},
    {1, 1, 1, 0}
};


struct ParsimonyCell
{   
    float cost;
    float leftcost;
    float rightcost;
    int leftbase;
    int rightbase;
    bool gap;
};


// assume binary tree
void parsimony_helper(Tree *tree, int node, int nseqs, char **seqs, 
                    ParsimonyCell *table)
{
    for (node=nseqs; node<tree->nnodes; node++) {
        int left = tree->nodes[node].children[0];
        int right = tree->nodes[node].children[1];       
        
        // process this node
        for (int a=0; a<4; a++) {
            int minleft = 0, minright = 0;
            float minleftcost = MAX_COST, minrightcost = MAX_COST;
            float leftsub = 0;
            float rightsub = 0;
            float leftmatch = 2;
            float rightmatch = 2;
            
            for (int b=0; b<4; b++) {
                float sub = subcost[a][b];
                float leftcost = table[matind(4, left, b)].cost + sub;
                float rightcost = table[matind(4, right, b)].cost + sub;
                
                if (leftcost < minleftcost ||
                    (leftcost == minleftcost &&
                     (rand() / float(RAND_MAX)) < (1/leftmatch)))
                {
                    minleftcost = leftcost;
                    minleft = b;
                    leftsub = sub;
                    if (leftcost == minleftcost)
                        leftmatch += 1;                    
                    else
                        leftmatch = 2;
                }
                
                
                if (rightcost < minrightcost ||
                    (rightcost == minrightcost && 
                           (rand() / float(RAND_MAX)) < (1/rightmatch)))
                {
                    minrightcost = rightcost;
                    minright = b;
                    rightsub = sub;
                    
                    if (rightcost == minrightcost)
                        rightmatch += 1;
                    else
                        rightmatch = 2;
                }
            }
            
            // save cost and pointers
            int k = matind(4, node, a);
            table[k].cost = minleftcost + minrightcost;
            table[k].leftcost = leftsub;
            table[k].rightcost = rightsub;
            table[k].leftbase = minleft;
            table[k].rightbase = minright;
            table[k].gap = table[matind(4, left, 0)].gap && \
                           table[matind(4, right, 0)].gap;
        }
    }
}


void getParsimonyCost(Tree *tree, int node, int base, 
                      ParsimonyCell *table, float *dists)
{
    if (tree->nodes[node].nchildren > 0) {
        int left = tree->nodes[node].children[0];
        int right = tree->nodes[node].children[1];
        
        dists[left] += table[matind(4, node, base)].leftcost;
        dists[right] += table[matind(4, node, base)].rightcost;
        
        // recurse
        getParsimonyCost(tree, left, table[matind(4, node, base)].leftbase, 
                         table, dists);
        getParsimonyCost(tree, right, table[matind(4, node, base)].rightbase, 
                         table, dists);
    }
}



void parsimony(int nnodes, int *ptree, int nseqs, char **seqs, float *dists)
{
    int seqlen = strlen(seqs[0]);
    
    // check seqs
    // CHANGE N's to gaps
    for (int i=0; i<nseqs; i++) {
        for (int j=0; j<seqlen; j++) {
            if (seqs[i][j] == 'N' || seqs[i][j] == 'n')
                seqs[i][j] = '-';
            assert(seqs[i][j] == '-' || dna2int[seqs[i][j]] != -1);
        }
    }
    

    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    
    // allocate dynamic table
    ParsimonyCell *table = new ParsimonyCell [nnodes * 4];
    int *gapless = new int [nnodes];
    
    // initalize distances
    for (int i=0; i<nnodes; i++) {
        dists[i] = 0.0;
        gapless[i] = 0;
    }
    
    
    for (int i=0; i<seqlen; i++) {
        
        // initialize leaves
        // iterate just over the leaves
        
        for (int j=0; j<nseqs; j++) {
            if (seqs[j][i] == '-') {
                // gap
                for (int k=0; k<4; k++) {             
                    table[matind(4, j, k)].cost = 0;
                    table[matind(4, j, k)].gap = true;
                }
            } else {
                for (int k=0; k<4; k++) {
                    table[matind(4, j, k)].cost = MAX_COST;
                    table[matind(4, j, k)].gap = false;
                }
                table[matind(4, j, dna2int[seqs[j][i]])].cost = 0;
            }
        }
        
        // populate cost table
        parsimony_helper(&tree, tree.root, nseqs, seqs, table);
        
        
        // find min cost at root
        float mincost = MAX_COST;
        int minbase= 0;
        int root = tree.root;
        dists[root] = 0;
        
        for (int a=0; a<4; a++) {
            if (table[matind(4, root, a)].cost < mincost) {
                mincost = table[matind(4, root, a)].cost;
                minbase = a;
            }
        }
        
        // add up dist
        getParsimonyCost(&tree, root, minbase, table, dists);
        
        // add up ungapped chars
        for (int j=0; j<nnodes; j++) {
            gapless[j] += table[matind(4, j, 0)].gap ? 0 : 1;
        }
    }
    
    // divide subsitutions by number of sites
    for (int i=0; i<nnodes; i++)
        dists[i] /= gapless[i];
        //float(seqlen);
    // TODO: need to protect against divide by zero
    
    
    // cleanup
    delete [] table;
    delete [] gapless;
}


//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s)
{
    if (s < 1e-10 || u < 1e-10) {
        return -INFINITY;
    }    
    return - logf(s * sqrt(2.0*PI)) - (x-u)*(x-u) / (2.0*s*s);
    //return log(1.0/(s * sqrt(2.0*PI)) * exp(- (x-u)*(x-u) / (2.0 * s*s)));
}



/* gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009172947146,    -86.50532032941677,
                          24.01409824083091,    -1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.5066282746310005*ser);
}



float gammalog(float x, float a, float b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}

//=============================================================================
// data structure manipulation

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
    
        nodes[i].parent = parent;
        
        if (parent != -1)
            nodes[parent].children[nodes[parent].nchildren++] = i;
    }
    
    // set root
    tree->root = nnodes - 1;
}


// create a tree object from a parent tree array
void tree2ptree(Tree *tree, int *ptree)
{
    Node *nodes = tree->nodes;
    int nnodes = tree->nnodes;
    
    for (int i=0; i<nnodes; i++)
        ptree[i] = nodes[i].parent;
}


//=============================================================================
// debug

void printIntArray(int *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%d ", array[i]);
    printf("\n");
}

void printFloatArray(float *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%f ", array[i]);
    printf("\n");
}


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}


// write out the newick notation of a tree
void printTree(Tree *tree, int node, int depth)
{
    if (node == -1) {
        if (tree->root != -1) {
            printTree(tree, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (tree->nodes[node].nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d", node);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=(\n", node);
            
            int nchildren = tree->nodes[node].nchildren;
            for (int i=0; i<nchildren - 1; i++) {
                printTree(tree, tree->nodes[node].children[i], depth+1);
                printf(",\n");
            }
            
            printTree(tree, tree->nodes[node].children[nchildren-1], depth+1);
            printf("\n");
            
            for (int i=0; i<depth; i++) printf("  ");
            printf(")");
        }
    }
}

// write out the newick notation of a tree
void writeNewick(Tree *tree, int node, int depth)
{
    if (node == -1) {
        if (tree->root != -1) {
            printTree(tree, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (tree->nodes[node].nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d", node);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("(\n");
            
            int nchildren = tree->nodes[node].nchildren;
            for (int i=0; i<nchildren - 1; i++) {
                writeNewick(tree, tree->nodes[node].children[i], depth+1);
                printf(",\n");
            }
            
            writeNewick(tree, tree->nodes[node].children[nchildren-1], depth+1);
            printf("\n");
            
            for (int i=0; i<depth; i++) printf("  ");
            printf(")");
        }
    }
}
