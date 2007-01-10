/*******************************************************************************
    Matt Rasmussen
    June 30, 2004
    
    cluto_tree
    
    supporting functions for working with trees and cluto

*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <cluto.h>
#include <assert.h>

#include "cluto_tree.h"


void ClutoReorder(
    int nrows, int ncols, int *rowptr, int *rowind, float *rowval,
    int treetype, int nparts, int *part, int *perm)
{
    int root;
    
    // build parent tree
    int *ptree = new int [2 * nrows];
    BuildParentTree(
      nrows, ncols, rowptr, rowind, rowval,
      treetype, nparts, part, ptree, &root);
    
    // build forward tree
    int **ftree = new int*[2 * nrows];    
    for (int i=0; i<2*nrows; i++)
        ftree[i] = new int[2];    
    BuildForwardTree(nrows, ptree, ftree);
    
    
    // use cluto reorder this tree in a way that orients similar 
    // subtrees next to each other
    CLUTO_V_TreeReorder(nrows, ncols, rowptr, rowind, rowval, 
        CLUTO_SIM_COSINE, CLUTO_ROWMODEL_NONE, CLUTO_COLMODEL_NONE, 1, 
        0, ptree, ftree);
    
    
    // build permutation
    int index = 0;
    ForwardTreeToPerm(ftree, root, perm, index);
    
    // add unclustered objects to permutation
    for (int i=0; i<nrows; i++)
        if (part[i] == -1)
            perm[index++] = i;
    assert(index == nrows);
    
    // clean up
    delete [] ptree;
    for (int i=0; i<2*nrows; i++)
        delete [] ftree[i];
    delete [] ftree;
}



void ClutoReorderTree(
    int nrows, int ncols, int *rowptr, int *rowind, float *rowval,
    int *ptree, int *perm)
{
    int root;
    
    // calc root position
    // last -1 in ptree counting from the back of the array
    for (root = 2 * nrows - 2; ptree[root] == -1; root--);
    root++;

    
    // build forward tree
    int **ftree = new int*[2 * nrows];    
    for (int i=0; i<2*nrows; i++)
        ftree[i] = new int[2];    
    BuildForwardTree(nrows, ptree, ftree);
    
    
    // use cluto reorder this tree in a way that orients similar 
    // subtrees next to each other
    CLUTO_V_TreeReorder(nrows, ncols, rowptr, rowind, rowval, 
        CLUTO_SIM_COSINE, CLUTO_ROWMODEL_NONE, CLUTO_COLMODEL_NONE, 1, 
        0, ptree, ftree);
    
    
    // build permutation
    int index = 0;
    ForwardTreeToPerm(ftree, root, perm, index);
    
    // clean up
    delete [] ptree;
    for (int i=0; i<2*nrows; i++)
        delete [] ftree[i];
    delete [] ftree;
}


void BuildParentTree(
    int nrows, int ncols, int *rowptr, int *rowind, float *rowval,
    int treetype, int nparts, int *part, int *ptree, int *root)
{ 
    // allocate space for tree
    int *ptree_top = NULL;
    int *ptree_temp;
    float *tsims = new float [2 * nrows];
    float *gains = new float [2 * nrows];
    
    for (int i=0; i<2*nrows; i++)
        ptree[i] = -1;
    
    // allocate extra space if top tree is used
    if (treetype == CLUTO_TREE_TOP) {
        ptree_top = new int [2 * nparts];
        ptree_temp = ptree_top;
    } else {
        ptree_temp =  ptree;
    }
    
    fprintf(stderr, "building tree...\n");
    
    // use cluto to build initial tree for either all nodes (FULL)
    // or just the clusters (TOP).  If a TOP tree is made then the
    // full tree is simulated with the proceeding code.
    CLUTO_V_BuildTree(
        nrows, ncols, rowptr, rowind, rowval, 
        CLUTO_SIM_COSINE, CLUTO_CLFUN_I2, 
        CLUTO_ROWMODEL_NONE, CLUTO_COLMODEL_NONE, 1, 
        treetype, 0, 
        nparts, part, ptree_temp, tsims, gains);
    
    delete [] tsims;
    delete [] gains;
    
    fprintf(stderr, "reordering tree...\n");
    
    // calc root position
    // last -1 in ptree counting from the back of the array
    for (*root = 2 * nrows - 2; ptree[*root] == -1; (*root)--);
    (*root)++;
    
    // build pseudo cluster trees for 'top' tree type
    if (treetype == CLUTO_TREE_TOP) {
        BuildPseudoParentTree(nrows, nparts, part, ptree, ptree_top, root);
        delete [] ptree_top;
    }
}


void BuildPseudoParentTree(
    int nrows, int nparts, int *part, int *ptree, int *ptree_top, int *root)
{
    // copy top tree structure to top of full ptree
    int offset = 2*nrows - 2*nparts;


    // build pseduo cluster trees
    int *parents = new int [nparts];
    int *counts = new int [nparts];
    int *child = new int [nparts];
    int *big_part = new int [2 * nrows];    // partitions for all nodes

    // initialize child counts to zero
    for (int i=0; i<nparts; i++) {
        counts[i] = 0;
        parents[i] = 2*nrows;
    }

    // copy over leaf partition ids from part to big_part
    for (int i=0; i<nrows; i++) 
        big_part[i] = part[i];
    for (int i=nrows; i<2*nrows-1; i++)
        big_part[i] = -1;

    // set first parent to be first non-leaf node
    int parent = nrows;
    
    // iterate through nodes until i equals the next parent
    // at that point, all the nodes in parents[] should be the roots
    // of their respective clusters.  Once we have those roots, we can connect
    // them together with the ptree_top
    for (int i=0; i<parent; i++) {
        // get partition of current node
        int p = big_part[i];

        // do not add unpartitioned objects to the tree
        if (p == -1) {
            ptree[i] = -1;
            continue;
        }

        // do not make a node its own parent
        // this subtree of the tree is done growing
        if (i >= parents[p])
            continue;            

        // if no parent already assigned for this partition
        // then get a new node to be the parent
        // parent becomes part of partition p
        if (counts[p] == 0) {
            parents[p] = parent;
            big_part[parent] = p;
            parent++;
        }

        // assign node i a parent and increment child count for parent 
        counts[p]++;
        ptree[i] = parents[p];
        child[p] = i;
        assert(i < parents[p]);

        // if parent has two children, reset counts[p] to 0
        if (counts[p] == 2) {
            counts[p] = 0;
        }
    }

    // handle special case when one of the parents has one one child
    // this happends when a cluster has only one object
    for (int i=0; i<nparts; i++) {
        if (counts[i] == 1) {
            // replace parent by last (and only) child of the parent
            parents[i] = child[i];
        }
    }

    // now connect cluster parents in parents[] to top ptree
    for (int i=0; i<nparts; i++)
        ptree[parents[i]] = ptree_top[i] + parent - nparts;

    // paste top half of top ptree on to final ptree
    for (int i=nparts, j=parent; i<2*nparts-1; i++, j++)
        ptree[j] = ptree_top[i] + parent - nparts;

    *root = parent + nparts - 2;
    ptree[*root] = -1;

    // clean up TOP tree related data structures
    delete [] parents;
    delete [] counts;
    delete [] child;
    delete [] big_part;
}


void BuildForwardTree(int nrows, int *ptree, int **ftree)
{
    // Initialize forward tree    
    for (int i=0; i<2*nrows; i++) {
        ftree[i][0] = -1;
        ftree[i][1] = -1;
    }

    // Setup forward tree child count array
    int *fnum = new int[2*nrows];
    for (int i=0; i<2*nrows; i++)
        fnum[i] = 0;

    // Build forward tree
    for (int i=0; i<2*nrows-2; i++) 
        if (ptree[i] != -1)
            ftree[ptree[i]][fnum[ptree[i]]++] = i;
}



void ForwardTreeToPerm(int **ftree, int n, int *perm, int &index)
{
    if (ftree[n][0] == -1) {
        // record leaf node
        perm[index++] = n;
    } else {
        // inernal node digs deeper into children
        ForwardTreeToPerm(ftree, ftree[n][0], perm, index);
        ForwardTreeToPerm(ftree, ftree[n][1], perm, index);
    }
}


void ForwardTreeToParentTree(int nrows, int **ftree, int *ptree)
{
    for (int i=0; i<2*nrows; i++) {
        if (ftree[i][0] != -1)
            ptree[ftree[i][0]] = i;
        if (ftree[i][1] != -1)
            ptree[ftree[i][1]] = i;
    }
}

