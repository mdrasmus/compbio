/*******************************************************************************
    Matt Rasmussen
    June 30, 2004
    
    cluto_reorder
    
    command line utility for accessing cluto's reordering operations on
    clustering solutions.

*******************************************************************************/


#include <stdio.h>
#include <cluto.h>
#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>

#include "cluto_io.h"
#include "cluto_tree.h"

using namespace std;

#define USAGE "usage: cluto_reorder_tree <mat file> <tree file>\n"


int main(int argc, char **argv)
{
    /* parse args */
    if (argc < 3) {
        fprintf(stderr, USAGE);
        return 1;
    }
    
    char *matfile  = argv[1];
    char *treefile = argv[2];

    
    // open input files
    FILE *treestream = fopen(treefile, "r");
    FILE *matstream  = fopen(matfile, "r");

    if (matstream == NULL || treestream == NULL) {
        fprintf(stderr, "cannot open input files\n");
        return 1;
    }
    
    fprintf(stderr, "reading matrix '%s'...\n", matfile);
    
    // read in matrix
    int nrows, ncols, nnz;
    int *rowptr, *rowind;
    float *rowval;
    if (!ReadMatrix(matstream, &nrows, &ncols, &nnz, &rowptr, &rowind, &rowval))
    {
        fprintf(stderr, "error reading matrix.\n");
        return 2;
    }
    
    fprintf(stderr, "nrows = %d, ncols = %d, nnz = %d\n", nrows, ncols, nnz);
    
    
    int *ptree = new int [2 * nrows];

    fprintf(stderr, "reading tree '%s'...\n", treefile);
    
    // read in part
    int size = nrows;
    int nparts;
    if (!ReadParentTree(treestream, nrows, ptree)) {
        fprintf(stderr, "error reading tree file.\n");
        return 3;
    }
    
    
    
    // perform reordering
    int *perm = new int [nrows]; 
    ClutoReorderTree(nrows, ncols, rowptr, rowind, rowval,
                     ptree, perm);
    
    // print out permutation
    for (int i=0; i<nrows; i++)
        printf("%d\n", perm[i]);

    // clean up
    delete [] perm;
    delete [] ptree;
    delete [] rowptr;
    delete [] rowind;
    delete [] rowval;
    
    fprintf(stderr, "done.\n");
}


