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

#define USAGE "usage: cluto_reorder <mat file> <part file> [top|full]\n"


int main(int argc, char **argv)
{
    /* parse args */
    if (argc < 3) {
        fprintf(stderr, USAGE);
        return 1;
    }
    
    char *matfile  = argv[1];
    char *partfile = argv[2];
    int treetype   = CLUTO_TREE_TOP;
    
    if (argc == 4) {
        if (string(argv[3]) == "top")
            treetype = CLUTO_TREE_TOP;
        else
            treetype = CLUTO_TREE_FULL;
    }
    
    
    
    // open input files
    FILE *partstream = fopen(partfile, "r");
    FILE *matstream  = fopen(matfile, "r");

    if (matstream == NULL || partstream == NULL) {
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
    
    int *part = new int [nrows];

    fprintf(stderr, "reading part '%s'...\n", partfile);
    
    // read in part
    int size = nrows;
    int nparts;
    if (!ReadPart(partstream, part, size, &nparts)) {
        fprintf(stderr, "error readin part file.\n");
        return 3;
    }
    
    
    // perform reordering
    int *perm = new int [nrows];    
    ClutoReorder(nrows, ncols, rowptr, rowind, rowval, 
                 treetype, nparts, part, perm);
    
    // print out permutation
    for (int i=0; i<nrows; i++)
        printf("%d\n", perm[i]);

    // clean up
    delete [] perm;
    delete [] part;
    delete [] rowptr;
    delete [] rowind;
    delete [] rowval;
    
    fprintf(stderr, "done.\n");
}


