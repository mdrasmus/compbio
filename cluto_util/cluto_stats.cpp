/*******************************************************************************
    Matt Rasmussen
    June 30, 2004
    
    cluto_stats
    
    command line utility for accessing cluto's statistics about a clustering
    solution   

*******************************************************************************/


#include <stdio.h>
#include <cluto.h>
#include <iostream>
#include <sstream>
#include <string>

#include "cluto_io.h"

using namespace std;

#define USAGE "usage: clusterstats <mat file> <part file> <nfeatures>\n\
output: \n\
	internal ids\n\
	internal precents\n\
	external ids\n\
	external precents\n"




int main(int argc, char **argv)
{
    /* parse args */
    if (argc < 4) {
        fprintf(stderr, USAGE);
        return 1;
    }
    
    char *matfile  = argv[1];
    char *partfile = argv[2];
    int nfeatures  = atoi(argv[3]);
    
    FILE *partstream = fopen(partfile, "r");
    FILE *matstream  = fopen(matfile, "r");

    if (matstream == NULL || partstream == NULL) {
        fprintf(stderr, "cannot open input files\n");
        return 1;
    }
    
    int nrows, ncols, nnz;
    int *rowptr, *rowind;
    float *rowval;
    if (!ReadMatrix(matstream, &nrows, &ncols, &nnz, &rowptr, &rowind, &rowval))
    {
        fprintf(stderr, "error reading matrix.\n");
        return 2;
    }
    
    int *part = new int [nrows];
    
    int size = nrows;
    int nparts = 1;
    if (!ReadPart(partstream, part, size, &nparts)) {
        fprintf(stderr, "error readin part file.\n");
        return 3;
    }
    
    // prepare output arrays
    int   *internalids  = new int [nfeatures * nparts];
    float *internalwgts = new float [nfeatures * nparts];
    int   *externalids  = new int [nfeatures * nparts];
    float *externalwgts = new float [nfeatures * nparts];
    
    // calc stats    
    CLUTO_V_GetClusterFeatures(
        nrows, ncols, rowptr, rowind, rowval, 
        CLUTO_SIM_COSINE, CLUTO_ROWMODEL_NONE, CLUTO_COLMODEL_NONE, 1, 
        nparts, part, nfeatures,
        internalids, internalwgts, externalids, externalwgts);
    
    for (int i=0; i<nparts; i++) {
        for (int j=0; j<nfeatures; j++) {
            printf("%d ", internalids[i*nfeatures + j]);
        }
        printf("\n");
        for (int j=0; j<nfeatures; j++) {
            printf("%f ", internalwgts[i*nfeatures + j]);
        }
		printf("\n");
		for (int j=0; j<nfeatures; j++) {
            printf("%d ", externalids[i*nfeatures + j]);
        }
        printf("\n");
        for (int j=0; j<nfeatures; j++) {
            printf("%f ", externalwgts[i*nfeatures + j]);
        }
        printf("\n\n");
    }
}

