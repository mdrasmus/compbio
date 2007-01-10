/*******************************************************************************
    Matt Rasmussen
    June 30, 2004
    
    cluto_graphify
    
    command line utility for accessing cluto's similarity graph creation
    functions.

*******************************************************************************/

#include <stdio.h>
#include <cluto.h>
#include <iostream>
#include <sstream>
#include <string>

#include "cluto_io.h"

using namespace std;

#define USAGE "\
usage: cluto_graphify <smat/mat file> <graph file> \n\
          [-n <nnbrs = 10>] [-type <type>] [-sim <type>]\n\
Graph Types (-type)\n\
  sd - symetric direct\n\
  ad - asymetric direct\n\
Similarity (-sim)\n\
  cos   - cosine similarity\n\
  edist - Euclidean distance\n"


int main(int argc, char **argv)
{
    // my change is right here

    /* parse args */
    if (argc < 2) {
        fprintf(stderr, USAGE);
        return 1;
    }
    
    char *smatfile  = argv[1];
    char *graphfile  = argv[2];
    int nnbrs = 10;
    int graphType = CLUTO_GRMODEL_ASYMETRIC_DIRECT;
    int simType = CLUTO_SIM_COSINE;
    char *gname = "ad";
    char *sname = "cos";

    // parse optional arguments
    for (int i=3; i<argc; i++) {
        // nearest neighbors
        if (string(argv[i]) == "-n" && i < argc - 1) {
            nnbrs = atoi(argv[++i]);
        } 
        // graph type
        else if (string(argv[i]) == "-type" && i < argc - 1) {
            i++;
            gname = argv[i];
            if (string(argv[i]) == "sd") {
                graphType = CLUTO_GRMODEL_SYMETRIC_DIRECT;
            } else if (string(argv[i]) == "ad") {
                graphType = CLUTO_GRMODEL_ASYMETRIC_DIRECT;
            } else {
               fprintf(stderr, "unknown graph type '%s'\n", argv[i]);
               return 1;
            }
        }
        // similarity type
        else if (string(argv[i]) == "-sim" && i < argc - 1) {
            i++;
            sname = argv[i];
            if (string(argv[i]) == "cos") {
                simType = CLUTO_SIM_COSINE;
            } else if (string(argv[i]) == "edist") {
                simType = CLUTO_SIM_EDISTANCE;
            } else {
               fprintf(stderr, "unknown sim type '%s'\n", argv[i]);
               return 1;
            }
        } else {
            fprintf(stderr, "unknown option '%s'\n", argv[i]);
            return 1;
        }
    }
    
    
    // open input files
    FILE *smatstream  = fopen(smatfile, "r");
    FILE *graphstream = fopen(graphfile, "w");

    if (smatstream == NULL || graphstream == NULL) {
        fprintf(stderr, "cannot open input files\n");
        return 1;
    }
    
    printf("reading input matrix '%s'...\n", smatfile);
    
    int nrows, ncols, nnz;
    int *rowptr, *rowind;
    float *rowval;
    if (string(smatfile).find(".smat", 0) < strlen(smatfile) - 4) {
        if (!ReadSMAT(smatstream, &nrows, &ncols, &nnz, &rowptr, &rowind, &rowval))
        {
            fprintf(stderr, "error reading matrix '%s'.\n", smatfile);
            return 2;
        }
    } else {
        if (!ReadMatrix(smatstream, &nrows, &ncols, &nnz, &rowptr, &rowind, &rowval))
        {
            fprintf(stderr, "error reading matrix '%s'.\n", smatfile);
            return 2;
        }
    }
    
    // graph data structure
    int *growptr, *growind;
    float *growval;
    
    printf("creating graph with nnbrs = %d, graph = '%s', sim = '%s'\n", 
            nnbrs, gname, sname);
    
    // calc stats    
    CLUTO_V_GetGraph(
        nrows, ncols, rowptr, rowind, rowval, 
        simType, CLUTO_ROWMODEL_NONE, CLUTO_COLMODEL_NONE, 1,
        graphType, nnbrs, 0, &growptr, &growind, &growval);

    printf("writing output graph '%s'...\n", graphfile);
    
    WriteSMAT(graphstream, nrows, nrows, growptr[nrows], growptr, growind, growval);
    
    delete [] growptr;
    delete [] growind;
    delete [] growval;
}

