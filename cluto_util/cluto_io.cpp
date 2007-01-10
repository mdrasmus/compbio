/*******************************************************************************
    Matt Rasmussen
    June 30, 2004
    
    cluto_io
    
    supporting functions for working with cluto io

*******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <cluto.h>
#include <iostream>
#include <sstream>
#include <string>

#include "cluto_io.h"

#define MAX_LINE 4000000

using namespace std;

bool ReadMatrix(FILE *stream, int *nrows, int *ncols, int *nnz, 
    int **rowptr, int **rowind, float **rowval)
{
    char *line =  new char [MAX_LINE];
    if (fgets(line, MAX_LINE, stream) == NULL) {
        delete [] line;
        return false;
    }
    
    // parse header
    if (sscanf(line, "%d %d %d", nrows, ncols, nnz) != 3) {
        int parsed = sscanf(line, "%d %d", nrows, ncols);
        delete [] line;
        
        if (parsed == 2) {
            return ReadDenseMatrix(stream, nrows, ncols, nnz,
                                   rowptr, rowind, rowval, true);
        } else if (parsed == 1) {
            *ncols = *nrows;
            return ReadDenseMatrix(stream, nrows, ncols, nnz,
                                   rowptr, rowind, rowval, true);
        } else {
            return false;
        }
    }
    
	 // allocate space
    *rowptr = new int [*nrows + 1];
    *rowind = new int [*nnz];
    *rowval = new float [*nnz];

    (*rowptr)[0] = 0;
    
    for (int i=0, k=0; i<*nrows; i++) {
        if (fgets(line, MAX_LINE, stream) == NULL) {
            delete [] line;
            return false;
        }

        stringstream tmpstream;
        tmpstream << line;

        int tmp;
        if (sscanf(line, "%d", &tmp) == 1) {
            while (1) {
		          int ind;
		          float val;
                tmpstream >> ind;
                tmpstream >> val;

		          if (tmpstream.eof()) {
			           break;
		          }

		          (*rowind)[k] = ind - 1;
		          (*rowval)[k] = val;
		          k++;
           }
        }

        (*rowptr)[i+1] = k;
    }
    
    delete [] line;
    
    return true;
}



bool ReadDenseMatrix(FILE *stream, int *nrows, int *ncols, int *nnz, 
    int **rowptr, int **rowind, float **rowval, bool usedim=false)
{
    if (!usedim) {
        char *line =  new char [MAX_LINE];
        if (fgets(line, MAX_LINE, stream) == NULL) {
            delete [] line;
            return false;
        }

        // parse header
        if (sscanf(line, "%d %d", nrows, ncols) != 2) {
            delete [] line;
            return false;
        }
        
        delete [] line;
    }
    
	 // allocate space
    *nnz = (*nrows) * (*ncols);
    *rowptr = NULL;
    *rowind = NULL;
    *rowval = new float [*nnz];

    int k = 0;
    
    for (int k=0; k<*nnz; k++) {
        float val;
        
        if (fscanf(stream, "%f", &val) != 1) {
            delete [] *rowval;
            return false;
        }
        
        (*rowval)[k] = val;
    }
    
    return true;
}


bool ReadSMAT(FILE *stream, int *nrows, int *ncols, int *nnz, 
    int **rowptr, int **rowind, float **rowval)
{
    // read header
    if (fscanf(stream, "%d %d %d", nrows, ncols, nnz) != 3) {
        fprintf(stderr, "error in first line\n");
        return false;
    }
    
    int m = *nrows;
    int n = *ncols;
    int z = *nnz;    
    
    // allocate cluto data structures
    *rowptr = new int [m + 1];
    *rowind = new int [z];
    *rowval = new float [z];    
    
    // allocate temps
    int *row = new int [z];
    int *col = new int [z];
    float *val = new float [z];
    int *row_nnz = new int [m];
    
    // init row nnz to zeros
    for (int i=0; i<m; i++)
        row_nnz[i] = 0;
    
    // read in all data
    for (int i=0; i<z; i++) {
        if (fscanf(stream, "%d %d %f", &row[i], &col[i], &val[i]) != 3) {
            fprintf(stderr, "format error line %d\n", i+1);
            return false;
        }
        
        // count nnz per row
        row_nnz[row[i]]++;
    }
    
    // allocate cluto structures
    int *rptr = (*rowptr);
    int *rind = (*rowind);
    float *rval = (*rowval);

    // copy data from temps to cluto structures
    
    // build rowptr
    rptr[0] = 0;
    for (int i=1; i<m+1; i++)
        rptr[i] = rptr[i-1] + row_nnz[i-1];
    
    // copy ind and val
    for (int i=0; i<z; i++) {
        int r = row[i];
        int j = rptr[r] + row_nnz[r] - 1;
        rind[j] = col[i];
        rval[j] = val[i];
        row_nnz[r]--;
    }
    
    // verify my code
    for (int i=0; i<m; i++)
        assert(row_nnz[i] == 0);
    
    // clean up
    delete [] row;
    delete [] col;
    delete [] val;
    delete [] row_nnz;
    
    return true;
}


bool WriteSMAT(FILE *stream, int nrows, int ncols, int nnz, 
    int *rowptr, int *rowind, float *rowval)
{
    // write header
    fprintf(stream, "%d %d %d\n", nrows, ncols, nnz);
    
    // write data
    for (int i=0; i<nrows; i++)
        for (int j=rowptr[i]; j<rowptr[i+1]; j++)
            fprintf(stream, "%d %d %f\n", i, rowind[j], rowval[j]);
    
    return true;
}




bool ReadPart(FILE *stream, int *part, int size, int *nparts)
{
    char line[100];
    *nparts = 0;
    
    for (int i=0; i<size && !feof(stream); i++) {
        fgets(line, 100, stream);
        if (sscanf(line, "%d", &part[i]) != 1)
            return false;
        if (part[i] > *nparts)
            *nparts = part[i];
    }
    
    // nparts = max_i(part[i]) + 1
    (*nparts)++;
    
    return true;
}


bool ReadParentTree(FILE *stream, int nrows, int *ptree)
{
    char line[100];
    
    for (int i=0; i<2*nrows && !feof(stream); i++) {
        fgets(line, 100, stream);
        if (sscanf(line, "%d", &ptree[i]) != 1)
            return false;
    }
    
    return true;
}
