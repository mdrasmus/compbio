#ifndef CLUTO_READ_H
#define CLUTO_READ_H

#define MAX_LINE 4000000

bool ReadMatrix(FILE *stream, int *nrows, int *ncols, int *nnz, 
    int **rowptr, int **rowind, float **rowval);
bool ReadDenseMatrix(FILE *stream, int *nrows, int *ncols, int *nnz, 
    int **rowptr, int **rowind, float **rowval, bool usedim);
bool ReadSMAT(FILE *stream, int *nrows, int *ncols, int *nnz, 
    int **rowptr, int **rowind, float **rowval);
bool WriteSMAT(FILE *stream, int nrows, int ncols, int nnz, 
    int *rowptr, int *rowind, float *rowval);


bool ReadPart(FILE *stream, int *part, int size, int *npart);
bool ReadParentTree(FILE *stream, int nrows, int *ptree);

#endif
