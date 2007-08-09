
#include "matlab_interface.h"


//=============================================================================
// MATLAB parsing

bool getInt(const mxArray *arr, int *value)
{
    int nrows = mxGetM(arr);
    int ncols = mxGetN(arr);
    
    if (!mxIsDouble(arr) || nrows != 1 || ncols != 1)
        return false;
    *value = (int) mxGetPr(arr)[0];
    return true;
}

bool getFloat(const mxArray *arr, float *value)
{
    int nrows = mxGetM(arr);
    int ncols = mxGetN(arr);
    
    if (!mxIsDouble(arr) || nrows != 1 || ncols != 1)
        return false;
    *value = mxGetPr(arr)[0];
    return true;
}


bool getFloatArray(const mxArray *arr, float **value, int *size)
{
    int nrows = mxGetM(arr);
    int ncols = mxGetN(arr);
    
    if (!mxIsDouble(arr) || nrows != 1)
        return false;
    
    *size = ncols;
    *value = new float [ncols];
    double *ptr = mxGetPr(arr);
    
    for (int i=0; i<ncols; i++) {
        (*value)[i] = ptr[i];
    }
    
    return true;    
}


bool getIntArray(const mxArray *arr, int **value, int *size)
{
    int nrows = mxGetM(arr);
    int ncols = mxGetN(arr);
    
    if (!mxIsDouble(arr) || nrows != 1)
        return false;
    
    *size = ncols;
    *value = new int [ncols];
    double *ptr = mxGetPr(arr);
    
    for (int i=0; i<ncols; i++) {
        (*value)[i] = (int) ptr[i];
    }
    
    return true;
}

