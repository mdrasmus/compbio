
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


bool getString(const mxArray *arr, char **str)
{
    int nrows = mxGetM(arr);
    int ncols = mxGetN(arr);
    
    if (nrows != 1)
        return false;
    
    *str = new char [ncols+1];    
    if (mxGetString(arr, *str, ncols+1) == 0)
        return true;
    else
        return false;
}


bool getStringArray(const mxArray *arr, char ***strings, int *size)
{
    int nrows = mxGetM(arr);
    int ncols = mxGetN(arr);
           
    *size = nrows;
    *strings = new char* [nrows];
    
    mxChar *ptr = mxGetChars(arr);
    if (!ptr) {
        *size = 0;
        return false;
    }
    
    for (int i=0; i<nrows; i++) {
        (*strings)[i] = new char [ncols+1];
        
        for (int j=0; j<ncols; j++) {
            (*strings)[i][j] = (char) ptr[j*nrows + i];
        }
        (*strings)[i][ncols] = '\0';
    }
    
    return true;
}


void freeStringArray(char **array, int size)
{
    if (!array)
        return;

    for (int i=0; i<size; i++) {
        if (array[i])
            delete [] array[i];
    }
    delete [] array;
}
