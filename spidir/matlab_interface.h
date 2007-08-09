#ifndef SPIDIR_MATLAB_INTERFACE_H
#define SPIDIR_MATLAB_INTERFACE_H

// matlab headers
#include <mex.h>

bool getInt(const mxArray *arr, int *value);
bool getFloat(const mxArray *arr, float *value);
bool getFloatArray(const mxArray *arr, float **value, int *size);
bool getIntArray(const mxArray *arr, int **value, int *size);


#endif // SPIDIR_MATLAB_INTERFACE_H
