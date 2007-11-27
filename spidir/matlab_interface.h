#ifndef SPIDIR_MATLAB_INTERFACE_H
#define SPIDIR_MATLAB_INTERFACE_H

#include <string>

// matlab headers
#include <mex.h>

bool getInt(const mxArray *arr, int *value);
bool getFloat(const mxArray *arr, float *value);
bool getFloatArray(const mxArray *arr, float **value, int *size);
bool getIntArray(const mxArray *arr, int **value, int *size);
bool getString(const mxArray *arr, char **str);
bool getStringArray(const mxArray *arr, char ***strings, int *size);
void freeStringArray(char **array, int size);


#endif // SPIDIR_MATLAB_INTERFACE_H
