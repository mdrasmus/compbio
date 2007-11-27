/*=============================================================================
    SPIDIR - Matlab interface

=============================================================================*/


// spidir headers
#include "../matlab_interface.h"
#include "../likelihood.h"
#include "../phylogeny.h"
#include "../mldist.h"
#include "../common.h"
#include "../Tree.h"
#include "../ExtendArray.h"
#include "../Matrix.h"

using namespace spidir;

// matlab headers
#include <mex.h>


//=============================================================================
// MATLAB gateway function

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    // check argument counts    
    if (nrhs != 1 || !mxIsDouble(prhs[0])) 
    {
        mexErrMsgTxt("wrong number of arguments.\n"
                     "usage:\n\n"
                     "  distmat  -- float matrix\n");
    }
    else if (nlhs != 2) 
    {
        mexErrMsgTxt("must have 2 output arguments\n"
                     "  ptree -- (int array) parent tree array\n"
                     "  dists -- (float array) branch lengths\n");
    }

    

    //================================================
    // parse args
    const mxArray *matrix = prhs[0];
    double *values = mxGetPr(matrix);
    
    int nrows = mxGetM(matrix);
    int ncols = mxGetN(matrix);
    
    
    if (nrows != ncols)
        mexErrMsgTxt("distance matrix must be square");
    
    
    // read in distance matrix
    int ngenes = nrows;    
    int nnodes = 2*ngenes - 1;    
    Matrix<float> distmat(ngenes, ngenes);
    for (int i=0; i<ngenes; i++) {
        for (int j=0; j<ngenes; j++) {
            distmat[i][j] = values[j*nrows + i];
        }
    }
    
    
    //=================================================    
    // call C code    
    ExtendArray<float> dists(nnodes);
    ExtendArray<int> ptree(nnodes);
    
    neighborjoin(ngenes, distmat.getMatrix(), ptree, dists);
    
    // return parent tree (ptree)
    plhs[0] = mxCreateDoubleMatrix(1, nnodes, mxREAL);
    double *ptr = mxGetPr(plhs[0]);
    for (int i=0; i<nnodes; i++)
        ptr[i] = ptree[i];
    
    // return branch lengths (dists)
    plhs[1] = mxCreateDoubleMatrix(1, nnodes, mxREAL);
    ptr = mxGetPr(plhs[1]);
    for (int i=0; i<nnodes; i++)
        ptr[i] = dists[i];
}



