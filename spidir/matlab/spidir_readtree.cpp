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
    if (nrhs != 1)
    {
        mexErrMsgTxt("wrong number of arguments.\n"
                     "usage:\n\n"
                     "  filename  -- string\n");
    }
    else if (nlhs != 2) 
    {
        mexErrMsgTxt("must have 2 output arguments\n"
                     "  ptree -- (int array) parent tree array\n"
                     "  dists -- (float array) branch lengths\n");
    }

    

    //================================================
    // parse args
    char *filename;
    if (!getString(prhs[0], &filename)) {
        mexErrMsgTxt("bad filename");
    }
    
    
        
    
    //=================================================    
    // read file
    Tree tree;
    if (!tree.readNewick(filename)) {
        mexErrMsgTxt("error in file format");
    }
        
    // convert to ptree
    ExtendArray<float> dists(tree.nnodes);
    ExtendArray<int> ptree(tree.nnodes);
    
    tree2ptree(&tree, ptree);
    tree.getDists(dists);
    
    //=================================================
    // return parent tree (ptree)
    
    plhs[0] = mxCreateDoubleMatrix(1, tree.nnodes, mxREAL);
    double *ptr = mxGetPr(plhs[0]);
    for (int i=0; i<tree.nnodes; i++)
        ptr[i] = ptree[i];
    
    // return branch lengths (dists)
    plhs[1] = mxCreateDoubleMatrix(1, tree.nnodes, mxREAL);
    ptr = mxGetPr(plhs[1]);
    for (int i=0; i<tree.nnodes; i++)
        ptr[i] = dists[i];
}



