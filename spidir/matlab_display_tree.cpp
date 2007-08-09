/*=============================================================================
    SPIDIR - Matlab interface for displayTree

=============================================================================*/



// matlab headers
#include <mex.h>

// spidir headers
#include "matlab_interface.h"
#include "common.h"
#include "Tree.h"
#include "ExtendArray.h"

using namespace spidir;



//=============================================================================
// MATLAB gateway function

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 3) 
    {
        mexErrMsgTxt("not enough arguments.\n"
                     "usage:\n\n"
                     "  ptree   -- int array\n"
                     "  dists   -- float array\n"
                     "  scaling -- float\n");
    } 
    else if (nlhs > 0) 
    {
        mexErrMsgTxt("Too many output arguments");
    }
    
    
    // debug output file
    FILE *out = fopen("debug.txt", "a");
    
    
        
    //================================================
    // parse args
    
    // gene parent tree (ptree)
    int nnodes;
    StackArray<int> ptree;
    if (!getIntArray(prhs[0], &ptree, &nnodes)) mexErrMsgTxt("bad ptree");
    
    // branch lengths (dists)
    StackArray<float> dists;
    if (!getFloatArray(prhs[1], &dists, &nnodes)) mexErrMsgTxt("bad dists");
    
    // scaling
    float xscale;
    if (!getFloat(prhs[2], &xscale)) mexErrMsgTxt("bad scaling");    
    
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    if (!tree.assertTree())
        mexErrMsgTxt("gene tree is invalid");
    tree.setDists(dists);
    
    displayTree(&tree, out, xscale);
    
    fclose(out);
}



