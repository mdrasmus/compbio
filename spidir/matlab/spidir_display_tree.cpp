/*=============================================================================
    SPIDIR - Matlab interface for displayTree

=============================================================================*/


// spidir headers
#include "../matlab_interface.h"
#include "../common.h"
#include "../Tree.h"
#include "../ExtendArray.h"

using namespace spidir;

// matlab headers
#include <mex.h>



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
    
    char leafname[101];
    
    for (int i=0; i<tree.nnodes; i++) {
        if (tree.nodes[i]->isLeaf()) {
            snprintf(leafname, 100, "%d", tree.nodes[i]->name);
            tree.nodes[i]->leafname = leafname;
        }
    }
    
    char **matrix = NULL;
    int nrows = 0;
    int ncols = 0;
    displayTreeMatrix(&tree, xscale, 2, 
                       &matrix, &nrows, &ncols);  
    
    for (int i=0; i<nrows; i++) 
        printf("%s\n", matrix[i]);
    
    for (int i=0; i<nrows; i++) 
        delete [] matrix[i];
    delete [] matrix;
}



