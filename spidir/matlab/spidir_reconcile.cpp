/*=============================================================================
    SPIDIR - Matlab interface

=============================================================================*/



// matlab headers
#include <mex.h>

// spidir headers
#include "../matlab_interface.h"
#include "../likelihood.h"
#include "../phylogeny.h"
#include "../common.h"
#include "../Tree.h"
#include "../ExtendArray.h"

using namespace spidir;



//=============================================================================
// MATLAB gateway function

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 3) 
    {
        mexErrMsgTxt("not enough arguments.\n"
                     "usage:\n\n"
                     "  ptree  -- int array\n"
                     "  pstree -- int array\n"
                     "  gene2species -- int array\n");
    } 
    else if (nlhs > 2) 
    {
        mexErrMsgTxt("Too many output arguments");
    }
        
    
        
    //================================================
    // parse args
    
    // gene parent tree (ptree)
    int nnodes;
    StackArray<int> ptree;
    if (!getIntArray(prhs[0], &ptree, &nnodes)) mexErrMsgTxt("bad ptree");
        
    // species parent tree (pstree)
    int nsnodes;
    StackArray<int> pstree;
    if (!getIntArray(prhs[1], &pstree, &nsnodes)) mexErrMsgTxt("bad pstree");
    
    // reconciliation (gene2species)
    StackArray<int> gene2species;
    int g2slen;
    if (!getIntArray(prhs[2], &gene2species, &g2slen)) mexErrMsgTxt("bad gene2species");    
    if (g2slen != nnodes)
        mexErrMsgTxt("gene2species is too short")
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    if (!tree.assertTree())
        mexErrMsgTxt("gene tree is invalid");
    
    
    // species tree
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    if (!stree.assertTree())
        mexErrMsgTxt("species tree is invalid");
    stree.setDepths();


    // reconcile gene tree to species tree
    ExtendArray<int> recon(nnodes);
    ExtendArray<int> events(nnodes);

    reconcile(&tree, &stree, gene2species, recon);
    labelEvents(&tree, recon, events);
    
    
    
    // return reconcilation (recon)
    plhs[0] = mxCreateDoubleMatrix(1, nnodes, mxREAL);
    double *ptr = mxGetPr(plhs[0]);
    for (int i=0; i<nnodes; i++)
        ptr[i] = recon[i];

    // return reconcilation (recon)
    plhs[1] = mxCreateDoubleMatrix(1, nnodes, mxREAL);
    ptr = mxGetPr(plhs[1]);
    for (int i=0; i<nnodes; i++)
        ptr[i] = events[i];

}

        



