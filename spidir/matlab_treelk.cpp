/*=============================================================================
    SPIDIR - Matlab interface

=============================================================================*/



// matlab headers
#include <mex.h>

// spidir headers
#include "matlab_interface.h"
#include "likelihood.h"
#include "phylogeny.h"
#include "common.h"
#include "Tree.h"
#include "ExtendArray.h"

using namespace spidir;



//=============================================================================
// MATLAB gateway function

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs != 9) 
    {
        mexErrMsgTxt("not enough arguments.\n"
                     "usage:\n\n"
                     "  ptree  -- int array\n"
                     "  dists  -- float array\n"
                     "  pstree -- int array\n"
                     "  gene2species -- int array\n"
                     "  mu           -- float array\n"
                     "  sigma        -- float array\n"
                     "  alpha        -- float\n"
                     "  beta         -- float\n"
                     "  generate     -- float \n");
    } 
    else if (nlhs > 1) 
    {
        mexErrMsgTxt("Too many output arguments");
    }
    
    
    // debug output file
    FILE *out = fopen("debug.txt", "a");
    
    
    // return value
    float logl = 0;
    
        
    //================================================
    // parse args
    
    // gene parent tree (ptree)
    int nnodes;
    StackArray<int> ptree;
    if (!getIntArray(prhs[0], &ptree, &nnodes)) mexErrMsgTxt("bad ptree");
    
    // branch lengths (dists)
    StackArray<float> dists;
    if (!getFloatArray(prhs[1], &dists, &nnodes)) mexErrMsgTxt("bad dists");
    
    // species parent tree (pstree)
    int nsnodes;
    StackArray<int> pstree;
    if (!getIntArray(prhs[2], &pstree, &nsnodes)) mexErrMsgTxt("bad pstree");
    
    // reconciliation (gene2species)
    StackArray<int> gene2species;
    if (!getIntArray(prhs[3], &gene2species, &nnodes)) mexErrMsgTxt("bad gene2species");
    
    // mean branch lengths (mu)
    StackArray<float> mu;    
    if (!getFloatArray(prhs[4], &mu, &nsnodes)) mexErrMsgTxt("bad mu");
    
    // standard deviation of branch lengths (sigma)
    StackArray<float> sigma;
    if (!getFloatArray(prhs[5], &sigma, &nsnodes)) mexErrMsgTxt("bad sigma");
    
    // gene-rate distribution (alpha)
    float alpha;
    if (!getFloat(prhs[6], &alpha)) mexErrMsgTxt("bad alpha");

    // gene-rate distribution (beta)
    float beta;
    if (!getFloat(prhs[7], &beta)) mexErrMsgTxt("bad beta");

    // gene-rate 
    float generate;
    if (!getFloat(prhs[8], &generate)) mexErrMsgTxt("bad generate");
    
    // predupprob 
    float predupprob = 1.0;
    
    // dupprob
    float dupprob = 1.0;
    
    // disterror
    float disterror = 0.0;
    
    // errorprob
    float errorprob = 0.0;
    
    
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    if (!tree.assertTree())
        mexErrMsgTxt("gene tree is invalid");
    tree.setDists(dists);

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

    // calculate likelihood
    logl = treelk(nnodes, ptree, dists,
                  nsnodes, pstree, 
                  recon, events,
                  mu, sigma, generate, disterror,
                  predupprob, dupprob, errorprob, alpha, beta);
    
    fclose(out);
    
    /* Create return value */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = logl;
}



