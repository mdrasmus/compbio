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

using namespace spidir;

// matlab headers
#include <mex.h>


//=============================================================================
// MATLAB gateway function

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    // check argument counts    
    if (nrhs != 5) 
    {
        mexErrMsgTxt("wrong number of arguments.\n"
                     "usage:\n\n"
                     "  ptree    -- int array\n"
                     "  seqs     -- string array\n"
                     "  bgfreq   -- float array \n"
                     "  tsvratio -- float\n"
                     "  maxiter  -- int\n");
    } 
    else if (nlhs != 2) 
    {
        mexErrMsgTxt("must have 2 output arguments\n"
                     "  logl  -- (float) log likelihood\n"
                     "  dists -- (float array) branch lengths\n");
    }

    //================================================
    // parse args
    
    // gene parent tree (ptree)
    int nnodes;
    StackArray<int> ptree;
    if (!getIntArray(prhs[0], &ptree, &nnodes)) mexErrMsgTxt("bad ptree");
    
    
    // sequences (seqs)
    char **seqs = NULL;    
    int nseqs = 0;
    if (!getStringArray(prhs[1], &seqs, &nseqs)) {
        freeStringArray(seqs, nseqs);
        mexErrMsgTxt("bad seqs");
    }
    
    // background base frequency (bgfreq)
    int nbases;    
    StackArray<float> bgfreq;
    if (!getFloatArray(prhs[2], &bgfreq, &nbases)) mexErrMsgTxt("bad bgfreq");
    
    // transition/transversion ratio
    float tsvratio;
    if (!getFloat(prhs[3], &tsvratio)) mexErrMsgTxt("bad tsvratio");
    
    // maximum iterations for optimizing Likelihood
    int maxiter;
    if (!getInt(prhs[4], &maxiter)) mexErrMsgTxt("bad maxiter");


    
    
    //=================================================    
    // call C code    
    ExtendArray<float> dists(nnodes); // = new float [nnodes];
    for (int i=0; i<nnodes; i++) 
        dists[i] = 0.0;
    
    
    float logl = findMLBranchLengthsHky(nnodes, ptree, nseqs, seqs, 
                                        dists, bgfreq, tsvratio, maxiter,
                                        true);
    
    freeStringArray(seqs, nseqs);
    
    // return logl
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[0]) = logl;
    
    // return branch lengths (dists)
    plhs[1] = mxCreateDoubleMatrix(1, nnodes, mxREAL);
    double *ptr = mxGetPr(plhs[1]);
    for (int i=0; i<nnodes; i++)
        ptr[i] = dists[i];
}



