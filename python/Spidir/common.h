#ifndef SPIDIR_COMMON_H
#define SPIDIR_COMMON_H

/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "spidir.h"




#define PI 3.1415926



extern int dna2int[256];



// spidir parameters
class SpidirParams
{
public:
    SpidirParams(int size, float *_mu, float *_sigma, float _alpha, float _beta) :
        nsnodes(size),
        alpha(_alpha),
        beta(_beta)
    {
        mu = new float [nsnodes];
        sigma = new float [nsnodes];
        
        for (int i=0; i<nsnodes; i++) {
            mu[i] = _mu[i];
            sigma[i] = _sigma[i];
        }
    }
    
    ~SpidirParams()
    {
        delete [] mu;
        delete [] sigma;
    }

    int nsnodes;
    float *mu;
    float *sigma;
    float alpha;
    float beta;
};



void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat);

//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s);
double gammln(double xx);
float gammalog(float x, float a, float b);


//=============================================================================
// debug

void printIntArray(int *array, int size);
void printFloatArray(float *array, int size);

#endif // SPIDIR_COMMON_H
