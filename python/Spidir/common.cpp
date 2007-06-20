/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "spidir.h"
#include "common.h"


#define DNA_A 0
#define DNA_C 1
#define DNA_G 2
#define DNA_T 3



int dna2int [256] = 
{
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 9
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 19
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 29
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 39
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 49
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 59
    -1, -1, -1, -1, -1,  0, -1,  1, -1, -1,   // 69
    -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,   // 79
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1,   // 89
    -1, -1, -1, -1, -1, -1, -1,  0, -1,  1,   // 99
    -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,   // 109
    -1, -1, -1, -1, -1, -1,  3, -1, -1, -1,   // 119
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 129
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 139
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 149
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 159
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 169
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 179
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 189
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 199
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 209
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 219
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 229
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 239
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 249
    -1, -1, -1, -1, -1, -1
};


// calculate the pairwise distances between sequences
// NOTE: simple version implemented first
void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat)
{
    for (int i=0; i<nseqs; i++) {
        distmat[i][i] = 0.0;
    
        for (int j=i+1; j<nseqs; j++) {
            float changes = 0.0;
            int len = 0;
            
            for (int k=0; k<seqlen; k++) {
                if (seqs[i][k] != '-' && seqs[j][k] != '-') {
                    len++;                
                    if (seqs[i][k] != seqs[j][k]) {
                        changes += 1;
                    }
                }
            }
            
            assert(len > 0);
            
            distmat[i][j] = changes / len;
            distmat[j][i] = changes / len;
        }
    }
}


//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s)
{
    if (s < 1e-10 || u < 1e-10) {
        return -INFINITY;
    }    
    return - logf(s * sqrt(2.0*PI)) - (x-u)*(x-u) / (2.0*s*s);
    //return log(1.0/(s * sqrt(2.0*PI)) * exp(- (x-u)*(x-u) / (2.0 * s*s)));
}



/* gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009172947146,    -86.50532032941677,
                          24.01409824083091,    -1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.5066282746310005*ser);
}



float gammalog(float x, float a, float b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}


//=============================================================================
// debug

void printIntArray(int *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%d ", array[i]);
    printf("\n");
}

void printFloatArray(float *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%f ", array[i]);
    printf("\n");
}


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}

