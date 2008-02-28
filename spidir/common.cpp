/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

// c++ headers
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// spidir headers
#include "common.h"

namespace spidir {


// stream for logging
static FILE *g_logstream = stderr;
static int g_loglevel = LOG_QUIET;


//=============================================================================
// Math



// probability density distribution of the Poisson
float poisson(int x, float lambda)
{
    if (x < 0 || lambda <= 0)
        return 0.0;
    
    float a = 0.0;
    for (float i=1.0; i<x+1; i+=1.0)
        a += log(lambda / i);
    return exp(-lambda + a);
}


/* gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
{
    double x,tmp,ser;
    const static double cof[6]={76.18009172947146,    -86.50532032941677,
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


// Normal distribution.
//
// mu is the mean, and sigma is the standard deviation.
//
float normalvariate(float mu, float sigma)
{
    // Uses Kinderman and Monahan method. Reference: Kinderman,
    // A.J. and Monahan, J.F., "Computer generation of random
    // variables using the ratio of uniform deviates", ACM Trans
    // Math Software, 3, (1977), pp257-260.

    const static float NV_MAGICCONST = 4 * exp(-0.5)/sqrt(2.0);
    float u1, u2, z, zz;

    do {
        u1 = frand();
        u2 = 1.0 - frand();
        z = NV_MAGICCONST*(u1-0.5)/u2;
        zz = z*z/4.0;
    } while (zz > -log(u2));
    
    return mu + z*sigma;
}


float gammavariate(float alpha, float beta)
{
    const static float LOG4 = 1.3862943611198906;
    const static float SG_MAGICCONST = 1.0 + log(4.5);
    
    assert(alpha > 0.0 && beta > 0.0);
    
    // convert beta
    beta = 1.0 / beta;
    
    if (alpha > 1.0) {
        // Uses R.C.H. Cheng, "The generation of Gamma
        // variables with non-integral shape parameters",
        // Applied Statistics, (1977), 26, No. 1, p71-74

        float ainv = sqrt(2.0 * alpha - 1.0);
        float bbb = alpha - LOG4;
        float ccc = alpha + ainv;

        while (1) {
            float u1 = frand();
            if (u1 < 1e-7 || u1 > .9999999)
                continue;
            float u2 = 1.0 - frand();
            float v = log(u1 / (1.0-u1)) / ainv;
            float x = alpha * exp(v);
            float z = u1*u1*u2;
            float r = bbb+ccc*v-x;
            if (r + SG_MAGICCONST - 4.5*z >= 0.0 || r >= log(z))
                return x * beta;
        }
        
    } else if (alpha == 1.0) {
        // expovariate(1)
        float u = 0;
        while (u <= 1e-7)
            u = frand();
        return -log(u) * beta;
        
    } else { 
        // alpha in (0, 1)
        // Uses ALGORITHM GS of Statistical Computing - Kennedy & Gentle
        float x;

        while (1) {
            float u = frand();
            float b = (M_E + alpha)/M_E;
            float p = b*u;
            if (p <= 1.0)
                x = pow(p, (1.0/alpha));
            else
                x = -log((b-p)/alpha);
            float u1 = frand();
            if (p > 1.0)
                if (u1 <= pow(x, (alpha - 1.0)))
                    break;
            else if (u1 <= exp(-x))
                break;
        }
        return x * beta;
    }
}


// Invert a permutation
void invertPerm(int *perm, int *inv, int size)
{
    for (int i=0; i<size; i++)
        inv[perm[i]] = i;
}




//=============================================================================
// input/output

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


bool inChars(char c, const char *chars)
{
   if (!chars)
      return false;
   for (;*chars; chars++)
      if (c == *chars) return true;
   return false;
}


bool chomp(char *str)
{
   int len = strlen(str);
   if (str[len-1] == '\n') {
      str[len-1] = '\0';
      return true;
   } else
      return false;
}


vector<string> split(const char *str, const char *delim, bool multiDelim)
{
    vector<string> tokens;   
    int i=0, j=0;
   
    while (str[i]) {
        // walk to end of next token
        for (; str[j] && !inChars(str[j], delim); j++);
        
        if (i == j)
            break;

        // save token
        tokens.push_back(string(&str[i], j-i));
        
        if (!str[j])
            break;
        j++;
        i = j;
    }
    
    return tokens;
}


string trim(const char *word)
{
    char buf[101];
    sscanf(word, "%100s", buf);
    return string(buf);
}


//=============================================================================
// Errors and Logging

void printError(const char *fmt, ...)
{
    va_list ap;   
    va_start(ap, fmt);
   
    fprintf(stderr, "error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
   
    va_end(ap);
}


void printLog(int level, const char *fmt, ...)
{
    if (level <= g_loglevel) {
        va_list ap;   
        va_start(ap, fmt);
        vfprintf(g_logstream, fmt, ap);
        va_end(ap);
    }
}


bool openLogFile(const char *filename)
{
    FILE *stream = fopen(filename, "w");
    
    if (g_logstream != NULL) {
        openLogFile(stream);
        return true;
    } else {
        return false;
    }
}

void openLogFile(FILE *stream)
{
    g_logstream = stream;
}


void setLogLevel(int level)
{
    g_loglevel = level;
}

bool isLogLevel(int level)
{
    return level <= g_loglevel;
}

void closeLogFile()
{
    fclose(g_logstream);
}

FILE *getLogFile()
{
    return g_logstream;
}



} // namespace spidir
