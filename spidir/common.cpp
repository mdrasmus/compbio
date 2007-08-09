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

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s)
{
    if (s == 0.0)
        return -INFINITY;
    return - log(s) - log(sqrt(2.0*PI)) - (x-u)*(x-u) / (2.0*s*s);
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