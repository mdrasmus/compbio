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
#include <string>
#include <vector>
#include <algorithm>
// #include <cstddef>

using namespace std;


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


template <class T>
T *resize(T *array, size_t oldsize, size_t newsize)
{
   T *tmp = new T[newsize];

    if (oldsize == 0)
        return tmp;

   copy(array, array + oldsize, tmp);
   
   delete [] array;
   return tmp;
}


/*
    len  -- the length of the data
    size -- the size of the allocated array
    
    easy, detachable wrapper for arrays
*/
template <class ValueType>
class ExtendArray
{
public:
    ExtendArray(ValueType *_data=NULL, int _len=0, int size=0,
                int _minsize=40) :
        data(_data),
        len(_len),
        datasize(size),
        minsize(_minsize)
    {
        if (data == NULL && datasize != 0) {
            data = new ValueType [datasize];
        }
    }
    
    ~ExtendArray()
    {
        if (data)
            delete [] data;
    }
    
    void detach()
    {
        data = NULL;
        len = 0;
        datasize = 0;
    }
    
    bool alloc(int needed=0)
    {
        int oldsize = datasize;
        
        if (datasize < minsize)
            datasize = minsize;
        
        while (datasize < needed)
            datasize *= 2;
        
        data = resize(data, oldsize, datasize);
        
        return true;
    }
    
    bool ensureSize(int needed)
    {
        if (needed <= datasize)
            return true;
        return alloc(needed);
    }
    
    void append(ValueType &val)
    {
        ensureSize(len + 1);
        data[len++] = val;
    }
    
    void extend(ValueType *vals, int nvals)
    {
        ensureSize(len + nvals);
        
        for (int i=0; i<nvals; i++)
            data[len++] = vals[i];
    }
    
    inline ValueType &operator[](const int i)
    {
        return data[i];
    }
    
    inline int size()
    {
        return len;
    }
    
    inline void setSize(int _size)
    {
        len = _size;
    }
    
    inline ValueType *getArray()
    {
        return data;
    }
    
protected:    
    ValueType *data;
    int len;
    int datasize;
    int minsize;
};


class Sequences
{
public:
    Sequences(int nseqs=0, int seqlen=0, char **seqs=NULL) :
        nseqs(nseqs),
        seqlen(seqlen)
    {
    }
    
    ~Sequences()
    {
        for (int i=0; i<seqs.size(); i++)
            delete [] seqs[i];
    }
    
    void append(string name, char *seq)
    {
        names.append(name);
        seqs.append(seq);
        nseqs++;
    }
    
    void setAlignLength()
    {
        seqlen = strlen(seqs[0]);
    }
    
    int nseqs;
    int seqlen;
    ExtendArray<char*> seqs;
    ExtendArray<string> names;
};

bool checkSequences(int nseqs, int seqlen, char **seqs);

void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat);
Sequences *readFasta(const char *filename);
bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names);


//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s);
double gammln(double xx);
float gammalog(float x, float a, float b);


//=============================================================================
// input/output

void printIntArray(int *array, int size);
void printFloatArray(float *array, int size);
int readLine(FILE *stream, char **line, int *size);
bool chomp(char *str);

#endif // SPIDIR_COMMON_H
