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

#include "ExtendArray.h"

using namespace std;



#define PI 3.1415926



extern int dna2int[256];



class BufferedReader
{
public:
    BufferedReader(FILE *stream) :
        m_stream(stream),
        m_line(0, 10000)
    {}
    
    char *readLine()
    {
        while (!feof(m_stream)) {
            int pos = m_line.size();
            char *ret = fgets(&(m_line.get()[pos]), 
                              m_line.capacity()-m_line.size(), m_stream);
            int readsize = strlen(&(m_line.get()[pos]));
            
            if (ret == NULL)
                return NULL;
            
            if (m_line.size() + readsize < m_line.capacity() - 1)
                return m_line.get();

            assert(m_line.increaseCapacity());
            m_line.setSize(m_line.size() + readsize);
        }
        
        return NULL;
    }
    
protected:
    FILE *m_stream;
    ExtendArray<char> m_line;
};



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
Sequences *readAlignFasta(const char *filename);
bool writeFasta(const char *filename, Sequences *seqs);
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
bool chomp(char *str);

#endif // SPIDIR_COMMON_H
