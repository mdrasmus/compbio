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

#include "ExtendArray.h"

using namespace std;



#define PI 3.1415926
#define MAX_FLOAT 1e10
#ifndef INFINITY
#   define INFINITY 1e1000
#endif

// indexing a matrix stored as a single array
#define matind(m, i, j) ((m)*(i) + (j))



// convert dna characters into standard numbers
extern int dna2int[256];

// convert standard numbers to dna characters
extern char *int2dna;

// base numbers
enum {
    DNA_A = 0,
    DNA_C = 1,
    DNA_G = 2,
    DNA_T = 3,
    DNA_PURINE,
    DNA_PRYMIDINE
};


// get the base type of a nucleotide
extern int dnatype[];





class BufferedReader
{
public:
    BufferedReader(FILE *stream=NULL, bool autoclose=true) :
        m_stream(stream),
        m_line(0, 10000),
        m_autoclose(autoclose)
    {}
    
    virtual ~BufferedReader()
    {
        if (m_autoclose && m_stream)
            fclose(m_stream);
    }
    
    
    bool open(const char *filename, char *mode, 
              const char *errmsg="cannot read file '%s'\n")
    {
        m_stream = fopen(filename, mode);
        
        if (!m_stream) {
            fprintf(stderr, errmsg, filename);        
            return false;
        }
        return true;
    }
    
    
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
    
    
    void close()
    {
        if (m_stream)
            fclose(m_stream);
        m_stream = NULL;
    }
    
    
protected:
    FILE *m_stream;
    ExtendArray<char> m_line;
    bool m_autoclose;
};



template <class T>
int findval(T *array, int size, const T &val)
{
    for (int i=0; i<size; i++)
        if (array[i] == val)
            return i;
    return -1;
}



// forward declaration
class SpeciesTree;

// spidir parameters
class SpidirParams
{
public:
    SpidirParams(int size, string *_names, 
                 float *_mu, float *_sigma, float _alpha, float _beta) :
        nsnodes(size),
        alpha(_alpha),
        beta(_beta)
    {
        names = new string [nsnodes];
        mu = new float [nsnodes];
        sigma = new float [nsnodes];
        
        for (int i=0; i<nsnodes; i++) {
            if (_names)
                names[i] = _names[i];
            mu[i] = _mu[i];
            sigma[i] = _sigma[i];
        }
    }
    
    ~SpidirParams()
    {
        delete [] names;
        delete [] mu;
        delete [] sigma;
    }
    
    // sorts the parameters to match newnames
    bool order(SpeciesTree *tree);    

    int nsnodes;
    string *names;
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



class Gene2speciesRule
{
public:
    Gene2speciesRule(int rule=PREFIX, string expr="", string species="") :
        rule(rule),
        expr(expr),
        species(species)
    {
    }
    
    enum {
        PREFIX,
        SUFFIX,
        EXACT
    };
    
    int rule;
    string expr;
    string species;
};


class Gene2species
{
public:
    Gene2species() :
        m_rules(0, 20)
    {}
    
    const static string NULL_SPECIES;
    
    bool read(const char *filename);
    string getSpecies(string gene);
    bool getMap(string *genes, int ngenes, string *species, int nspecies, 
                int *map);
    
protected:
    ExtendArray<Gene2speciesRule> m_rules;
    
};


inline float frand()
{
    return rand() / float(RAND_MAX);
}


bool checkSequences(int nseqs, int seqlen, char **seqs);

void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat);
Sequences *readFasta(const char *filename);
Sequences *readAlignFasta(const char *filename);
bool writeFasta(const char *filename, Sequences *seqs);
bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names);
SpidirParams *readSpidirParams(const char* filename);


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
bool inChars(char c, const char *chars);
bool chomp(char *str);
vector<string> split(const char *str, const char *delim, bool multiDelim = true);
void printError(const char *fmt, ...);
char readChar(FILE *stream, int &depth);
char readUntil(FILE *stream, string &token, char *stops, int &depth);
string trim(const char *word);
float readDist(FILE *infile, int &depth);



#endif // SPIDIR_COMMON_H
