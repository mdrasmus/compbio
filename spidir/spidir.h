#ifndef SPIDIR_H
#define SPIDIR_H

#include <string>

using namespace std;

namespace spidir {


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



//=============================================================================
// Distance Matrices

void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat);
bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names);





//=============================================================================
// SPIDIR parameters


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


SpidirParams *readSpidirParams(const char* filename);


} // namespace spidir

#endif
