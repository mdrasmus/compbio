#ifndef SPIDIR_SEQUENCES_H
#define SPIDIR_SEQUENCES_H

#include <string>

#include "common.h"
#include "ExtendArray.h"

using namespace std;

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


Sequences *readFasta(const char *filename);
Sequences *readAlignFasta(const char *filename);
bool writeFasta(const char *filename, Sequences *seqs);
bool checkSequences(int nseqs, int seqlen, char **seqs);


#endif
