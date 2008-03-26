#ifndef SPIDIR_SEQUENCES_H
#define SPIDIR_SEQUENCES_H

#include <string>

#include "common.h"
#include "ExtendArray.h"

using namespace std;

namespace spidir {

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
    
    void alloc(int _nseqs, int _seqlen)
    {
        // delete previous data
        for (int i=0; i<seqs.size(); i++)
            delete [] seqs[i];
        
        // allocate a blank sequence alignment
        seqs.clear();
        names.clear();
        seqlen = _seqlen;
        nseqs = 0;
        
        for (int i=0; i<_nseqs; i++) {
            char *seq = new char [seqlen+1];
            for (int j=0; j<seqlen; j++)
                seq[j] = ' ';
            seq[seqlen] = '\0';
            append("", seq);
        }
    }
    
    int nseqs;
    int seqlen;
    ExtendArray<char*> seqs;
    ExtendArray<string> names;
};


Sequences *readFasta(const char *filename);
Sequences *readAlignFasta(const char *filename);
void writeFasta(FILE *stream, Sequences *seqs);
bool writeFasta(const char *filename, Sequences *seqs);
bool checkSequences(int nseqs, int seqlen, char **seqs);
void resampleAlign(Sequences *aln, Sequences *aln2);

} // namespace spidir

#endif
