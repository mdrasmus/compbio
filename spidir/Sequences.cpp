
// spidir headers
#include "spidir.h"
#include "Sequences.h"

namespace spidir {


Sequences *readFasta(const char *filename)
{
    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "cannot read file '%s'\n", filename);
        return NULL;
    }
    
    BufferedReader reader(infile);
    char *line;
    
    Sequences *seqs = new Sequences();
    string key;
    ExtendArray<char> seq(0, 10000);

    
    while ((line = reader.readLine())) {
        chomp(line);
        
        if (line[0] == '>') {
            if (seq.size() > 0) {  
                // add new sequence
                seq.append('\0');
                seqs->append(key, seq.detach());
            }
        
            // new key found
            key = string(&line[1]);
        } else {
            seq.extend(line, strlen(line));
        }
        
    }
    
    // add last sequence
    if (seq.size() > 0) {
        seq.append('\0');
        seqs->append(key, seq.detach());
    }
    
    return seqs;
}


Sequences *readAlignFasta(const char *filename)
{
    Sequences *seq = readFasta(filename);
    if (!seq)
        return NULL;
    seq->setAlignLength();
    return seq;
}


bool writeFasta(const char *filename, Sequences *seqs)
{
    FILE *stream = NULL;
    
    if ((stream = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "cannot open '%s'\n", filename);
        return false;
    }

    for (int i=0; i<seqs->nseqs; i++) {
        fprintf(stream, ">%s\n", seqs->names[i].c_str());
        fprintf(stream, "%s\n", seqs->seqs[i]);
    }
    
    fclose(stream);
    return true;
}



// ensures that all characters in the alignment are sensible
// TODO: do not change alignment (keep Ns)
bool checkSequences(int nseqs, int seqlen, char **seqs)
{
    // check seqs
    // CHANGE N's to gaps
    for (int i=0; i<nseqs; i++) {
        for (int j=0; j<seqlen; j++) {
            if (seqs[i][j] == 'N' || seqs[i][j] == 'n')
                // treat Ns as gaps
                seqs[i][j] = '-';
            if (seqs[i][j] != '-' &&
                dna2int[(int) (unsigned char) seqs[i][j]] == -1)
            {
                // a unknown character is in the alignment
                return false;
            }
        }
    }
    
    return true;
}



} // namespace spidir
