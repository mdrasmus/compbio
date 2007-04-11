/********************************************************************************
    formatfa
    Matt Rasmussen
    Fri Jul 14 21:53:39 EDT 2006



    Creates convenient index files for fasta files
********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

#define MAXLINE 1000


// TODO:  verify that the fasta file has consistent wrapping


inline void chomp(char *line)
{
    int keylen = strlen(line);
    if (line[keylen - 1] == '\n')
        line[keylen - 1] = '\0';
}


bool indexFasta(string filename)
{
    FILE *infile;
    FILE *outfile;
    
    
    // create output filename from input
    string outfilename = filename + ".index";
    
    
    // open input and output files
    if ((infile = fopen(filename.c_str(), "r")) == NULL) {
        fprintf(stderr, "formatfa: cannot open file '%s'\n", filename.c_str());
        return false;
    }
    
    if ((outfile = fopen(outfilename.c_str(), "w")) == NULL) {
        fprintf(stderr, "formatfa: cannot open output file '%s'\n", outfilename.c_str());
        return false;
    }
    
    int lineno = 0;
    char *line = new char [MAXLINE];
    string key;
    long long int start = 0, end = 0, pos = 0;
    long long int linelen;
    
    // process file
    while (fgets(line, MAXLINE, infile)) {
        lineno++;
        linelen = strlen(line);
        pos += linelen;
        
        if (linelen == MAXLINE - 1) {
            fprintf(stderr, "formatfa: line %d too long (MAXLINE=%d)\n", 
                    lineno, MAXLINE);
            exit(1);
        }
        
        if (line[0] == '>') {
            // chomp
            chomp(line);
            
            end = pos - linelen;
            
            // print key and start position
            //fprintf(outfile, "%s\t%d\n", &line[1], ftell(infile));
            if (start != 0) {
                fprintf(outfile, "%s\t%lld\t%lld\n", key.c_str(), start, end);
            }
            key = &line[1];
            start = pos;
        }
        

    }
    
    // print last sequence
    if (start != 0) {
        end = pos;
        fprintf(outfile, "%s\t%lld\t%lld\n", key.c_str(), start, end);
    }
    
    
    // clean up
    delete [] line;
    fclose(infile);
    fclose(outfile);
    
    return true;
}


int main(int argc, char **argv)
{
    // check args
    if (argc < 2) {
        fprintf(stderr, "usage: formatfa <fasta file> ...\n");
        exit(1);
    }
    
    // process fff files
    for (int i=1; i<argc; i++) {
        if (!indexFasta(argv[i])) {
            fprintf(stderr, "formatfa: error processing '%s'\n", argv[i]);
            exit(1);
        }
        fprintf(stderr, "%s\n", argv[i]);
    }
    
    return 0;
}

