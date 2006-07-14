/********************************************************************************
    formatfa
    Matt Rasmussen    
    Fri Jun 30 13:17:02 EDT 2006


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
    
    char *line = new char [MAXLINE];
    
    
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
    
    // process file
    while (fgets(line, MAXLINE, infile)) {
        lineno++;
        
        if (strlen(line) == MAXLINE - 1) {
            fprintf(stderr, "formatfa: line %d too long\n", lineno);
            exit(1);
        }
        
        if (line[0] == '>') {
            // chomp
            chomp(line);
            fprintf(outfile, "%s\t%d\n", &line[1], ftell(infile));
        }
    }
    
    // clean up
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

