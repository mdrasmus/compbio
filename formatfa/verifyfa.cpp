/********************************************************************************
    verifyfa
    Matt Rasmussen    
    Fri Jul 14 21:53:30 EDT 2006



    Creates convenient index files for fasta files
********************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

#define MAXLINE 1000


// TODO:  verify that the fasta file has consistent wrapping


inline int chomp(char *line)
{
    int len = strlen(line);
    if (line[len - 1] == '\n') {
        line[len - 1] = '\0';
        return len - 1;
    } else {
        return len;
    }
}


bool verifyFasta(string filename)
{
    FILE *infile;
    
    // open input file
    if ((infile = fopen(filename.c_str(), "r")) == NULL) {
        fprintf(stderr, "  cannot open file '%s'\n", filename.c_str());
        return false;
    }
    
    
    int lineno = 0;
    int width = 0;
    int diffwidth = 0;
    bool lastline = false;
    char *line = new char [MAXLINE];
    bool ret = true;
    
    // process file
    while (fgets(line, MAXLINE, infile)) {
        lineno++;
        
        if (strlen(line) == MAXLINE - 1) {
            fprintf(stderr, "line %d too long\n", lineno);
            ret = false;
            break;
        }
        
        if (line[0] == '>') {
            width = -1;
            lastline = false;
        } else {
            int len = chomp(line);
            
            if (lastline) {
                fprintf(stderr, "line %d occurs after shorter line\n", 
                        lineno, width, len);
                ret = false;
                break;
            }
            
            if (width == -1) {
                // width is unset, set it
                width = len;
            } else if (len < width) {
                lastline = true;
            } else if (len > width) {
                fprintf(stderr, "line %d has width %d != %d\n", 
                        lineno, len, width);
                ret = false;
                break;
            }
        }
    }
    
    // clean up
    delete [] line;
    fclose(infile);
    
    return ret;
}


int main(int argc, char **argv)
{
    // check args
    if (argc < 2) {
        fprintf(stderr, "usage: verifyfa <fasta file> ...\n");
        exit(1);
    }
    
    // return code
    int ret = 0;
    
    // process fff files
    for (int i=1; i<argc; i++) {
        fprintf(stderr, "%s\t", argv[i]);
        if (!verifyFasta(argv[i])) {
            //fprintf(stderr, "ERROR\n", argv[i]);
            ret = 1;
        } else {
            fprintf(stderr, "OK\n", argv[i]);
        }
    }
    
    return ret;
}

