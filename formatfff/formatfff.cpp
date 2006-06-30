#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;


bool indexFff(string filename, int nstarts)
{
    FILE *infile;
    FILE *outfile;

    // parsing variables
    char name[51];
    char species[51];
    char chrom[51];
    char startstr[51];
    int len;
    char comment;
    
    
    // ensure correct filename names and construct output filename
    if (filename.substr(filename.size() - 4, 4) != ".fff") {
        fprintf(stderr, "formatfff: input file does not end with '.fff'\n");
        return false;
    }
    string outfilename = filename;
    outfilename[outfilename.size() - 1] = 'i';  // make new file a '*.ffi'
    
    
    // open input and output files
    if ((infile = fopen(filename.c_str(), "r")) == NULL) {
        fprintf(stderr, "formatfff: cannot open file '%s'\n", filename.c_str());
        return false;
    }
    
    if ((outfile = fopen(outfilename.c_str(), "w")) == NULL) {
        fprintf(stderr, "formatfff: cannot open output file '%s'\n", outfilename.c_str());
        return false;
    }
    
    
    // process file
    while (!feof(infile)) {
        // check for comments
        comment = fgetc(infile);
        
        if (comment == '\n') {
            // skip blank lines
            continue;
        } else if (comment == '#' || comment == ' ') {
            // skip comments (they start with '#' or ' ')
            while ((char) getc(infile) != '\n');
            continue;
        } else {
            // not a comment, replace char
            ungetc(comment, infile);
        }
            
        
        // parse frequent feature
        fscanf(infile, "%50[^\t]\t%50[^\t]\t%50[^\t]\t%d%*[\t\n]", 
               name, species, chrom, &len);
        
        // parse starts
        int i = nstarts;
        int pos = ftell(infile);
        int start2;
        int last = 0;
        while (fscanf(infile, "\t%50[-0-9]", startstr) == 1) {
            i++;
            
            // check for start sorting
            start2 = abs(atoi(startstr));
            if (start2 < last) {
                fprintf(stderr, "formatfff: feature starts are not sorted!\n");
                return false;
            }
            last = start2;
            
            
            if (i >= nstarts) {
                // output file position to index file
                int start = abs(atoi(startstr));
                fprintf(outfile, "%s\t%s\t%s\t%d\t%d\t%d\n", 
                        name, species, chrom, len, start, pos);
                i = 0;
            }
            
            // get file position
            pos = ftell(infile);
        }
    }
    
    // clean up
    fclose(infile);
    fclose(outfile);
    
    return true;
}


int main(int argc, char **argv)
{
    // number of feature per index line
    int nstarts = 1000;
    
    
    // check args
    if (argc < 2) {
        fprintf(stderr, "usage: formatfff [-n indexsize] <fff file>\n");
        exit(1);
    }
    
    // parse args
    int i;
    for (i=1; i<argc; i++) {
        if (argv[i] == "-n") {
            nstarts = atoi(argv[++i]);
        } else {
            // this must be a fff filename
            break;
        }
    }
    
    // process fff files
    for (; i<argc; i++) {
        if (!indexFff(argv[i], nstarts)) {
            fprintf(stderr, "formatfff: error processing '%s'\n", argv[i]);
            exit(1);
        }
        fprintf(stderr, "processed '%s'\n", argv[i]);
    }
    
    return 0;
}

