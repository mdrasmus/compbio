#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
    FILE *infile;
    char name[51];
    char species[51];
    char chrom[51];
    char startstr[51];
    int len;
    char comment;
    
    // number of feature per index line
    int nstarts = 1000;
    
    // check args
    if (argc < 2) {
        fprintf(stderr, "usage: formatfff <fff file>\n");
        exit(1);
    }
    
    // open input file
    if ((infile = fopen(argv[1], "r")) == NULL) {
        fprintf(stderr, "cannot open file '%s'\n", argv[1]);
        exit(1);
    }
    
    
    
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
        
        int i = nstarts;
        int pos = ftell(infile);
        while (fscanf(infile, " %50[-0-9]", startstr) == 1) {
            i++;
            
            if (i >= nstarts) {
                int start = abs(atoi(startstr));
                printf("%s\t%s\t%s\t%d\t%d\t%d\n", 
                       name, species, chrom, len, start, pos);
                i = 0;
            }
            
            pos = ftell(infile);
        }
    }
}

