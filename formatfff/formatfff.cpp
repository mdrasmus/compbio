#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
    FILE *infile = stdin;
    char name[31];
    char species[31];
    char chrom[31];
    char startstr[21];
    int len;
    
    int nstarts = 1000;
    
    while (!feof(infile)) {
        fscanf(infile, "%30[^\t]\t%30[^\t]\t%30[^\t]\t%d\t", 
               name, species, chrom, &len);
        
        
        int i = nstarts;
        int pos = ftell(infile);
        while (fscanf(infile, "%20[0-9] ", startstr) == 1) {
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
