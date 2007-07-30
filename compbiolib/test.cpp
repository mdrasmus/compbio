#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

extern "C" {
#include "dnatrie.h"
};

static char *int2dna = "ACGT";

#define MAX_SEQ 30


char *randseq(char *seq, int size)
{
    for (int i=0; i<size; i++)
        seq[i] = int2dna[int(rand() / float(RAND_MAX) * 4)];
    seq[size] = '\0';
    return seq;
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    
    char seq[MAX_SEQ+1];
    
    DnaTrie *tree = NULL;
    
    dnatrie_insert(&tree, "accctg", 0);
    dnatrie_insert(&tree, "acccta", 0);
    dnatrie_insert(&tree, "acactg", 0);
    dnatrie_insert(&tree, "acactt", 0);
    
    printf("query = %d (6)\n", dnatrie_query(tree, "accctg"));
    printf("query = %d (4)\n", dnatrie_query(tree, "accc"));
    printf("query = %d (1)\n", dnatrie_query(tree, "atcc"));
    printf("query = %d (6)\n", dnatrie_query(tree, "acacttaaaa"));
    
    
    
    // insert a bunch of sequences
    for (int i=0; i<1e6; i++) {
        dnatrie_insert(&tree, randseq(seq, MAX_SEQ), 0);
    }
    
    
    
}
