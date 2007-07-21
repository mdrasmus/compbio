/*=============================================================================

    SPIDIR - SPecies Informed DIstanced-based Reconstruction
    
    Matt Rasmussen
    Wed Jun 13 22:09:24 EDT 2007


=============================================================================*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "common.h"
#include "Tree.h"





int dna2int [256] = 
{
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 9
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 19
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 29
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 39
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 49
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 59
    -1, -1, -1, -1, -1,  0, -1,  1, -1, -1,   // 69
    -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,   // 79
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1,   // 89
    -1, -1, -1, -1, -1, -1, -1,  0, -1,  1,   // 99
    -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,   // 109
    -1, -1, -1, -1, -1, -1,  3, -1, -1, -1,   // 119
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 129
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 139
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 149
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 159
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 169
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 179
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 189
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 199
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 209
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 219
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 229
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 239
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,   // 249
    -1, -1, -1, -1, -1, -1                    // 255
};

char *int2dna = "ACGT";

int dnatype[] = { 
    DNA_PURINE,     // A
    DNA_PRYMIDINE,  // C
    DNA_PURINE,     // G
    DNA_PRYMIDINE   // T
};    


//=============================================================================
// Distance Matix

// calculate the pairwise distances between sequences
// NOTE: simple version implemented first
void calcDistMatrix(int nseqs, int seqlen, char **seqs, float **distmat)
{
    for (int i=0; i<nseqs; i++) {
        distmat[i][i] = 0.0;
    
        for (int j=i+1; j<nseqs; j++) {
            float changes = 0.0;
            int len = 0;
            
            for (int k=0; k<seqlen; k++) {
                if (seqs[i][k] != '-' && seqs[j][k] != '-') {
                    len++;                
                    if (seqs[i][k] != seqs[j][k]) {
                        changes += 1;
                    }
                }
            }
            
            assert(len > 0);
            
            distmat[i][j] = changes / len;
            distmat[j][i] = changes / len;
        }
    }
}





bool writeDistMatrix(const char *filename, int ngenes, float **dists, 
                     string *names)
{
    FILE *stream = NULL;
    
    if ((stream = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "cannot open '%s'\n", filename);
        return false;
    }
    
    // print number of genes
    fprintf(stream, "%d\n", ngenes);
    
    for (int i=0; i<ngenes; i++) {
        fprintf(stream, "%s ", names[i].c_str());
        
        for (int j=0; j<ngenes; j++) {
            fprintf(stream, "%f ", dists[i][j]);
        }
        fprintf(stream, "\n");
    }
    
    fclose(stream);
    return true;
}


//=============================================================================
// Spidir Parameters

SpidirParams *readSpidirParams(const char* filename)
{
    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "cannot read file '%s'\n", filename);
        return NULL;
    }
    
    const int MAX_NAME = 51;
    float param1, param2;
    float alpha = -1, beta = -1;
    ExtendArray<float> mu(0, 40);
    ExtendArray<float> sigma(0, 40);
    ExtendArray<string> names(0, 40);
    
    char name[MAX_NAME];
    
    while (!feof(infile)) {
        int ntokens = fscanf(infile, "%50s\t%f\t%f", name, &param1, &param2);
        if (ntokens <= 0)
            break;
        if (ntokens != 3) {
            return NULL;
        }
        
        if (!strcmp(name, "baserate")) {
            alpha = param1;
            beta = param2;
        } else {
            names.append(name);
            mu.append(param1);
            sigma.append(param2);
        }
    }
    fclose(infile);    
    
    return new SpidirParams(names.size(), names, mu, sigma, alpha, beta);
}


// get the preorder traversal of the species tree
void paramsOrder_helper(Node *node, ExtendArray<Node*> *nodeorder)
{
    nodeorder->append(node);
    for (int i=0; i<node->nchildren; i++) {
        paramsOrder_helper(node->children[i], nodeorder);
    }
}


// UNDER CONSTRUCTION
bool SpidirParams::order(SpeciesTree *stree)
{
    if (stree->nnodes != nsnodes) {
        printf("%d %d\n", stree->nnodes, nsnodes);
        return false;
    }
    
    ExtendArray<Node*> nodeorder(0, stree->nnodes);
    paramsOrder_helper(stree->root, &nodeorder);
    ExtendArray<int> perm(0, stree->nnodes);
    ExtendArray<int> invperm(0, stree->nnodes);
    
    // make interior node names
    ExtendArray<int> inodes(0, stree->nnodes);
    
    int inodename = 1;
    for (int i=0; i<stree->nnodes; i++) {
        if (nodeorder[i]->isLeaf()) {
            inodes.append(-1);
        } else {
            inodes.append(inodename++);
        }
    }
    
    
    // loop through preordered nodes to construct permutation
    for (int j=0; j<nsnodes; j++) {
        for (int i=0; i<stree->nnodes; i++) {
            if (nodeorder[i]->isLeaf()) {
                if (names[j] == nodeorder[i]->leafname) {
                    invperm.append(i);
                    break;
                }
            } else {
                if (atoi(names[j].c_str()) == inodes[i]) {
                    invperm.append(i);
                    break;
                }
            }
        }
    }
    
    // apply permutation
    invertPerm(invperm, perm, nsnodes);
    permute(names, perm, nsnodes);
    permute(mu, perm, nsnodes);
    permute(sigma, perm, nsnodes);
    
    return true;
}





//=============================================================================
// Math

// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s)
{
    if (s == 0.0)
        return -INFINITY;
    return - log(s) - log(sqrt(2.0*PI)) - (x-u)*(x-u) / (2.0*s*s);
    //return log(1.0/(s * sqrt(2.0*PI)) * exp(- (x-u)*(x-u) / (2.0 * s*s)));
}



/* gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009172947146,    -86.50532032941677,
                          24.01409824083091,    -1.231739572450155,
                          0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.5066282746310005*ser);
}



float gammalog(float x, float a, float b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}


// Invert a permutation
void invertPerm(int *perm, int *inv, int size)
{
    for (int i=0; i<size; i++)
        inv[perm[i]] = i;
}


//=============================================================================
// input/output

void printIntArray(int *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%d ", array[i]);
    printf("\n");
}

void printFloatArray(float *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%f ", array[i]);
    printf("\n");
}


bool inChars(char c, const char *chars)
{
   if (!chars)
      return false;
   for (;*chars; chars++)
      if (c == *chars) return true;
   return false;
}


bool chomp(char *str)
{
   int len = strlen(str);
   if (str[len-1] == '\n') {
      str[len-1] = '\0';
      return true;
   } else
      return false;
}


vector<string> split(const char *str, const char *delim, bool multiDelim)
{
    vector<string> tokens;   
    int i=0, j=0;
   
    while (str[i]) {
        // walk to end of next token
        for (; str[j] && !inChars(str[j], delim); j++);
        printf("> %d %d\n", i, j);

        if (i == j)
            break;

        // save token
        tokens.push_back(string(&str[i], j-i));
        
        if (!str[j])
            break;
        j++;
        i = j;
    }
    
    return tokens;
}


void printError(const char *fmt, ...)
{
   va_list ap;   
   va_start(ap, fmt);
   
   fprintf(stderr, "error: ");
   vfprintf(stderr, fmt, ap);
   fprintf(stderr, "\n");
   
   va_end(ap);
}


char readChar(FILE *stream, int &depth)
{
    char chr;
    do {
        if (fread(&chr, sizeof(char), 1, stream) != 1) {
            // indicate EOF
            return '\0';
        }
    } while (chr == ' ' || chr == '\n');
    
    // keep track of paren depth
    if (chr == '(') depth++;
    if (chr == ')') depth--;
    
    return chr;
}


char readUntil(FILE *stream, string &token, char *stops, int &depth)
{
    char chr;
    token = "";
    while (true) {
        chr = readChar(stream, depth);
        if (!chr)
            return chr;
        
        // compare char to stop characters
        for (char *i=stops; *i; i++) {
            if (chr == *i)
                return chr;
        }
        token += chr;
    }
}


string trim(const char *word)
{
    char buf[101];
    sscanf(word, "%100s", buf);
    return string(buf);
}

