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


#include "spidir.h"
#include "common.h"
#include "Tree.h"


#define DNA_A 0
#define DNA_C 1
#define DNA_G 2
#define DNA_T 3



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
    -1, -1, -1, -1, -1, -1
};


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


SpidirParams *readSpidirParams(const char* filename)
{
    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "cannot read file '%s'\n", filename);
        return NULL;
    }
    
    const int MAX_NAME = 51;
    float param1, param2;
    float alpha, beta;
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



template <class T>
void permute(T* array, int *perm, int size)
{
    ExtendArray<T> tmp(0, size);
    tmp.extend(array, size);
    
    // transfer permutation to temp array
    for (int i=0; i<size; i++)
        tmp[perm[i]] = array[i];
    
    // copy permutation back to original array
    for (int i=0; i<size; i++)
        array[i] = tmp[i];
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
                    perm.append(i);
                    break;
                }
            } else {
                if (atoi(names[j].c_str()) == inodes[i]) {
                    perm.append(i);
                    break;
                }
            }
        }
    }
    
    // apply permutation
    permute(names, perm, nsnodes);
    permute(mu, perm, nsnodes);
    permute(sigma, perm, nsnodes);
    
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


//=============================================================================
// Gene2species

const string Gene2species::NULL_SPECIES;

bool Gene2species::read(const char *filename)
{
    BufferedReader reader;
    if (!reader.open(filename, "r"))
        return false;

    char *line;
    string expr, species;
    char *ptr;
    while ((line = reader.readLine())) {
        //chomp(line);

        expr = strtok_r(line, "\t", &ptr);
        species = strtok_r(NULL, "\n", &ptr);

        if (expr[0] == '*') {
            // suffix
            m_rules.append(Gene2speciesRule(Gene2speciesRule::SUFFIX,
                                            expr.substr(1, expr.size()-1), 
                                            species));
        } else if (expr[expr.size() - 1] == '*') {
            // prefix
            m_rules.append(Gene2speciesRule(Gene2speciesRule::PREFIX,
                                            expr.substr(0, expr.size()-1), 
                                            species));
        } else {
            // exact match
            assert(0);
        }
    }

    return false;
}

string Gene2species::getSpecies(string gene)
{
    for (int i=0; i<m_rules.size(); i++) {
        switch (m_rules[i].rule) {
            case Gene2speciesRule::PREFIX:
                if (gene.find(m_rules[i].expr, 0) == 0)
                    return m_rules[i].species;
                break;

            case Gene2speciesRule::SUFFIX:
                if (gene.rfind(m_rules[i].expr, gene.size()-1) == 
                    gene.size() - m_rules[i].expr.size())
                    return m_rules[i].species;
                break;                

            case Gene2speciesRule::EXACT:
                break;
        }
    }

    return NULL_SPECIES;
}

bool Gene2species::getMap(string *genes, int ngenes, 
                          string *species, int nspecies, int *map)
{
    for (int i=0; i<ngenes; i++) {
        string sp = getSpecies(genes[i]);

        if (sp.size() == 0) {
            map[i] = -1;
        } else {
            map[i] = -1;
            for (int j=0; j<nspecies; j++) {
                if (sp == species[j])
                    map[i] = j;
            }
        }
    }

    return true;
}





/*

def makeGene2species(maps):
    # find exact matches and expressions
    exacts = {}
    exps = []
    for mapping in maps:
        if "*" not in mapping[0]:
            exacts[mapping[0]] = mapping[1]
        else:
            exps.append(mapping)
    
    # create mapping function
    def gene2species(gene):
        # eval expressions first in order of appearance
        for exp, species in exps:
            if exp[-1] == "*":
                if gene.startswith(exp[:-1]):
                    return species
            elif exp[0] == "*":
                if gene.endswith(exp[1:]):
                    return species
        
        if gene in exacts:
            return exacts[gene]
        
        raise Exception("Cannot map gene '%s' to any species" % gene)
    return gene2species


*/


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


bool chomp(char *str)
{
   int len = strlen(str);
   if (str[len-1] == '\n') {
      str[len-1] = '\0';
      return true;
   } else
      return false;
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


float readDist(FILE *infile, int &depth)
{
    float dist = 0;
    fscanf(infile, "%f", &dist);
    return dist;
}
