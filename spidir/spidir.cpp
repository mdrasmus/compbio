
#include "spidir.h"
#include "phylogeny.h"

namespace spidir {


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
// Distance Matrix

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


} // namespace spidir
