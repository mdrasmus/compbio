/*=============================================================================

    SPIDIR    
    Maximum Likelihood Branch Length Estimation
    
    mldist.cpp
    started: Sun Jul  1 13:11:02 EDT 2007

=============================================================================*/

// c++ headers
#include <math.h>

// spidir headers
#include "common.h"
#include "Matrix.h"
#include "mldist.h"
#include "parsimony.h"
#include "spidir.h"
#include "Tree.h"



namespace spidir {

/*=============================================================================
    From: Felsenstein. Inferring Phylogenies. p 202.

    NOTE: HKY is a special case of Tamura-Nei where 
        alpha_r / alpha_y = rho = pi_r / pi_y 

    NOTE: parameters are chosen such that 
        P(j != i | i, t=1, pi, ratio) = 1
        thus, time also corresponds to sub/site

    Definitions:
        i    = destination base
        j    = source base
        t    = time
        pi_i = background/prior distribution of base i
        R    = Transtition/Transversion ratio

    Parameterization:
        beta = 1 / (2.0 * pi_r * pi_y * (1+R))
        rho = pi_r / pi_y    

        alpha_y = (pi_r * pi_y * R - pi_a*pi_g - pi_c*pi_t) / 
                  (2.0*(1+R)*(pi_y*pi_a*pi_g*rho + pi_r*pi_c*pi_t))
        alpha_r = rho * alpha_y    

    Convenience variables:
        if (dnatype[i] == DNA_PURINE) {
            alpha_i = alpha_r;
            pi_ry = pi_r;
        } else {
            alpha_i = alpha_y;
            pi_ry = pi_y;
        }
        int delta_ij =  int(i == j);
        int e_ij = int(dnatype[i] == dnatype[j]);    

    Formula:
        prob(j | i, t, pi, R) = 
            exp(-(alpha_i + beta)*t) * delta_ij + 
            exp(-beta*t) * (1 - exp(-alpha_i*t)) * (pi_j*e_ij/pi_ry) + 
            (1 - exp(-beta*t)) * pi_j

*/
class HkyModel
{
public:
    HkyModel(const float *bgfreq, float ratio) :
        ratio(ratio)
    {
        // set background base frequencies
        for (int i=0; i<4; i++)
            pi[i] = bgfreq[i];

        pi_r = pi[DNA_A] + pi[DNA_G];
        pi_y = pi[DNA_C] + pi[DNA_T];
        rho = pi_r / pi_y;


        // determine HKY parameters alpha_r, alpha_y, and beta
        b = 1.0 / (2.0 * pi_r * pi_y * (1.0+ratio));
        a_y = (pi_r*pi_y*ratio - 
               pi[DNA_A]*pi[DNA_G] - 
               pi[DNA_C]*pi[DNA_T]) / 
              (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + 
                              pi_r*pi[DNA_C]*pi[DNA_T]));
        a_r = rho * a_y;
    }


    // transition probability P(j | i, t)
    inline float operator()(int j, int i, float t)
    {
        swap(i, j);
        
        // convenience variables
        // NOTE: it is ok to assign pi_ry, because it is only used when
        // dnatype[i] == dnatype[j]
        float a_i, pi_ry;
        switch (dnatype[i]) {
            case DNA_PURINE:
                a_i = a_r;
                pi_ry = pi_r;  
                break;
            case DNA_PRYMIDINE:
                a_i = a_y;
                pi_ry = pi_y;
                break;
            default:
                assert(0);
        }
        int delta_ij = int(i == j);
        int e_ij = int(dnatype[i] == dnatype[j]);

        // return transition probability
        float ait = expf(-a_i*t);
        float ebt = expf(-b*t);
        
        float prob = ait*ebt * delta_ij + 
                     ebt * (1.0 - ait) * (pi[j]*e_ij/pi_ry) + 
                     (1 - ebt) * pi[j];
        return prob;        
    }
    
    // parameters
    float ratio;
    float pi[4];
    float pi_r;
    float pi_y;
    float rho;
    float b;
    float a_y;
    float a_r;
};


void makeHkyMatrix(const float *bgfreq, float ratio, float t, float *matrix)
{
    HkyModel model(bgfreq, ratio);
    
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            matrix[4*i+j] = model(i, j, t);
}




// This is a functor that represents the derivative of distLikelihood
// Finding the root of the derivative will find the maximum of distLikelihood
template <class Model>
class DistLikelihoodFunc
{
public:
    DistLikelihoodFunc(float *probs1, float *probs2, int seqlen, 
                       const float *bgfreq, Model *model, float step=.001,
                       int _sample=200, float *_probs3=NULL) :
        probs1(probs1),
        probs2(probs2),
        seqlen(seqlen),
        bgfreq(bgfreq),
        model(model),
        step(step),
        sample(_sample),
        bases(_sample),
        transmat1(16),
        transmat2(16),
        probs3(_probs3)
    {
        if (sample == 0) {
            sample = seqlen;
            bases.setCapacity(seqlen);
            
            for (int i=0; i<sample; i++) 
                bases[i] = i;
        } else {
            // choose a random sample of bases for estimating branch length
            for (int i=0; i<sample; i++) 
                bases[i] = irand(seqlen);
        }
        
        // allocate if needed the precomputed probability terms
        if (!probs3) {
            probs3 = new float [seqlen*4*4];
            ownProbs3 = true;
        } else
            ownProbs3 = false;
        
        // interate over sequence
        for (int kk=0; kk<sample; kk++) {
            int k = bases[kk];
            
            float terms1[4];
            float terms2[4];
            
            for (int i=0; i<4; i++)
                terms1[i] = bgfreq[i] * probs1[matind(4, k, i)];
            
            for (int j=0; j<4; j++)
                terms2[j] = probs2[matind(4, k, j)];
            
            // integrate over all transitions
            float *ptr = &probs3[matind(16, k, 0)];
            for (int ij=0; ij<16; ij++) {
                int i = ij >> 2;
                int j = ij & 3;
                ptr[ij] = terms1[i] * terms2[j];
            }
        }
    }
    
    ~DistLikelihoodFunc()
    {
        // free probability table if we own it
        if (ownProbs3)
            delete [] probs3;
    }


    float operator()(float t)
    {
        // trivial case
        if (t < 0)
            return INFINITY;

        float totprob = 1.0;

        
        // build transition matrices
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                transmat1[4*i+j] = (*model)(i, j, t+step);
                transmat2[4*i+j] = (*model)(i, j, t);
            }
        }
        
        // iterate over sites
        for (int kk=0; kk<sample; kk++) {
            int k = bases[kk];
            float prob1 = 0.0, prob2 = 0.0;

            // integrate over all transitions
            float *termptr = &probs3[matind(16, k, 0)];
            for (int ij=0; ij<16; ij++) {
                prob1 += termptr[ij] * transmat1[ij];
                prob2 += termptr[ij] * transmat2[ij];
            }
            totprob *= prob1 / prob2;
        }
        
        return logf(totprob) / step;
    }
    
    float *probs1;
    float *probs2;
    int seqlen;
    const float *bgfreq;
    Model *model;
    float step;
    int sample;
    ExtendArray<int> bases;
    ExtendArray<float> transmat1;
    ExtendArray<float> transmat2;
    float *probs3;
    bool ownProbs3;
};


// Find the maximum likelihood estimate (MLE) of the distance (time) between 
// two sequences (represented probabilistically as probs1 and probs2)
//
//  bgfreq = background frequency of bases
//  model = sequence evolution model
//
template <class Model>
float mleDistance(float *probs1, float *probs2, int seqlen, 
                  const float *bgfreq, Model &model, 
                  float t0=.001, float t1=1, float step=.0001,
                  float *probs3=NULL, int samples=500)
{
    DistLikelihoodFunc<Model> df(probs1, probs2, seqlen, bgfreq, &model, step,
                                 samples, probs3);
    
    const int maxiter = 20;    
    return bisectRoot(df, t0, t1, maxiter);
}


//=============================================================================
// Dynamic programming of conditional likelihood


/*

    From: Inferring Phylogenies. Felsenstein. p 254.
    
    Baisc recursion of conditional likelihoods
        L_k^(i)(s) = (sum_x P(x|s,t_l) L_l^(i)(x)) (sum_y P(y,s,t_m) L_m^(i)(y))
    
    where k is a node
          l and m are children of k
          i indexes sites
          t_l and t_m are the branch lengths of l and m
    
    I define my table to be:
        lktable[i][matind(4, j, k)] = log(L_i^(j)(k))

*/


// conditional likelihood recurrence
template <class Model>
inline void calcLkTable(ExtendArray<float*> &lktable, int seqlen, Model &model,
                        int a, int b, int c, float adist, float bdist)
{
    static Matrix<float> atransmat(4, 4);
    static Matrix<float> btransmat(4, 4);
    
    // build transition matrices
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            atransmat[i][j] = model(i, j, adist);
            btransmat[i][j] = model(i, j, bdist);
        }
    }
    
    float *lktablea = lktable[a];
    float *lktableb = lktable[b];
    float *lktablec = lktable[c];
    
    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        float terma[4];
        float termb[4];       
        
        for (int x=0; x<4; x++) {
            terma[x] = lktablea[matind(4, j, x)]; //expf(lktablea[matind(4, j, x)]);
            termb[x] = lktableb[matind(4, j, x)]; //expf(lktableb[matind(4, j, x)]);
        }
    
        
        for (int k=0; k<4; k++) {
            //float prob1 = 0.0, prob2 = 0.0;                

            //for (int x=0; x<4; x++) {
            float *aptr = atransmat[k];
            float *bptr = btransmat[k];
            
            float prob1 = aptr[0] * terma[0] +
                          aptr[1] * terma[1] +
                          aptr[2] * terma[2] +
                          aptr[3] * terma[3];
            float prob2 = bptr[0] * termb[0] +
                          bptr[1] * termb[1] +
                          bptr[2] * termb[2] +
                          bptr[3] * termb[3];
            

            lktablec[matind(4, j, k)] = prob1 * prob2; //logf(prob1 * prob2);
        }
        
        /*
        for (int k=0; k<4; k++) {
            float prob1 = 0.0, prob2 = 0.0;                

            for (int x=0; x<4; x++) {                
                prob1 += atransmat[k][x] * terma[x];
                prob2 += btransmat[k][x] * termb[x];
            }

            lktable[c][matind(4, j, k)] = logf(prob1 * prob2);
        }*/
        
    }
}


// initialize the condition likelihood table
template <class Model>
void initCondLkTable(ExtendArray<float*> &lktable, Tree *tree, 
                     int nseqs, int seqlen, char **seqs, Model &model)
{
    // recursively calculate cond. lk. of internal nodes
    ExtendArray<Node*> nodes;
    getTreePostOrder(tree, &nodes);
    
    for (int l=0; l<nodes.size(); l++) {
        Node *node = nodes[l];
        int i = node->name;
        
        if (node->isLeaf()) {
            // initialize leaves from sequence
        
            // iterate over sites
            for (int j=0; j<seqlen; j++) {
                int base = dna2int[int(seqs[i][j])];

                if (base == -1) {
                    // handle gaps
                    lktable[i][matind(4, j, 0)] = 1.0; // 0.0;
                    lktable[i][matind(4, j, 1)] = 1.0; // 0.0;
                    lktable[i][matind(4, j, 2)] = 1.0; // 0.0;
                    lktable[i][matind(4, j, 3)] = 1.0; // 0.0;
                } else 
                    lktable[i][matind(4, j, base)] = 1.0; // 0.0;
            }
        } else {
            // compute internal nodes from children
            Node *node1 = node->children[0];
            Node *node2 = node->children[1];
            
            calcLkTable(lktable, seqlen, model, node1->name, node2->name, i,
                        node1->dist, node2->dist);
        }
    }
}


// calculate P(D | T, B)
template <class Model>
float getTotalLikelihood(ExtendArray<float*> &lktable, Tree *tree, 
                         int seqlen, Model &model, const float *bgfreq)
{
    // find total likelihood
    calcLkTable(lktable, seqlen, model, 
                tree->root->children[0]->name, 
                tree->root->children[1]->name, 
                tree->root->name,
                tree->root->children[0]->dist, 
                tree->root->children[1]->dist);
    
    float lk = 0.0;
    for (int k=0; k<seqlen; k++) {
        float prob = 0.0;
        for (int x=0; x<4; x++) {
            prob += bgfreq[x] * 
                    lktable[tree->root->name][matind(4, k, x)];
                    //expf(lktable[tree->root->name][matind(4, k, x)]);
        }
        lk += logf(prob);
    }
    return lk;
}



//=============================================================================
// find MLE branch lengths


// TODO: need to think about more carefully
void getRootOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node=NULL)
{
    if (!node) {
        // start at tree root (but don't include root)
        getRootOrder(tree, nodes, tree->root->children[0]);
        getRootOrder(tree, nodes, tree->root->children[1]);
    } else {
        // record pre-process
        if (node->parent == tree->root) {
            nodes->append(tree->root->children[0]);
            nodes->append(tree->root->children[1]);
        } else {
            nodes->append(node);
            nodes->append(node->parent);
        }

        // recurse
        for (int i=0; i<node->nchildren; i++)
            getRootOrder(tree, nodes, node->children[i]);
        
        /*
        // record post-process
        if (!node->isLeaf()) {
            if (node->parent == tree->root) {
                nodes->append(tree->root->children[0]);
                nodes->append(tree->root->children[1]);
            } else {
                nodes->append(node);
                nodes->append(node->parent);
            }
        }*/
    }
}

template <class Model>
float calcSeqProb(Tree *tree, int nseqs, char **seqs, 
                  const float *bgfreq, Model &model)
{
    int seqlen = strlen(seqs[0]);
    
    // allocate conditional likelihood dynamic programming table
    ExtendArray<float*> lktable(tree->nnodes);
    for (int i=0; i<tree->nnodes; i++) {
        lktable[i] = new float [4 * seqlen];
        for (int j=0; j<4*seqlen; j++)
            lktable[i][j] = 0.0; //-INFINITY;
    }
    
    // initialize the condition likelihood table
    initCondLkTable(lktable, tree, nseqs, seqlen, seqs, model);
    float logl = getTotalLikelihood(lktable, tree, seqlen, model, bgfreq);
    
    // cleanup
    for (int i=0; i<tree->nnodes; i++)
        delete [] lktable[i];
    
    return logl;
}


float calcHkySeqProb(Tree *tree, int nseqs, char **seqs, 
                     const float *bgfreq, float ratio)
{
    HkyModel hky(bgfreq, ratio);
    return calcSeqProb(tree, nseqs, seqs, bgfreq, hky);
}


// NOTE: assumes binary Tree
template <class Model>
float findMLBranchLengths(Tree *tree, int nseqs, char **seqs, 
                          const float *bgfreq, Model &model,
                          int maxiter=10)
{
    const float converge = logf(2.0);
    const int samples = 0; // fixed for now  
    //int convergenum = 1000; //2* tree->nnodes;

    int seqlen = strlen(seqs[0]);
    float lastLogl = -INFINITY, logl = -INFINITY;
    
    
    // allocate conditional likelihood dynamic programming table
    ExtendArray<float*> lktable(tree->nnodes);
    for (int i=0; i<tree->nnodes; i++) {
        lktable[i] = new float [4 * seqlen];
        for (int j=0; j<4*seqlen; j++)
            lktable[i][j] = 0.0; //-INFINITY;
    }
    
    // allocate auxiliary likelihood table
    ExtendArray<float> probs3(seqlen*4*4);
    
    // initialize the condition likelihood table
    initCondLkTable(lktable, tree, nseqs, seqlen, seqs, model);
    
    // remember original rooting for restoring later
    Node *origroot1 = tree->root->children[0];
    Node *origroot2 = tree->root->children[1];
    
    
    // determine rooting order
    ExtendArray<Node*> rootingOrder(0, 2*tree->nnodes);
    getRootOrder(tree, &rootingOrder);
    //getRootOrder(tree, &rootingOrder);
    //getRootOrder(tree, &rootingOrder);
    
    // add additional random roots
    //for (int i=rootingOrder.size(); i<maxiter; i++) {
    //    rootingOrder.append(tree->nodes[irand(tree->nnodes)]);
    //}
    
    
    // iterate over branches improving each likelihood
    for (int j=0; j<maxiter; j++) {
        printLog(LOG_HIGH, "hky: iter %d\n", j);    
        for (int i=0; i<rootingOrder.size(); i+=2) {
            // remembering old children of root
            Node *oldnode1 = tree->root->children[0];
            Node *oldnode2 = tree->root->children[1];

            // choose new root
            tree->reroot(rootingOrder[i], rootingOrder[i+1]);
            if (tree->root->children[0] == oldnode1 &&
                tree->root->children[1] == oldnode2)
                continue;

            // rebuild invaild likelihood values
            // determine starting node to rebuild
            Node *ptr;
            if (oldnode1->parent == oldnode2)
                ptr = oldnode1;
            else if (oldnode2->parent == oldnode1)
                ptr = oldnode2;
            else
                assert(0);

            // walk up to root of tree, rebuilding conditional likelihoods
            for (; ptr->parent; ptr = ptr->parent) {
                if (!ptr->isLeaf())
                    calcLkTable(lktable, seqlen, model, 
                                ptr->children[0]->name, 
                                ptr->children[1]->name, 
                                ptr->name,
                                ptr->children[0]->dist, 
                                ptr->children[1]->dist);
            }

            // get total probability before branch length change
            float loglBefore = getTotalLikelihood(lktable, tree, 
                                                  seqlen, model, bgfreq);
            logl = -INFINITY;
            Node *node1 = tree->root->children[0];
            Node *node2 = tree->root->children[1];

            // find new MLE branch length for root branch
            float initdist = node1->dist + node2->dist;
            float mle = mleDistance(lktable[node1->name], 
                                    lktable[node2->name], seqlen, 
                                    bgfreq, model, 
                                    max(initdist*.95, 0.0), 
                                    max(initdist*1.05, 0.001),
                                    .001, probs3, samples);
            node1->dist = mle / 2.0;
            node2->dist = mle / 2.0;

            // get total probability after bracnh change
            logl = getTotalLikelihood(lktable, tree, 
                                      seqlen, model, bgfreq);

            // don't accept a new branch length if it lowers total likelihood
            if (logl < loglBefore) {
                // revert
                node1->dist = initdist / 2.0;
                node2->dist = initdist / 2.0;
                logl = loglBefore;
            }
            
            printLog(LOG_HIGH, "hky: lk=%f\n", logl);        
        }
        
        // determine whether logl has converged
        float diff = fabs(logl - lastLogl);
        if (diff < converge) {
            printLog(LOG_HIGH, "hky: diff = %f < %f\n", diff, converge);
            break;
        } else {
            printLog(LOG_HIGH, "hky: diff = %f > %f\n", diff, converge);
        }
        lastLogl = logl;
    }
    
    // restore original rooting
    if (origroot1->parent != tree->root ||
        origroot2->parent != tree->root)
    {
        if (origroot1->parent == origroot2) {
            tree->reroot(origroot1);
        } else if  (origroot2->parent == origroot1) {
            tree->reroot(origroot2);
        } else
            assert(0);
    }
    
    
    // cleanup
    for (int i=0; i<tree->nnodes; i++)
        delete [] lktable[i];
    
    return logl;
}


float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                             const float *bgfreq, float ratio, int maxiter)
{
    HkyModel hky(bgfreq, ratio);
    return findMLBranchLengths(tree, nseqs, seqs, bgfreq, hky, maxiter);
}


float findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                             float *dists, const float *bgfreq, float ratio, 
                             int maxiter, bool parsinit)
{
    //int seqlen = strlen(seqs[0]);
        
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    
    if (parsinit)
        parsimony(&tree, nseqs, seqs);
    
    float logl = findMLBranchLengthsHky(&tree, nseqs, seqs, bgfreq, 
                                        ratio, maxiter);
    tree.getDists(dists);
    
    return logl;
}


} // namespace spidir
