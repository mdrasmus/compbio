/*=============================================================================

    SPIDIR    
    Branch Length Estimation
    
    branchlen.cpp
    started: Sun Jul  1 13:11:02 EDT 2007

=============================================================================*/

// c++ headers
#include <math.h>

// spidir headers
#include "parsimony.h"
#include "common.h"
#include "Matrix.h"
#include "Tree.h"
#include "branchlen.h"
#include "spidir.h"




namespace spidir {

/*=============================================================================
    From: Felsenstein. Inferring Phylogenies. p 202.

    NOTE: HKY is a special case of Tamura-Nei where 
        alpha_r / alpha_y = rho = pi_r / pi_y 

    NOTE: parameters are chosen such that 
        P(i != j | j, t=1, pi, ratio) = 1
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

    prob(i | j, t, pi, R) = 
        exp(-(alpha_i + beta)*t) * delta_ij + 
        exp(-beta*t) * (1 - exp(-alpha_i*t)) * (pi_j*e_ij/pi_ry) + 
        (1 - exp(-beta*t)) * pi_j

*/
class HkyModel
{
public:
    HkyModel(float *bgfreq, float ratio) :
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
        a_y = (pi_r*pi_y*ratio - pi[DNA_A]*pi[DNA_G] - pi[DNA_C]*pi[DNA_T]) / 
          (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + pi_r*pi[DNA_C]*pi[DNA_T]));
        a_r = rho * a_y;

    }


    // transition probability P(i | j, t)
    inline float operator()(int i, int j, float t)
    {
        //  convenience variables
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
                     ebt * (1 - ait) * (pi[j]*e_ij/pi_ry) + 
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




template <class Model>
float distLikelihood(float *probs1, float *probs2, 
                      int seqlen, Model &model, float *bgfreq, float t)
{
    // trivial case
    if (t < 0)
        return 0;
    
    float logl = 0;
    float logbgfreq[4];
    logbgfreq[0] = log(bgfreq[0]);
    logbgfreq[1] = log(bgfreq[1]);
    logbgfreq[2] = log(bgfreq[2]);
    logbgfreq[3] = log(bgfreq[3]);
    
    // interate over sequence
    for (int k=0; k<seqlen; k++) {
        float prob = 0.0;
        
        // integrate over all transitions
        for (int i=0; i<4; i++) {
            float logterm = logbgfreq[i] + probs1[matind(4, k, i)];
            for (int j=0; j<4; j++)
                prob += expf(logterm + probs2[matind(4, k, j)]) * model(i, j, t);
        }
        logl += logf(prob);
    }
    
    return logl;
}


// Find a root of a function func(x) using the secant method
// x0 and x1 are initial estimates of the root
template <class Func>
float secantRoot(Func &f, float x0, float x1, int maxiter, 
                 float minx=.000001, float esp=.002)
{
    float f0 = f(x0);
    for (int i=0; i<maxiter; i++) {
        if (fabs((x1 - x0)*2.0 / (x0+x1)) < esp)
            return x0;
        float f1 = f(x1);
        float x2 = x1 - (x1 - x0) * f1 / (f1 - f0);
        //printf("%f %f | %f %f\n", x0, x1, f0, f1);
        
        x0 = x1;
        x1 = (x2 > minx) ? x2 : minx;
        f0 = f1;
    }

    return x1;
}


// Find a root of a function func(x) using the bisection method
// This is less efficient but is more robust than Newton's or Secant
// x0 and x1 are initial estimates of the root
template <class Func>
float bisectRoot(Func &f, float x0, float x1, int maxiter, 
                 float minx=.0001, float maxx=10.0, float esp=.02)
{
    // we expect f(x0) > 0 and f(x1) < 0
    
    // move x0 left until f(x0) > 0
    float f0 = f(x0);
    while (f0 < 0) {
        //printf("x0=%f f0=%f\n", x0, f0);
        x0 /= 2.0;
        f0 = f(x0);
        if (x0 < minx)
            return x0;
    }
    
    // move x1 right until f(x1) < 0
    float f1 = f(x1);
    while (f1 > 0) {
        //printf("x1=%f f1=%f\n", x1, f1);
        x1 *= 2.0;
        f1 = f(x1);
        if (x1 > maxx)
            return x1;
    }
    
    for (int i=0; i<maxiter; i++) {
        if (fabs((x1 - x0)*2.0 / (x0+x1)) < esp)
            return x0;
        
        float x2 = (x0 + x1) / 2.0;
        float f2 = f(x2);
        
        if (f2 > 0) {
            x0 = x2;
            f0 = f2;
        } else {
            x1 = x2;
            f1 = f2;
        }
        
        //printf("%d: %.4f %.4f | %.2f %.2f\n", i, x0, x1, f0, f1);
    }

    return x1;
}



// This is a functor that represents the derivative of distLikelihood
// Finding the root of the derivative will find the maximum of distLikelihood
template <class Model>
class DistLikelihoodFunc
{
public:
    DistLikelihoodFunc(float *probs1, float *probs2, int seqlen, 
                       float *bgfreq, Model *model, float step=.001,
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
                bases[i] = int(frand() * seqlen);
        }
        
        
        // save log back ground frequency
        logbgfreq[0] = log(bgfreq[0]);
        logbgfreq[1] = log(bgfreq[1]);
        logbgfreq[2] = log(bgfreq[2]);
        logbgfreq[3] = log(bgfreq[3]);        
        
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
                terms1[i] = expf(logbgfreq[i] + probs1[matind(4, k, i)]);
            
            for (int j=0; j<4; j++)
                terms2[j] = expf(probs2[matind(4, k, j)]);
            
            float *ptr = &probs3[matind(16, k, 0)];
            for (int ij=0; ij<16; ij++) {
                int i = ij >> 2;
                int j = ij & 3;
                ptr[ij] = terms1[i] * terms2[j];
            }
            
            
            /*
            // integrate over all transitions
            for (int ij=0; ij<16; ij++) {
                int i = ij / 4;
                int j = ij % 4;
                float logterm = logbgfreq[i] + probs1[matind(4, k, i)];
                probs3[matind(16, k, i*4 + j)] = 
                        expf(logterm + probs2[matind(4, k, j)]);
            }
            */
            
            /*
            // integrate over all transitions
            for (int i=0; i<4; i++) {
                float logterm = logbgfreq[i] + probs1[matind(4, k, i)];
                for (int j=0; j<4; j++) {
                    probs3[matind(16, k, i*4 + j)] = 
                        expf(logterm + probs2[matind(4, k, j)]);
                }
            }
            */
        }
    }
    
    ~DistLikelihoodFunc()
    {
        if (ownProbs3)
            delete [] probs3;
    }


    float operator()(float t)
    {
        // trivial case
        if (t < 0)
            return INFINITY;

        float logl = 0;

        
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
                float term2 = termptr[ij];
                prob1 += term2 * transmat1[ij];
                prob2 += term2 * transmat2[ij];
            }
            logl += logf(prob1 / prob2);
        }
        
        /*
        // interate over sequence
        for (int kk=0; kk<sample; kk++) {
            int k = bases[kk];
            float prob1 = 0.0, prob2 = 0.0;

            // integrate over all transitions
            for (int i=0; i<4; i++) {
                float logterm = logbgfreq[i] + probs1[matind(4, k, i)];
                for (int j=0; j<4; j++) {
                    float term2 = expf(logterm + probs2[matind(4, k, j)]);
                    prob1 += term2 * transmat1[i][j];
                    prob2 += term2 * transmat2[i][j];
                }
            }
            logl += logf(prob1 / prob2);
        }
        */
        
        return logl / step;
    }
    
    float *probs1;
    float *probs2;
    int seqlen;
    float *bgfreq;
    Model *model;
    float step;
    int sample;
    ExtendArray<int> bases;
    ExtendArray<float> transmat1;
    ExtendArray<float> transmat2;
    float *probs3;
    bool ownProbs3;
    float logbgfreq[4];
};


// Find the maximum likelihood estimate (MLE) of the distance (time) between 
// two sequences (represented probabilistically as probs1 and probs2)
//
//  bgfreq = background frequency of bases
//  model = sequence evolution model
//
template <class Model>
float mleDistance(float *probs1, float *probs2, int seqlen, 
                  float *bgfreq, Model &model, 
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
    
    // iterate over sites
    for (int j=0; j<seqlen; j++) {
        float terma[4];
        float termb[4];
        
        for (int x=0; x<4; x++) {
            terma[x] = expf(lktablea[matind(4, j, x)]);
            termb[x] = expf(lktableb[matind(4, j, x)]);
        }
    
        for (int k=0; k<4; k++) {
            float prob1 = 0.0, prob2 = 0.0;                

            for (int x=0; x<4; x++) {                
                prob1 += atransmat[k][x] * terma[x];
                prob2 += btransmat[k][x] * termb[x];
            }

            lktable[c][matind(4, j, k)] = logf(prob1 * prob2);
        }
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
                    lktable[i][matind(4, j, 0)] = 0.0;
                    lktable[i][matind(4, j, 1)] = 0.0;
                    lktable[i][matind(4, j, 2)] = 0.0;
                    lktable[i][matind(4, j, 3)] = 0.0;
                } else 
                    lktable[i][matind(4, j, base)] = 0.0;
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
                         int seqlen, Model &model, float *bgfreq)
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
                    expf(lktable[tree->root->name][matind(4, k, x)]);
        }
        lk += logf(prob);
    }
    return lk;
}



//=============================================================================
// find MLE branch lengths

// NOTE: assumes binary Tree
template <class Model>
float findMLBranchLengths(Tree *tree, int nseqs, char **seqs, 
                         float *bgfreq, Model &model,
                         int maxiter=10)
{
    int seqlen = strlen(seqs[0]);
    float lastLogl = -INFINITY, logl = -INFINITY;
    float converge = logf(2.0);
    
    
    // allocate conditional likelihood dynamic programming table
    ExtendArray<float*> lktable(tree->nnodes);
    for (int i=0; i<tree->nnodes; i++) {
        lktable[i] = new float [4 * seqlen];
        for (int j=0; j<4*seqlen; j++)
            lktable[i][j] = -INFINITY;
    }
    
    // allocate auxiliary likelihood table
    ExtendArray<float> probs3(seqlen*4*4);
    
    // initialize the condition likelihood table
    initCondLkTable(lktable, tree, nseqs, seqlen, seqs, model);
    
    Node *origroot1 = tree->root->children[0];
    Node *origroot2 = tree->root->children[1];
    int convergenum = 10;
    const int samples = 0; // fixed for now
    
    // iterate over branches improving each likelihood
    for (int i=0; i<maxiter; i++) {
        printLog("hky: iter %d\n", i);
        
        Node *newroot = tree->nodes[int(frand() * tree->nnodes)];
        
        // remembering old children of root
        Node *oldnode1 = tree->root->children[0];
        Node *oldnode2 = tree->root->children[1];

        // reroot tree, skip if choosen node does not change root
        if (newroot == tree->root ||
            newroot == oldnode1 ||
            newroot == oldnode2)
            continue;
        tree->reroot(newroot);


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
                                max(initdist, (float)0.0), 
                                max(initdist*1.1, 0.001),
                                .001, probs3, samples);
        node1->dist = mle / 2.0;
        node2->dist = mle / 2.0;

        // get total probability after bracnh change
        logl = getTotalLikelihood(lktable, tree, 
                                  seqlen, model, bgfreq);

        if (logl < loglBefore) {
            // revert
            node1->dist = initdist / 2.0;
            node2->dist = initdist / 2.0;
            logl = loglBefore;
        }

        printLog("hky: lk=%f\n", logl);
        
        // determine whether logl has converged
        float diff = fabs(logl - lastLogl);
        if (i > 0 && diff < converge) {
            printLog("hky: diff = %f < %f\n", diff, converge);
            convergenum--;
            if (convergenum < 0)
                i = maxiter;
        } else {
            printLog("hky: diff = %f > %f\n", diff, converge);
        }
        lastLogl = logl;
    }
    
    // restore original rooting
    if (origroot1->parent != tree->root) {
        if (origroot1->parent == origroot2)
            tree->reroot(origroot1);
        else if  (origroot2->parent == origroot1)
            tree->reroot(origroot2);
        else
            assert(0);
    }
    
    
    // cleanup
    for (int i=0; i<tree->nnodes; i++)
        delete [] lktable[i];
    
    return logl;
}


float findMLBranchLengthsHky(Tree *tree, int nseqs, char **seqs, 
                            float *bgfreq, float ratio, int maxiter)
{
    HkyModel hky(bgfreq, ratio);
    return findMLBranchLengths(tree, nseqs, seqs, bgfreq, hky, maxiter);
}


float findMLBranchLengthsHky(int nnodes, int *ptree, int nseqs, char **seqs, 
                          float *dists, float *bgfreq, float ratio, int maxiter,
                          bool parsinit)
{
    int seqlen = strlen(seqs[0]);
        
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
