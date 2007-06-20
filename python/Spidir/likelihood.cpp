//=============================================================================
//  SPIDIR - Likelihood calculation


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "spidir.h"
#include "common.h"



// Branch distribution parameters for one branch
class BranchParam
{
public:
    BranchParam(float mu=-1.0, float sigma=-1.0) :
        mu(mu),
        sigma(sigma)
    {}
    
    bool isNull()
    {
        return mu == -1.0;
    }
    
    float mu;
    float sigma;
};

BranchParam NULL_PARAM;


// Reconciliation parameters
class ReconParams
{
public:
    ReconParams(int nnodes, bool *freebranches, int unfold, float unfolddist) :
        nnodes(nnodes),
        freebranches(freebranches),
        unfold(unfold),
        unfolddist(unfolddist)
    {
        startparams = new BranchParam [nnodes];
        midparams = new BranchParam [nnodes];
        endparams = new BranchParam [nnodes];
        
        startfrac = new int [nnodes];
        endfrac = new int [nnodes];
        midpoints = new float [nnodes];
    }
    
    ~ReconParams()
    {
        delete [] startparams;
        delete [] midparams;
        delete [] endparams;
        
        delete [] startfrac;
        delete [] endfrac;
        delete [] midpoints;
    }
    
    
    int nnodes;
    BranchParam *startparams;
    BranchParam * midparams;
    BranchParam *endparams;
    int *startfrac;
    int *endfrac;
    float *midpoints;
    bool *freebranches;
    int unfold;
    float unfolddist;
};


// solves x^3 + ax^2 + bx + x = 0 for x
//   in our case there is only one positive root; find it
float maxCubicRoot(float a, float b, float c)
{
    float const accuracy = .001;
    float x = 0.0;
    float y;
    
    // first estimate should be negative    
    assert (x*x*x + a*x*x + b*x + c <= 0);
    
    // increase x until y is positive
    x = .01;
    do {
        x *= 2.0;
        y = x*x*x + a*x*x + b*x + c;
    } while (y < 0);
    
    // binary search to find root
    float min = x / 2.0;
    float max = x;
    
    while (max - min > accuracy) {
        x = (max + min) / 2.0;
        y = x*x*x + a*x*x + b*x + c;
        
        if (y == 0)
            return x;
        else if (y > 0)
            max = x;
        else
            min = x;
    }
    
    return x;
}


float mleGenerate(int count, float *dists, float *means, float *sdevs,
                  float alpha, float beta)
{
    float a, b, c;
    
    a = (1.0 - alpha) / beta;

    // b = sum(means[i] * lens[i] / sdevs[i]**2
    //         for i in range(len(lens))) / beta    
    // c = - sum(lens[i] ** 2 / sdevs[i] ** 2
    //          for i in range(len(lens))) / beta    
    
    b = 0.0;
    c = 0.0;    
    for (int i=0; i<count; i++) {
        b += means[i] * dists[i] / (sdevs[i] * sdevs[i]);
        c += dists[i] * dists[i] / (sdevs[i] * sdevs[i]);
    }
    b /= beta;
    c = -c / beta;
    
    return maxCubicRoot(a, b, c);
}


// help set depths for every node
// depth is distance to next subtree root
void estimateGenerate_helper(Tree *tree, Node *node, float *depths, int *sroots,
                             int *recon, int *events, bool free)
{
    int name = node->name;
    
    if (events[name] == EVENT_SPEC)
        free = false;
    
    if (node != tree->root) {
        int parent = node->parent->name;
        
        if (free) {
            // mark freebranches with -1
            depths[name] = -1;
            sroots[name] = sroots[parent];
        } else {
            if (events[parent] == EVENT_DUP) {
                depths[name] = depths[parent] + node->dist;
                sroots[name] = sroots[parent];
            } else {
                depths[name] = node->dist;
                sroots[name] = recon[parent];
            }
        }
    }
    
    for (int i=0; i<node->nchildren; i++)
        estimateGenerate_helper(tree, node->children[i], 
                                depths, sroots, recon, events, free);    
}


float estimateGenerate(Tree *tree, SpeciesTree *stree, 
                       int *recon, int *events, 
                       float *mu, float *sigma, float alpha, float beta)
{
    float *depths = new float [tree->nnodes];
    int *sroots = new int [tree->nnodes];   // species roots
    
    depths[tree->root->name] = 0;
    sroots[tree->root->name] = recon[tree->root->name];
    
    // determine if top subtree is free
    bool free = (recon[tree->root->name] == stree->root->name && 
                 events[tree->root->name] == EVENT_DUP);
    
    estimateGenerate_helper(tree, tree->root, depths, sroots, 
                            recon, events, free);
    
    //printFloatArray(depths, tree->nnodes);
    //printIntArray(sroots, tree->nnodes);
    
    float *dists = new float [tree->nnodes];
    float *means = new float [tree->nnodes];
    float *sdevs = new float [tree->nnodes];
    
    int count = 0;
    
    for (int i=0; i<tree->nnodes; i++) {
        if (events[i] != EVENT_DUP && i != tree->root->name) {
            // we are at suvtree leaf
            
            // figure out species branches that we cross
            // get total mean and variance of this path            
            float u = 0.0;
            float s2 = 0.0;   
            int snode = recon[i];
            
            // branch is also free if we do not cross any more species
            // don't estimate baserates from extra branches
            if (snode != sroots[i] && depths[i] != -1) {
                while (snode != sroots[i] && snode != stree->root->name) {
                    u += mu[snode];
                    s2 += sigma[snode]*sigma[snode];
                    snode = stree->nodes[snode].parent->name;
                }
                if (fabs(s2) < .0000001) {
                    printf("gene tree\n");
                    printTree(tree);
                    
                    printf("species tree\n");
                    printTree(stree);
                    
                    printf("i         = %d\n", i);
                    printf("snode     = %d\n", snode);
                    printf("recon[i]  = %d\n", recon[i]);
                    printf("events[i] = %d\n", events[i]);
                    printf("events[p] = %d\n", events[tree->nodes[i].parent->name]);
                    printf("sroots[i] = %d\n", sroots[i]);
                    printf("depths[i] = %f\n", depths[i]);
                    assert(0);
                }
                
                // save dist and params
                dists[count] = depths[i];
                means[count] = u;
                sdevs[count] = sqrt(s2);
                count ++;
            }
        }
    }
    
    
    float generate = mleGenerate(count, dists, means, sdevs, alpha, beta);
    
    delete [] depths;
    delete [] sroots;
    delete [] dists;
    delete [] means;
    delete [] sdevs;    
    
    return generate;
}



// calculate the likelihood of rare events such as gene duplication
float rareEventsLikelihood(int nnodes, int *recon, int *events, int nsnodes,
                           float predupprob, float dupprob)
{
    float logl = 0.0;
    int sroot = nsnodes - 1;
    
    for (int i=0; i<nnodes; i++) {
        if (events[i] == EVENT_DUP) {
            if (recon[i] == sroot)
                logl += logf(predupprob);
            else
                logl += logf(dupprob);
        }
    }

    return logl;
}


// Reconcile a branch to the species tree
void reconBranch(int node, int *ptree, int *pstree, int *recon, int *events, 
                 float *mu, float *sigma,
                 ReconParams *reconparams)
{
    // set fractional branches
    if (recon[node] == recon[ptree[node]]) {
        // start reconciles to a subportion of species branch
        if (events[node] == EVENT_DUP)
            // only case k's are dependent
            reconparams->startfrac[node] = FRAC_DIFF;   // k[node] - k[node.parent]
        else
            reconparams->startfrac[node] = FRAC_PARENT; // 1.0 - k[node.parent]
        
        reconparams->startparams[node] = BranchParam(mu[recon[node]],
                                                     sigma[recon[node]]);

        // there is only one frac
        reconparams->endfrac[node] = FRAC_NONE;
        reconparams->endparams[node] = NULL_PARAM;
    } else {
        if (events[ptree[node]] == EVENT_DUP) {
            // start reconciles to last part of species branch
            reconparams->startfrac[node] = FRAC_PARENT; // 1.0 - k[node.parent]
            int snode = recon[ptree[node]];
            reconparams->startparams[node] = BranchParam(mu[snode], sigma[snode]);
        } else {
            reconparams->startfrac[node] = FRAC_NONE;
            reconparams->startparams[node] = NULL_PARAM;
        }

        if (events[node] == EVENT_DUP) {
            // end reconciles to first part of species branch
            reconparams->endfrac[node] = FRAC_NODE; // k[node]
            reconparams->endparams[node] = BranchParam(mu[recon[node]],
                                                       sigma[recon[node]]);
        } else {
            // end reconcile to at least one whole species branch
            reconparams->endfrac[node] = FRAC_NONE;
            reconparams->endparams[node] = NULL_PARAM;
        }
    }
    
    // set midparams
    if (recon[node] == recon[ptree[node]]) {
        // we begin and end on same branch
        // there are no midparams
        reconparams->midparams[node] = NULL_PARAM;
    } else {
        // we begin and end on different branches
        float totmean = 0.0;
        float totvar = 0.0;
        int snode;

        // determine most recent species branch which we fully recon to
        if (events[node] == EVENT_DUP)
            snode = pstree[recon[node]];
        else
            snode = recon[node];

        // walk up species spath until starting species branch
        // starting species branch is either fractional or NULL
        int parent_snode = recon[ptree[node]];
        while (snode != parent_snode) {
            totmean += mu[snode];
            totvar += sigma[snode] * sigma[snode];
            snode = pstree[snode];
        }

        reconparams->midparams[node] = BranchParam(totmean, sqrt(totvar));
    }
}


// Generate a random sample of duplication points
void setRandomMidpoints(int root, int *ptree,
                        int *subnodes, int nsubnodes, 
                        int *recon, int *events, 
                        ReconParams *reconparams)
{
    const float esp = .0001;
    
    reconparams->midpoints[ptree[root]] = 1.0;
    
    for (int i=0; i<nsubnodes; i++) {
        int node = subnodes[i];
        float lastpoint;
        
        if (events[node] == EVENT_DUP && ptree[node] != -1) {
            if (recon[node] == recon[ptree[node]])
                // if im the same species branch as my parent 
                // then he is my last midpoint
                lastpoint = reconparams->midpoints[ptree[node]];
            else
                // i'm the first on this branch so the last midpoint is zero
                lastpoint = 0.0;
            
            // pick a midpoint uniformly after the last one
            float remain = 1.0 - lastpoint;
            reconparams->midpoints[node] = lastpoint + esp * remain +
                                           (1.0-esp) * remain *
                                           (rand() / float(RAND_MAX));
        } else {
            // genes or speciations reconcile exactly to the end of the branch
            // gene tree roots also reconcile exactly to the end of the branch
            reconparams->midpoints[node] = 1.0;
        }
    }
}


// Calculate branch likelihood
float branchlk(float dist, int node, int *ptree, ReconParams *reconparams)
{
    float totmean = 0.0;
    float totvar = 0.0;
    float sigma;
    BranchParam bparam;
    
    
    float *k = reconparams->midpoints;
    
    int startfrac = reconparams->startfrac[node];
    bparam = reconparams->startparams[node];
    if (startfrac == FRAC_DIFF) {    
        totmean += (k[node] - k[ptree[node]]) * bparam.mu;
        totvar  += (k[node] - k[ptree[node]]) * bparam.sigma * bparam.sigma;
    } else if (startfrac == FRAC_PARENT) {
        totmean += (1.0 - k[ptree[node]]) * bparam.mu;
        totvar  += (1.0 - k[ptree[node]]) * bparam.sigma * bparam.sigma;
    }
    // startfrac == FRAC_NONE, do nothing
    
    bparam = reconparams->midparams[node];
    if (!bparam.isNull()) {
        totmean += bparam.mu;
        totvar  += bparam.sigma * bparam.sigma;
    }
    
    int endfrac = reconparams->endfrac[node];
    bparam = reconparams->endparams[node];    
    if (endfrac == FRAC_PARENT) {    
        totmean += (1.0 - k[ptree[node]]) * bparam.mu;
        totvar  += (1.0 - k[ptree[node]]) * bparam.sigma * bparam.sigma;
    } else if (endfrac == FRAC_NODE) {
        totmean += k[node] * bparam.mu;
        totvar  += k[node] * bparam.sigma * bparam.sigma;
    }
    // endfrac == FRAC_NONE, do nothing
    
    // unhandle partially-free branches and unfold
    if (reconparams->unfold == node)
        dist += reconparams->unfolddist;
    
    // skip a branch if it is partially free
    if (reconparams->freebranches[node]) {
        if (dist > totmean)
            dist = totmean;
    }
    
    return normallog(dist, totmean, sqrt(totvar));
}


// Calculate the likelihood of a subtree
float subtreelk(int nnodes, int *ptree, int **ftree, float *dists, int root,
                int nsnodes, int *pstree, 
                int *recon, int *events,
                float *mu, float *sigma, float generate,
                ReconParams *reconparams,
                int nsamples=100)
{
    float logl = 0.0;
    int sroot = nsnodes - 1;

    
    if (events[root] != EVENT_DUP) {
        // single branch, no integration needed
                
        if (recon[root] != sroot) {
            // no midpoints, no integration needed
            reconBranch(root, ptree, pstree, recon, events, 
                        mu, sigma, reconparams);
            logl = branchlk(dists[root] / generate, 
                            root, ptree, reconparams);
        }
    } else {
        // multiple branches, integrate
        
        // set reconparams by traversing subtree
        TreeWalker walk = TreeWalker(nnodes, ftree, root);
        int node;
        int *nodes =  new int [nnodes];
        int nodesi = 0;
        
        while ((node = walk.next()) != -1) {
            reconBranch(node, ptree, pstree, 
                        recon, events, mu, sigma, reconparams);
            nodes[nodesi++] = node;
            
            if (events[node] == EVENT_DUP) {
                walk.recurse(node);
            }
        }
        
        // choose number of samples based on number of nodes to integrate over
        nsamples = int(100*log(nodesi)) + 50;
        if (nsamples > 400) nsamples = 400;
        
        // perform integration by sampling
        float prob = 0.0;
        for (int i=0; i<nsamples; i++) {
            float sampleLogl = 0.0;
            
            // propose a setting of midpoints
            setRandomMidpoints(root, ptree, nodes, nodesi,
                               recon, events, reconparams);
            
            // loop through all branches in subtree
            for (int j=0; j<nodesi; j++) {
                int node = nodes[j];
                
                if (recon[node] != sroot) {
                    sampleLogl += branchlk(dists[node] / generate, 
                                           node, ptree, reconparams);
                }
            }
            
            prob += expf(sampleLogl) / nsamples;
        }
        
        logl = log(prob);
        
        // cleanup
        delete [] nodes;
    }
    
    return logl;
}


// Calculate the likelihood of a tree
float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate, float disterror,
             float predupprob, float dupprob, float errorlogl,
             float alpha, float beta)
{
    float logl = 0.0;
    int root = nnodes - 1;
    int sroot = nsnodes - 1;

    // make forward tree
    int **ftree;
    makeFtree(nnodes, ptree, &ftree);
    
    
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();
    
    // estimate generate
    float generate2 = estimateGenerate(&tree, &stree, 
                       recon, events,
                       mu, sigma, alpha, beta);
    
    //printf("generate: %f %f\n", generate, generate2);
    generate = generate2;
    
    /*
      find free branches

      A branch is an free branch if (1) its parent node reconciles to the species
      tree root, (2) it parent node is a duplication, (3) and the node it self
      reconciles not to the species tree root.

      A top branch unfolds, if (1) its parent is the root, (2) its parent
      reconciles to the species tree root, (3) its parent is a duplication, and
      (4) itself does not reconcile to the species tree root. There is atmost one
      unfolding branch per tree.

      If a branch is free, augment its length to min(dist, mean)
    */
    
    bool *freebranches = new bool [nnodes];
    int unfold = -1;
    float unfolddist = 0.0;
    
    for (int i=0; i<nnodes; i++) {
        if (recon[ptree[i]] == sroot &&
            events[ptree[i]] == EVENT_DUP &&
            recon[i] != sroot)
        {
            freebranches[i] = true;
        } else {
            freebranches[i] = false;
        }
    }
    
    // find unfolding branch
    if (ftree[root][0] != -1 &&
        recon[root] == sroot &&
        events[root] == EVENT_DUP)
    {
        if (recon[ftree[root][0]] != sroot) {
            unfold = ftree[root][0];
            unfolddist = dists[ftree[root][1]];
        } else {
            unfold = ftree[root][1];
            unfolddist = dists[ftree[root][1]];
        }
    }   
    
    ReconParams reconparams = ReconParams(nnodes, freebranches, 
                                          unfold, unfolddist);
    
    // loop through independent subtrees
    for (int i=0; i<nnodes; i++) {
        if (events[i] == EVENT_SPEC || i == root) {
            for (int j=0; j<2; j++) {
                logl += subtreelk(nnodes, ptree, ftree, dists, ftree[i][j],
                                  nsnodes, pstree, 
                                  recon, events,
                                  mu, sigma, generate,
                                  &reconparams);
            }
        }
    }
    
    
    // rare events
    logl += rareEventsLikelihood(nnodes, recon, events, nsnodes,
                                 predupprob, dupprob);
    
    
    // error cost
    logl += errorlogl * disterror;
    
        
    // generate probability
    if (alpha > 0 && beta > 0)
        logl += gammalog(generate, alpha, beta);
   
    
    // clean up
    freeFtree(nnodes, ftree);
    delete [] freebranches;
    
    return logl;
}
