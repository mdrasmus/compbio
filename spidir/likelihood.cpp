//=============================================================================
//  SPIDIR - Likelihood calculation


// c++ headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

// spidir headers
#include "common.h"
#include "likelihood.h"
#include "phylogeny.h"
#include "ExtendArray.h"


namespace spidir {



// fractional branches
enum {
    FRAC_NONE,
    FRAC_DIFF,
    FRAC_PARENT,
    FRAC_NODE
};


// Branch distribution parameters for one branch
class BranchParams
{
public:
    BranchParams(float _mu=-1.0, float _sigma=-1.0) :
        mu(_mu),
        sigma(_sigma)
    {}
    
    bool isNull()
    {
        return mu == -1.0;
    }
    
    float mu;
    float sigma;
};

BranchParams NULL_PARAM;


// Reconciliation parameters
class ReconParams
{
public:
    ReconParams(int nnodes) :
        nnodes(nnodes),
        unfold(-1),
        unfolddist(0)
    {
        startparams = new BranchParams [nnodes];
        midparams = new BranchParams [nnodes];
        endparams = new BranchParams [nnodes];
        
        startfrac = new int [nnodes];
        endfrac = new int [nnodes];
        midpoints = new float [nnodes];
        
        freebranches = new bool [nnodes];
    }
    
    ~ReconParams()
    {
        delete [] startparams;
        delete [] midparams;
        delete [] endparams;
        
        delete [] startfrac;
        delete [] endfrac;
        delete [] midpoints;
        
        delete [] freebranches;
    }
    
    
    int nnodes;
    BranchParams *startparams;
    BranchParams * midparams;
    BranchParams *endparams;
    int *startfrac;
    int *endfrac;
    float *midpoints;
    bool *freebranches;
    int unfold;
    float unfolddist;
};


//=============================================================================
// gene rate estimation


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

int floatcmp(const void *a, const void *b)
{
    float fa = *((float*)a);
    float fb = *((float*)b);
    
    if (fa < fb)
        return -1;
    else if (fa > fb)
        return 1;
    else
        return 0;
}


float mleGenerate(int count, float *dists, float *means, float *sdevs,
                  float alpha, float beta)
{
    float a, b, c;
    float threshold = 0;
    
    a = (1.0 - alpha) / beta;

    // b = sum(means[i] * lens[i] / sdevs[i]**2
    //         for i in range(len(lens))) / beta    
    // c = - sum(lens[i] ** 2 / sdevs[i] ** 2
    //          for i in range(len(lens))) / beta    
    

    ExtendArray<float> dists2(0, count);
    dists2.extend(dists, count);
    qsort((void*) dists2.get(), dists2.size(), sizeof(float),
                  floatcmp);
    //printFloatArray(dists2, dists2.size());
    //printFloatArray(dists, dists2.size());
    int limit = int(count * .5) + 1;
    if (limit < 4) limit = 4;
    threshold = dists2[limit];
    //printf("threshold %f\n", threshold);
    
    
    b = 0.0;
    c = 0.0;    
    for (int i=0; i<count; i++) {
        if (dists[i] > threshold && sdevs[i] > 0.0001) {
            b += means[i] * dists[i] / (sdevs[i] * sdevs[i]);
            c += dists[i] * dists[i] / (sdevs[i] * sdevs[i]);
        }
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
                       int *recon, int *events, SpidirParams *params)
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
    
    // make quick access to params
    float *mu = params->mu;
    float *sigma = params->sigma;
    
    int count = 0;
    
    for (int i=0; i<tree->nnodes; i++) {
        if (events[i] != EVENT_DUP && i != tree->root->name) {
            // we are at subtree leaf
            
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
                    snode = stree->nodes[snode]->parent->name;
                }
                assert(fabs(s2) > .0000001);
                                
                // save dist and params
                dists[count] = depths[i];
                means[count] = u;
                sdevs[count] = sqrt(s2);
                count ++;
            }
        }
    }
    
    
    float generate = mleGenerate(count, dists, means, sdevs, 
                                 params->alpha, params->beta);
    
    delete [] depths;
    delete [] sroots;
    delete [] dists;
    delete [] means;
    delete [] sdevs;    
    
    return generate;
}


//=============================================================================
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


//=============================================================================
// likelihood functions

// Reconcile a branch to the species tree
void reconBranch(int node, int *ptree, int *pstree, int *recon, int *events, 
                 SpidirParams *params,
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
        
        reconparams->startparams[node] = BranchParams(params->mu[recon[node]],
                                                      params->sigma[recon[node]]);

        // there is only one frac
        reconparams->endfrac[node] = FRAC_NONE;
        reconparams->endparams[node] = NULL_PARAM;
    } else {
        if (events[ptree[node]] == EVENT_DUP) {
            // start reconciles to last part of species branch
            reconparams->startfrac[node] = FRAC_PARENT; // 1.0 - k[node.parent]
            int snode = recon[ptree[node]];
            reconparams->startparams[node] = BranchParams(params->mu[snode],
                                                          params->sigma[snode]);
        } else {
            reconparams->startfrac[node] = FRAC_NONE;
            reconparams->startparams[node] = NULL_PARAM;
        }

        if (events[node] == EVENT_DUP) {
            // end reconciles to first part of species branch
            reconparams->endfrac[node] = FRAC_NODE; // k[node]
            reconparams->endparams[node] = BranchParams(params->mu[recon[node]],
                                                        params->sigma[recon[node]]);
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
        
        int parent_snode;
        if (ptree[node] != -1)
            parent_snode = recon[ptree[node]];
        else
            parent_snode = -1;
        while (snode != parent_snode && pstree[snode] != -1) {
            totmean += params->mu[snode];
            totvar += params->sigma[snode] * params->sigma[snode];
            snode = pstree[snode];
        }
        
        reconparams->midparams[node] = BranchParams(totmean, sqrt(totvar));
    }
}


// Generate a random sample of duplication points
void setRandomMidpoints(int root, int *ptree,
                        int *subnodes, int nsubnodes, 
                        int *recon, int *events, 
                        ReconParams *reconparams)
{
    const float esp = .0001;
    
    // this should not be here
    //reconparams->midpoints[ptree[root]] = 1.0;
    
    for (int i=0; i<nsubnodes; i++) {
        int node = subnodes[i];
        
        if (events[node] == EVENT_DUP && ptree[node] != -1) {
            float lastpoint;
            
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
                                           (1.0-2*esp) * remain * frand();
        } else {
            // genes or speciations reconcile exactly to the end of the branch
            // gene tree roots also reconcile exactly to the end of the branch
            reconparams->midpoints[node] = 1.0;
        }
    }
}


BranchParams getBranchParams(int node, int *ptree, ReconParams *reconparams)
{
    float totmean = 0.0;
    float totvar = 0.0;
    BranchParams bparam;
    
    
    float *k = reconparams->midpoints;
    
    int startfrac = reconparams->startfrac[node];
    bparam = reconparams->startparams[node];
    if (startfrac == FRAC_DIFF) {    
        totmean += (k[node] - k[ptree[node]]) * bparam.mu;
        totvar  += (k[node] - k[ptree[node]]) * bparam.sigma * bparam.sigma;
    } else if (startfrac == FRAC_PARENT) {
        float kp = 0;
        
        if (!reconparams->freebranches[node])
            kp = k[ptree[node]];
    
        totmean += (1.0 - kp) * bparam.mu;
        totvar  += (1.0 - kp) * bparam.sigma * bparam.sigma;
    }
    // startfrac == FRAC_NONE, do nothing
    
    bparam = reconparams->midparams[node];
    if (!bparam.isNull()) {
        if (bparam.mu < 0) {
            fprintf(stderr, "bparam %f\n", bparam.mu);
            fprintf(stderr, "node %d\n", node);
            assert(0);
        }
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
    
    return BranchParams(totmean, sqrt(totvar));
}


// Calculate branch likelihood
float branchlk(float dist, int node, int *ptree, ReconParams *reconparams)
{
    BranchParams bparam = getBranchParams(node, ptree, reconparams);
    float totmean = bparam.mu;
    float totsigma = bparam.sigma;
    
    // handle partially-free branches and unfold
    if (reconparams->unfold == node)
        dist += reconparams->unfolddist;
    
    // augment a branch if it is partially free
    if (reconparams->freebranches[node]) {
        if (dist > totmean)
            dist = totmean;
    }
    
    float logl = normallog(dist, totmean, totsigma);
    assert(!isnan(logl));
    return logl;
}


void getSubtree(int **ftree, int node, int *events, ExtendArray<int> *subnodes)
{
    subnodes->append(node);

    // recurse
    if (events[node] == EVENT_DUP) {
        if (ftree[node][0] != -1)
            getSubtree(ftree, ftree[node][0], events, subnodes);
        if (ftree[node][1] != -1)
            getSubtree(ftree, ftree[node][1], events, subnodes);
    }
}


void getSubtree(Node *node, int *events, ExtendArray<Node*> *subnodes)
{
    subnodes->append(node);

    // recurse
    if (events[node->name] == EVENT_DUP) {
        for (int i=0; i<node->nchildren; i++) 
            getSubtree(node->children[i], events, subnodes);
    }
}


// Calculate the likelihood of a subtree
float subtreelk(int nnodes, int *ptree, int **ftree, float *dists, int root,
                int nsnodes, int *pstree, 
                int *recon, int *events, SpidirParams *params,
                float generate,
                ReconParams *reconparams,
                int nsamples=100)
{
    float logl = 0.0;
    int sroot = nsnodes - 1;

    
    if (events[root] != EVENT_DUP) {
        // single branch, no integration needed
                
        if (recon[root] != sroot) {
            // no midpoints, no integration needed
            reconBranch(root, ptree, pstree, recon, events, params,
                        reconparams);
            logl = branchlk(dists[root] / generate, 
                            root, ptree, reconparams);
        }
    } else {
        // multiple branches, integrate
        
        // set reconparams by traversing subtree
        ExtendArray<int> subnodes(0, nnodes);
        getSubtree(ftree, root, events, &subnodes);
        
        
        for (int i=0; i<subnodes.size(); i++) {
            reconBranch(subnodes[i], ptree, pstree, 
                        recon, events, params, reconparams);
        }
        
        // choose number of samples based on number of nodes to integrate over
        nsamples = int(500*logf(subnodes.size())) + 500;
        if (nsamples > 2000) nsamples = 2000;
        //nsamples = int(100*log(nodesi)) + 50;
        //if (nsamples > 600) nsamples = 600;

        
        // perform integration by sampling
        double prob = 0.0;
        for (int i=0; i<nsamples; i++) {
            double sampleLogl = 0.0;
            
            // propose a setting of midpoints
            reconparams->midpoints[root] = 1.0; // TODO: need to understand why this is here
            setRandomMidpoints(root, ptree, subnodes, subnodes.size(),
                               recon, events, reconparams);
            
            // loop through all branches in subtree
            for (int j=0; j<subnodes.size(); j++) {
                int node = subnodes[j];
                
                if (recon[node] != sroot) {
                    sampleLogl += branchlk(dists[node] / generate, 
                                           node, ptree, reconparams);
                }
            }
            
            prob += exp(sampleLogl);
        }
        
        logl = log(prob  / nsamples);
    }

    assert(!isnan(logl));
    return logl;
}


void determineFreeBranches(Tree *tree, SpeciesTree *stree, 
                           int *recon, int *events, float generate,
                           int *unfold, float *unfolddist, bool *freebranches)
{
    /*
      find free branches

      A branch is a (partially) free branch if (1) its parent node reconciles to
      the species tree root, (2) it parent node is a duplication, (3) and the
      node it self reconciles not to the species tree root.

      A top branch unfolds, if (1) its parent is the root, (2) its parent
      reconciles to the species tree root, (3) its parent is a duplication, and
      (4) itself does not reconcile to the species tree root. There is atmost one
      unfolding branch per tree.

      If a branch is free, augment its length to min(dist, mean)
    */
    
    int sroot = stree->root->name;
    
    *unfold = -1;
    *unfolddist = 0.0;
    
    for (int i=0; i<tree->nnodes; i++) {
        if (tree->nodes[i]->parent &&
            recon[tree->nodes[i]->parent->name] == sroot &&
            events[tree->nodes[i]->parent->name] == EVENT_DUP &&
            recon[i] != sroot)
        {
            freebranches[i] = true;
        } else {
            freebranches[i] = false;
        }
    }
    
    // find unfolding branch
    if (tree->root->nchildren >= 2 &&
        recon[tree->root->name] == sroot &&
        events[tree->root->name] == EVENT_DUP)
    {
        if (recon[tree->root->children[0]->name] != sroot) {
            *unfold = tree->root->children[0]->name;
            *unfolddist = tree->root->children[1]->dist / generate;
        } else {
            *unfold = tree->root->children[1]->name;
            *unfolddist = tree->root->children[0]->dist / generate;
        }
    }
}


float treelk(Tree *tree,
             SpeciesTree *stree,
             int *recon, int *events, SpidirParams *params,
             float generate, float disterror,
             float predupprob, float dupprob, float errorlogl)
{
    float logl = 0.0; // log likelihood
    
    // make parent trees
    int *ptree = new int [tree->nnodes];
    tree2ptree(tree, ptree);
    
    int *pstree = new int [stree->nnodes];
    tree2ptree(stree, pstree);
    
    // make forward tree
    int **ftree;
    makeFtree(tree->nnodes, ptree, &ftree);
    
    float *dists = new float [tree->nnodes];
    tree->getDists(dists);
    
    // estimate generate
    if (generate <= 0)
        generate = estimateGenerate(tree, stree, 
                       recon, events, params);
    if (generate == 0.0)
        generate = estimateGenerate(tree, stree, 
                       recon, events, params);
    
    
    // determine reconciliation parameters
    ReconParams reconparams = ReconParams(tree->nnodes);
    determineFreeBranches(tree, stree, recon, events, generate,
                          &reconparams.unfold, 
                          &reconparams.unfolddist, 
                          reconparams.freebranches);
    //printf("UNFOLD: %d %f\n", unfold, unfolddist);    

    
        
    // loop through independent subtrees
    for (int i=0; i<tree->nnodes; i++) {
        if (events[i] == EVENT_SPEC || i == tree->root->name) {
            if (events[i] == EVENT_SPEC) {
                for (int j=0; j<2; j++) {
                    int node = tree->nodes[i]->children[j]->name;
                    float slogl = subtreelk(tree->nnodes, ptree, ftree, dists, 
                                            node,
                                            stree->nnodes, pstree, 
                                            recon, events, params,
                                            generate,
                                            &reconparams);
                    logl += slogl;
                }
            } else {
                float slogl = subtreelk(tree->nnodes, ptree, ftree, dists, 
                                        i,
                                        stree->nnodes, pstree, 
                                        recon, events, params,
                                        generate,
                                        &reconparams);
                logl += slogl;
            }
        }
    }
    
    // rare events
    logl += rareEventsLikelihood(tree->nnodes, recon, events, stree->nnodes,
                                 predupprob, dupprob);
    
    
    // error cost
    logl += errorlogl * disterror;
    
        
    // generate probability
    if (params->alpha > 0 && params->beta > 0)
        logl += gammalog(generate, params->alpha, params->beta);
    
    // clean up
    delete [] ptree;
    delete [] pstree;
    delete [] dists;
    freeFtree(tree->nnodes, ftree);
    
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
    // create tree objects
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();    

    
    SpidirParams params = SpidirParams(nsnodes, NULL, mu, sigma, alpha, beta);
    
    return treelk(&tree, &stree,
                  recon, events, &params, 
                  generate, disterror,
                  predupprob, dupprob, errorlogl);
}


//=============================================================================
// branch length generation

float genBranch(float generate, float mean, float sdev)
{
    float blen = 0;
    
    while (blen <= 0)
        blen = generate * normalvariate(mean, sdev);
    
    return blen;
}



// NOTE: currently generates zero branch length for freebranches...
void genSubtree(Tree *tree, Node *root,
                SpeciesTree *stree,
                int *ptree, int *pstree,
                int *recon, int *events, SpidirParams *params,
                float generate,
                ReconParams *reconparams)
{
    int sroot = stree->root->name;

    
    if (events[root->name] != EVENT_DUP) {
        // single branch
                
        if (recon[root->name] != sroot) {
            // no midpoints
            reconBranch(root->name, ptree, pstree, recon, events, params,
                        reconparams);
            BranchParams bparam = getBranchParams(root->name, ptree, reconparams);            
            root->dist = genBranch(generate, bparam.mu, bparam.sigma);
        }
    } else {
        // multiple branches
                
        // set reconparams by traversing subtree
        ExtendArray<Node*> subnodes(0, tree->nnodes);
        getSubtree(tree->root, events, &subnodes);
        
        ExtendArray<int> subnames(subnodes.size());
        for (int i=0; i<subnodes.size(); i++)
            subnames[i] = subnodes[i]->name;
        
        
        for (int i=0; i<subnodes.size(); i++) {
            reconBranch(subnodes[i]->name, ptree, pstree, 
                        recon, events, params, reconparams);
        }
        
        
        // propose a setting of midpoints
        reconparams->midpoints[root->name] = 1.0; // TODO: need to understand why this is here
        setRandomMidpoints(root->name, ptree, subnames, subnodes.size(),
                           recon, events, reconparams);

        // loop through all branches in subtree
        for (int j=0; j<subnodes.size(); j++) {
            Node *node = subnodes[j];

            if (recon[node->name] != sroot) {
                BranchParams bparam = getBranchParams(node->name, ptree, reconparams);
                node->dist = genBranch(generate, bparam.mu, bparam.sigma);
            }
        }
    }
}


void generateBranchLengths(Tree *tree,
                           SpeciesTree *stree,
                           int *recon, int *events,
                           SpidirParams *params)
{
    // generate a generate
    float generate = gammavariate(params->alpha, params->beta);
    
    
    // determine reconciliation parameters
    ReconParams reconparams = ReconParams(tree->nnodes);
    determineFreeBranches(tree, stree, recon, events, generate,
                          &reconparams.unfold, 
                          &reconparams.unfolddist, 
                          reconparams.freebranches);
    
    
    // make array formats
    ExtendArray<int> ptree(tree->nnodes);
    ExtendArray<int> pstree(stree->nnodes);
    tree2ptree(tree, ptree);
    tree2ptree(stree, pstree);
        
    
    // loop through independent subtrees
    for (int i=0; i<tree->nnodes; i++) {
        if (events[i] == EVENT_SPEC || i == tree->root->name) {
            for (int j=0; j<2; j++) {
                Node *node = tree->nodes[i]->children[j];
                genSubtree(tree, node,
                           stree, ptree, pstree,
                           recon, events, params,
                           generate,
                           &reconparams);
            }
        }
    }
}


void generateBranchLengths(int nnodes, int *ptree, 
                           int nsnodes, int *pstree,
                           int *recon, int *events,
                           float *mu, float *sigma,
                           float alpha, float beta,
                           float *dists)
{
    // create gene tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);
    
    // create species tree object
    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();  
    
    // build parameters object
    SpidirParams params(nsnodes, NULL, mu, sigma, alpha, beta);
    
    generateBranchLengths(&tree, &stree, recon, events, &params);
    
    // record distances into array
    tree.getDists(dists);
}
                         


} // namespace spidir
