//=============================================================================
//  SPIDIR - tree search



#include "common.h"
#include "Matrix.h"
#include "phylogeny.h"
#include "likelihood.h"
#include "parsimony.h"
#include "search.h"
#include "mldist.h"


namespace spidir {



//=============================================================================
// Nearest Neighbor Interchange Topology Proposal

/*

    Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      uncle    node1
               /  \
         child[0]  child[1]
    
    special case with rooted branch:
    
              node2
             /     \
        node2'      node1
       /     \     /     \
      uncle   * child[0] child[1]
    
      
    I do not need to renumber nodes is if child[change] and uncle are 
    both leaves or both all internal.    
    Renumbering is simply a swap of the two nodes otherwise
*/
void proposeNni(Tree *tree, Node *node1, Node *node2, Node *child)
{
    int uncle = 0;

    // ensure node1 is the child of node2
    if (node1->parent != node2) {
        Node *tmp = node1; node1 = node2; node2 = tmp;
    }
    assert(node1->parent == node2);
    
    
    // try to see if edge is one branch (not root edge)
    if (tree->isRooted() && node2 == tree->root) {
        // special case of specifying root edge
        if (node2->children[0] == node1)
            node2 = node2->children[1];
        else
            node2 = node2->children[0];
        
        // if edge is not an internal edge, give up
        if (node2->nchildren < 2)
            return;
    }
    
    if (node1->parent == tree->root &&
        node2->parent == tree->root)
    {
        uncle = 0;
        
        if (node2->children[0]->nchildren < 2 && 
            node2->children[1]->nchildren < 2) {
            // can't do NNI on this branch
            return;
        }
    } else {
        // find uncle
        uncle = 0;
        if (node2->children[uncle] == node1)
            uncle = 1;
    }
    
    int change = (node1->children[0] == child) ? 0 : 1;
    assert(node1->children[change] == child);
    
    // swap parent pointers
    //node1->children[change]->parent = node2;
    child->parent = node2;
    node2->children[uncle]->parent = node1;
    
    // swap child pointers
    Node *tmp = node2->children[uncle];
    node2->children[uncle] = child;
    node1->children[change] = tmp;    
    
    //node2->children[uncle] = node1->children[change];
    //node1->children[change] = tmp;
}


/*

    Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      nodeb    node1
               /  \
         nodea     * 

*/
void proposeNni(Tree *tree, Node *nodea, Node *nodeb)
{
    Node *node1 = nodea->parent;
    Node *node2 = nodeb->parent;
    
    // assert that node1 and node2 are incident to the same branch
    assert(node1->parent == node2 ||
           node2->parent == node1);
    
    // find child indexes
    int a = (node1->children[0] == nodea) ? 0 : 1;
    assert(node1->children[a] == nodea);

    int b = (node2->children[0] == nodeb) ? 0 : 1;
    assert(node2->children[b] == nodeb);
    
    // swap parent pointers
    nodea->parent = node2;
    nodeb->parent = node1;
    
    // swap child pointers
    node2->children[b] = nodea;
    node1->children[a] = nodeb;
}


void proposeRandomNni(Tree *tree, Node **a, Node **b)
{
    // find edges for NNI
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->isLeaf() || 
             tree->nodes[choice]->parent == NULL);
    
    Node *node1 = tree->nodes[choice];
    Node *node2 = tree->nodes[choice]->parent;
    *a = node1->children[irand(2)];
    *b = (node2->children[0] == node1) ? node2->children[1] :
                                         node2->children[0];
}


NniProposer::NniProposer(SpeciesTree *stree, int *gene2species,
                         int niter) :
    nodea(NULL),
    nodeb(NULL),
    nodec(NULL),
    noded(NULL),
    oldroot1(NULL),
    oldroot2(NULL),
    stree(stree),
    gene2species(gene2species),
    niter(niter),
    iter(0),
    correctTree(NULL),
    correctSeen(false)
{}
    
    
void NniProposer::propose(Tree *tree)
{
    const float rerootProb = 1.0;
    const float doubleNniProb = 0.1;
    
    // increase iteration
    iter++;
    
    nodea = nodeb = nodec = noded = NULL;

    // propose new tree
    proposeRandomNni(tree, &nodea, &nodeb);
    proposeNni(tree, nodea, nodeb);

    if (frand() < doubleNniProb) {
        printLog(LOG_MEDIUM, "search: double NNI\n");
        proposeRandomNni(tree, &node3, &node4, &nodec, &noded);
        proposeNni(tree, nodec, noded);
    }
    
    // reroot tree if stree is given
    if (frand() < rerootProb) {
        oldroot1 = tree->root->children[0];
        oldroot2 = tree->root->children[1];
        
        if (stree != NULL) {
            reconRoot(tree, stree, gene2species);
        }
    } else {
        oldroot1 = NULL;
        oldroot2 = NULL;
    }
    
    if (correctTree) {
        if (tree->sameTopology(correctTree))
            correctSeen = true;
    }

    assert(tree->assertTree());
}

void NniProposer::revert(Tree *tree)
{
    // reject, undo topology change

    if (oldroot1) {
        if (oldroot1->parent == oldroot2)
            tree->reroot(oldroot1);
        else
            tree->reroot(oldroot2);
    }
    
    //printf("NNI %d %d %d %d\n", node1->name, node1->parent->name, 
    //       node2->name, node2->nchildren);
    if (nodec)
        proposeNni(tree, nodec, noded);
    if (nodea)
        proposeNni(tree, nodea, nodeb);
}


bool NniProposer::more()
{
    return iter < niter;
}

//=============================================================================
// Subtree pruning and regrafting (SPR)


/*
    a = subtree
    e = newpos
    
    BEFORE
            ....
        f         d
       /           \
      c             e
     / \           ...
    a   b
   ... ...

    AFTER

        f         d
       /           \
      b             c
     ...           / \
                  a   e
                 ... ...

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary

*/
void proposeSpr(Tree *tree, Node *subtree, Node *newpos)
{
    Node *a = subtree;
    Node *e = newpos;
    
    Node *c = a->parent;
    Node *f = c->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    const int ci = (f->children[0] == c) ? 0 : 1;
    Node *d = e->parent;
    const int ei = (d->children[0] == e) ? 0 : 1;

    d->children[ei] = c;
    c->children[bi] = e;
    f->children[ci] = b;
    b->parent = f;
    c->parent = d;
    e->parent = c;
}

/*
    What if e == f  (also equivalent to NNI) this is OK

    BEFORE
    
          d
         / \
        e  ...
       / \
      c  ...         
     / \           
    a   b
   ... ...

    AFTER
          d
         / \
        c
       / \
      a   e
     ... / \
        b  ...
       ...
       
  What if d == f  (also equivalent to NNI) this is OK
  
    BEFORE
          
        f
       / \
      c   e
     / \  ...
    a   b
   ... ...

    AFTER
          
        f
       / \
      b   c  
     ... / \ 
        a   e
       ... ...  
*/



/*
    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary
*/
void proposeRandomSpr(Tree *tree, Node **subtree, Node **newpos)
{
    assert(tree->nnodes >= 5);

    // find subtree (a) to cut off (any node that is not root or child of root)
    int choice;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->parent == NULL ||
             tree->nodes[choice]->parent->parent == NULL);
    *subtree = tree->nodes[choice];
    Node *a = *subtree;
    
    // find sibling (b) of a
    Node *c = a->parent;
    const int bi = (c->children[0] == a) ? 1 : 0;
    Node *b = c->children[bi];
    
    // choose newpos (e)
    Node *e = NULL;
    do {
        choice = irand(tree->nnodes);
        Node *e = tree->nodes[choice];
        
        // test if e is a valid choice
        if (e->parent == NULL || e == a || e == c || e == b)
            continue;
        
        // also test if e is a descendent of a
        bool under_a = false;
        for (Node *ptr = e->parent; ptr != NULL; ptr = ptr->parent) {
            if (ptr == a) {
                under_a = true;
                break;
            }
        }            
        
        if (under_a)
            continue;
    } while (false);
    *newpos = e;
}


//=============================================================================
// Fitting branch lengths

ParsimonyFitter::ParsimonyFitter(int nseqs, int seqlen, char **seqs) :
    nseqs(nseqs),
    seqlen(seqlen),
    seqs(seqs)
{}


float ParsimonyFitter::findLengths(Tree *tree)
{
    parsimony(tree, nseqs, seqs);
    return 0.0;
}



HkyFitter::HkyFitter(int nseqs, int seqlen, char **seqs, 
                     float *bgfreq, float tsvratio, int maxiter,
                     bool useLogl) :
    nseqs(nseqs),
    seqlen(seqlen),
    seqs(seqs),
    bgfreq(bgfreq),
    tsvratio(tsvratio),
    maxiter(maxiter),
    useLogl(useLogl)
{
}

float HkyFitter::findLengths(Tree *tree)
{ 
    float logl = findMLBranchLengthsHky(tree, nseqs, seqs, bgfreq, 
                                        tsvratio, maxiter);
    if (useLogl)
        return logl;
    else
        return 0.0;
}


//=============================================================================
// Likelihood function

SpidirBranchLikelihoodFunc::SpidirBranchLikelihoodFunc(
    int nnodes, SpeciesTree *stree, 
    SpidirParams *params, 
    int *gene2species,
    float predupprob, float dupprob, bool estGenerate, bool onlyduploss) :
    
    nnodes(nnodes),
    stree(stree),
    params(params),
    gene2species(gene2species),
    recon(nnodes),
    events(nnodes),
    predupprob(predupprob),
    dupprob(dupprob),
    estGenerate(estGenerate),
    onlyduploss(onlyduploss)
{}


float SpidirBranchLikelihoodFunc::likelihood(Tree *tree) {
    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    float generate;
    
    if (estGenerate)
        generate = -1;
    else
        generate = -99;
    
    eventslk = rareEventsLikelihood(tree, stree, recon, 
                                    events, predupprob, dupprob);

    
    return treelk(tree, stree,
                  recon, events, params,
                  generate, predupprob, dupprob, onlyduploss);
}



float SpidirBranchLikelihoodFunc::likelihood2(Tree *tree) {
    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    
    return rareEventsLikelihood(tree, stree, recon, 
                                events, predupprob, dupprob);
}


//=============================================================================


// propose initial tree by Neighbor Joining
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs)
{
    int nnodes = nseqs * 2 - 1;

    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    Matrix<float> distmat(nseqs, nseqs);

    calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
    neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);

    Tree *tree = new Tree(nnodes);
    ptree2tree(nnodes, ptree, tree);
    tree->setLeafNames(genes);

    return tree;
}


// propose initial tree and root by species
Tree *getInitialTree(string *genes, int nseqs, int seqlen, char **seqs,
                     SpeciesTree *stree, int *gene2species)
{
    Tree *tree = getInitialTree(genes, nseqs, seqlen, seqs);
    reconRoot(tree, stree, gene2species);

    return tree;
}



//=============================================================================
// MCMC search

// NOTE: not used any more
// DEPRECATED,
Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                string *genes, int nseqs, int seqlen, char **seqs,
                int niter, float predupprob, float dupprob)
{
        
    NniProposer nniProposer(stree, gene2species, niter);
    ParsimonyFitter parsFitter(nseqs, seqlen, seqs);
    
    int nnodes = nseqs * 2 - 1;
    SpidirBranchLikelihoodFunc lkfunc(nnodes, stree, params, gene2species,
                                      predupprob, dupprob, -99);
    
    return searchMCMC(initTree, 
                      genes, nseqs, seqlen, seqs,
                      &lkfunc,
                      &nniProposer,
                      &parsFitter);
}



void printSearchStatus(Tree *tree, SpeciesTree *stree, int *gene2species,
                       int *recon=NULL, int *events=NULL)
{
    if (stree) {
        bool cleanupRecon = false;
        if (!recon) {
            recon = new int [tree->nnodes];
            cleanupRecon = true;
        }

        bool cleanupEvents = false;    
        if (!events) {
            events = new int [tree->nnodes];
            cleanupEvents = true;
        }    

    
    
        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        int losses = countLoss(tree, stree, recon);        
        
        // count dups
        int dups = 0;
        for (int i=0; i<tree->nnodes; i++)
            if (events[i] == EVENT_DUP)
                dups++;
        
        printLog(LOG_LOW, "search: dups = %d\n", dups);
        printLog(LOG_LOW, "search: loss = %d\n", losses);
        
        if (cleanupRecon)
            delete [] recon;
        if (cleanupEvents)
            delete [] events;        
    }                
    
    if (isLogLevel(LOG_LOW))
        displayTree(tree, getLogFile());
    
    if (isLogLevel(LOG_MEDIUM)) {
        printLog(LOG_MEDIUM, "tree: ");
        tree->writeNewick(getLogFile(), NULL, 0, true);
        printLog(LOG_MEDIUM, "\n\n");
    }
}


Tree *searchMCMC(Tree *initTree, 
                 string *genes, int nseqs, int seqlen, char **seqs,
                 BranchLikelihoodFunc *lkfunc,
                 TopologyProposer *proposer,
                 BranchLengthFitter *fitter)
{
    Tree *toptree = NULL;
    float toplogl = -1e10, logl=-1e10, nextlogl;
    Tree *tree = initTree;
    //int nnodes = nseqs * 2 - 1;
        
    
    // determine initial tree
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs);
    
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // init likelihood score
    parsimony(tree, nseqs, seqs); // get initial branch lengths TODO: replace?
    logl = fitter->findLengths(tree);
    logl += lkfunc->likelihood(tree);
    toplogl = logl;
    toptree = tree->copy();
    
    float speed = 0;
    
    typedef ExtendArray<int> TopologyKey;
    typedef pair<Tree*,float> TreeLogl;
    TopologyKey key(tree->nnodes);

    HashTable<TopologyKey, TreeLogl, HashTopology> hashtrees(2000, 
                                                   TreeLogl(NULL, 0));

    
    // MCMC loop
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposer->propose(tree);
        
        // hash topology
        //tree->hashkey(key);
        
        // TURN OFF hashing
        // have we seen this topology before?
        //TreeLogl &tmp = hashtrees[key];
        //if (0 && tmp.first != NULL) {
        //    // retrieve previously seen tree and logl
        //    nextlogl = tmp.second;
        //} else {    
            // calculate likelihood
            float seqlk = fitter->findLengths(tree);
            float branchlk = lkfunc->likelihood(tree);
            nextlogl = seqlk + branchlk;
            
            // hash result
            //tmp = TreeLogl(tree->copy(), nextlogl);
        //}
        
        
        // acceptance rule
        if (nextlogl > logl ||
            nextlogl - logl + speed > log(frand()))
        {
            printLog(LOG_MEDIUM, "search: accept %f  %f\n", nextlogl, logl);
            // accept
            logl = nextlogl;
            speed /= 2.0;

            // keep track of toptree            
            if (logl > toplogl) {
                float logl2 = lkfunc->likelihood2(tree);
            
                printLog(LOG_LOW, "search: lnl   = %f\n", logl);
                printLog(LOG_LOW, "search: seqlk = %f\n", seqlk);
                printLog(LOG_LOW, "search: lenlk = %f\n", branchlk);
                printLog(LOG_LOW, "search: lnl2  = %f\n", logl2);
                printSearchStatus(tree, lkfunc->getSpeciesTree(), 
                                  lkfunc->getGene2species(), recon, events);
                
                delete toptree;
                speed = 0.0;
                toptree = tree->copy();
                toplogl = logl;
            } else {
                if (isLogLevel(LOG_MEDIUM))
                    displayTree(toptree, getLogFile());
            }
        } else {
            printLog(LOG_MEDIUM, "search: reject %f < %f\n", nextlogl, toplogl);
            if (isLogLevel(LOG_MEDIUM))
                displayTree(tree, getLogFile());
            
            //speed = (speed + 1) * 1.3;
            
            // reject, undo topology change
            proposer->revert(tree);
        }
    }
    
    
    return toptree;
}


Tree *searchNni(Tree *initTree, 
                string *genes, int nseqs, int seqlen, char **seqs,
                BranchLikelihoodFunc *lkfunc,
                TopologyProposer *proposer,
                BranchLengthFitter *fitter)
{
    Tree *toptree = NULL;
    float toplogl = -1e10, nextlogl;
    Tree *tree = initTree;
    //int nnodes = nseqs * 2 - 1;
        
    
    // determine initial tree
    if (initTree == NULL)
        tree = getInitialTree(genes, nseqs, seqlen, seqs);
    
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    
    // init likelihood score
    parsimony(tree, nseqs, seqs); // get initial branch lengths TODO: replace?
    toplogl = fitter->findLengths(tree);
    toplogl += lkfunc->likelihood(tree);
    toptree = tree->copy();
    
    
    // MCMC loop
    for (int i=0; proposer->more(); i++) {
        printLog(LOG_LOW, "search: iter %d\n", i);
    
        // propose new tree 
        proposer->propose(tree);
        
        // calculate likelihood
        float seqlk = fitter->findLengths(tree);
        float branchlk = lkfunc->likelihood(tree);
        nextlogl = seqlk + branchlk;
        
        
        // acceptance rule
        if (nextlogl > toplogl)
        {
            float logl2 = lkfunc->likelihood2(tree);
        
            printLog(LOG_MEDIUM, "search: accept %f\n", nextlogl);
            printLog(LOG_LOW, "search: lnl   = %f\n", nextlogl);
            printLog(LOG_LOW, "search: seqlk = %f\n", seqlk);
            printLog(LOG_LOW, "search: lenlk = %f\n", branchlk);          
            printLog(LOG_LOW, "search: lnl2  = %f\n", logl2);
                        
            // accept
            toplogl = nextlogl;
            
            printSearchStatus(tree, lkfunc->getSpeciesTree(), 
                  lkfunc->getGene2species(), recon, events);
            delete toptree;
            toptree = tree->copy();
        } else {
            // reject, undo topology change
            proposer->revert(tree);
        }
    }
    
    
    return toptree;
}

} // namespace spidir
