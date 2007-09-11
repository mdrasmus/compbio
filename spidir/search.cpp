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
void proposeNni(Tree *tree, Node *node1, Node *node2, int change)
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
    
    // swap parent pointers
    node1->children[change]->parent = node2;
    node2->children[uncle]->parent = node1;
    
    // swap child pointers
    Node *tmp = node2->children[uncle];
    node2->children[uncle] = node1->children[change];
    node1->children[change] = tmp;
}



void proposeRandomNni(Tree *tree, Node **node1, Node **node2, int *change)
{
    // find edges for NNI
    int choice = tree->root->name;
    do {
        choice = irand(tree->nnodes);
    } while (tree->nodes[choice]->isLeaf() || 
             tree->nodes[choice]->parent == NULL);
    
    *node1 = tree->nodes[choice];
    *node2 = tree->nodes[choice]->parent;
    *change = irand(2);
}


NniProposer::NniProposer(SpeciesTree *stree, int *gene2species,
                         int niter) :
    node1(NULL),
    node2(NULL),
    node3(NULL),
    node4(NULL),
    stree(stree),
    gene2species(gene2species),
    niter(niter),
    iter(0)
{}
    
    
void NniProposer::propose(Tree *tree)
{
    const float rerootProb = 1.0;
    const float doubleNniProb = 0.2;
    
    // increase iteration
    iter++;
    
    node1 = NULL;
    node3 = NULL;
    


        
    // propose new tree
    proposeRandomNni(tree, &node1, &node2, &change1);
    proposeNni(tree, node1, node2, change1);

    if (frand() < doubleNniProb) {
        proposeRandomNni(tree, &node3, &node4, &change2);        
        proposeNni(tree, node3, node4, change2);
    }

    // TODO: need random reroot or recon root.
    if (frand() < rerootProb) {
        oldroot = tree->root->children[0];
        
        if (stree != NULL) {
            reconRoot(tree, stree, gene2species);
        } else {
            /*
            // do random reroot if species tree is not given
            int choice = irand(tree->nnodes);
            tree->reroot(tree->nodes[choice]);

            int choice1 = irand(2);
            int choice2 = irand(2);
            if (tree->root->children[choice1]->nchildren == 2)
                tree->reroot(tree->root->children[choice1]->children[choice2]);
            else
                tree->reroot(tree->root->children[!choice1]->children[choice2]);
            */
        }
    } else {
        oldroot = NULL;
    }

    assert(tree->assertTree());
}

void NniProposer::revert(Tree *tree)
{
    // reject, undo topology change
    
    //printf("NNI %d %d %d %d\n", node1->name, node1->parent->name, 
    //       node2->name, node2->nchildren);
    if (node3)
        proposeNni(tree, node3, node3->parent, change2);
    if (node1)
        proposeNni(tree, node1, node1->parent, change1);

    if (oldroot)
        tree->reroot(oldroot);

}


bool NniProposer::more()
{
    return iter < niter;
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
    float predupprob, float dupprob) :
    
    nnodes(nnodes),
    stree(stree),
    params(params),
    gene2species(gene2species),
    recon(nnodes),
    events(nnodes),
    predupprob(predupprob),
    dupprob(dupprob)
{}


float SpidirBranchLikelihoodFunc::likelihood(Tree *tree) {
    // reconcile tree to species tree
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    return treelk(tree, stree,
                  recon, events, params,
                  -1, 0,
                  predupprob, dupprob, 0);
}




//=============================================================================
// MCMC search

Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                string *genes, int nseqs, int seqlen, char **seqs,
                int niter, float predupprob, float dupprob)
{
        
    NniProposer nniProposer(stree, gene2species, niter);
    ParsimonyFitter parsFitter(nseqs, seqlen, seqs);
    
    int nnodes = nseqs * 2 - 1;
    SpidirBranchLikelihoodFunc lkfunc(nnodes, stree, params, gene2species,
                                      predupprob, dupprob);
    
    return searchMCMC(initTree, 
                      genes, nseqs, seqlen, seqs,
                      &lkfunc,
                      &nniProposer,
                      &parsFitter);
}
                


Tree *searchMCMC(Tree *initTree, 
                 string *genes, int nseqs, int seqlen, char **seqs,
                 BranchLikelihoodFunc *lkfunc,
                 TopologyProposer *proposer,
                 BranchLengthFitter *fitter)
{
    Tree *toptree = NULL;
    float toplogl = -1e10, logl=-1e10, nextlogl;
    Tree *tree = NULL;
    int nnodes = nseqs * 2 - 1;
        
    
    // determine initial tree
    if (initTree != NULL) {
        // use specified initial tree
        tree = initTree;
        logl = 0.0;
    } else {
        // initialized tree with NJ
        ExtendArray<int> ptree(nnodes);
        ExtendArray<float> dists(nnodes);
        Matrix<float> distmat(nseqs, nseqs);
        
        calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
        neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);
        tree = new Tree(nnodes);
        ptree2tree(nnodes, ptree, tree);
        tree->setLeafNames(genes);
        
        // reconroot        
        // tree = phylo.reconRoot(tree, stree, gene2species)
        parsimony(tree, nseqs, seqs);
        logl = fitter->findLengths(tree);
    }
    
    
    // init likelihood score
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
        tree->hashkey(key);
        
        // TURN OFF hashing
        // have we seen this topology before?
        TreeLogl &tmp = hashtrees[key];
        if (0 && tmp.first != NULL) {
            // retrieve previously seen tree and logl
            nextlogl = tmp.second;
        } else {    
            // calculate lieklihood
            nextlogl = fitter->findLengths(tree);
            nextlogl += lkfunc->likelihood(tree);
        
            // hash result
            tmp = TreeLogl(tree->copy(), nextlogl);
        }
        
        
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
                printLog(LOG_LOW, "search: logl = %f\n", logl);            
                if (isLogLevel(LOG_LOW))
                    displayTree(toptree, getLogFile());
            
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


} // namespace spidir
