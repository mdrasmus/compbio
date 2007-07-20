//=============================================================================
//  SPIDIR - tree search


#include "common.h"
#include "Tree.h"
#include "Matrix.h"
#include "likelihood.h"
#include "parsimony.h"
#include "search.h"

NniProposer nniProsposer;


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
        //assert(node2->nchildren >= 2);
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
    /*
    if random.random() < conf["rerootprob"]:
        nodes = tree.nodes.values()
        newnode = nodes[random.randint(0, len(nodes)-1)]
        tree2 = treelib.reroot(tree2, newnode.name)
    */
    
    // find edges for NNI
    int choice = tree->root->name;
    do {
        choice = int((rand() / float(RAND_MAX)) * tree->nnodes);
    } while (tree->nodes[choice]->isLeaf() || 
             tree->nodes[choice]->parent == NULL);
    
    *node1 = tree->nodes[choice];
    *node2 = tree->nodes[choice]->parent;
    *change = int(rand() / float(RAND_MAX) * 2);
}


NniProposer::NniProposer() :
    node1(NULL),
    node2(NULL),
    node3(NULL),
    node4(NULL)
{}
    
    

void NniProposer::propose(Tree *tree)
{
    // propose new tree
    proposeRandomNni(tree, &node1, &node2, &change1);
    proposeNni(tree, node1, node2, change1);

    if (frand() < .5) {
        proposeRandomNni(tree, &node3, &node4, &change2);        
        proposeNni(tree, node3, node4, change2);
    }

    // TODO: need random reroot or recon root.
    oldroot = tree->root->children[0];
    //int choice = int((rand() / float(RAND_MAX)) * tree->nnodes);
    //tree->reroot(tree->nodes[choice]);
    int choice1 = int(frand() * 2);
    int choice2 = int(frand() * 2);
    if (tree->root->children[choice1]->nchildren == 2)
        tree->reroot(tree->root->children[choice1]->children[choice2]);
    else
        tree->reroot(tree->root->children[!choice1]->children[choice2]);

    assert(tree->assertTree());
}

void NniProposer::revert(Tree *tree)
{
    // reject, undo topology change
    tree->reroot(oldroot);
    //printf("NNI %d %d %d %d\n", node1->name, node1->parent->name, 
    //       node2->name, node2->nchildren);
    if (node3)
        proposeNni(tree, node3, node3->parent, change2);
    proposeNni(tree, node1, node1->parent, change1);    
}



Tree *searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                int nseqs, int seqlen, char **seqs,
                int niter, 
                TopologyProposer *proposer)
{
    Tree *toptree = NULL;
    float toplogl = -1e10, logl=-1e10;
    Tree *tree = NULL;
    int nnodes = nseqs * 2 - 1;
    
    float predupprob=.001, dupprob=1.0, errorlogl=0;
    
    
    // init with NJ    
    if (initTree != NULL) {
        // need to do a copy
        tree = initTree;
    } else {
        ExtendArray<int> ptree(nnodes);
        ExtendArray<float> dists(nnodes);
        Matrix<float> distmat(nseqs, nseqs);
        
        calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
        neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);
        tree = new Tree(nnodes);
        ptree2tree(nnodes, ptree, tree);
        
        // reconroot        
        // tree = phylo.reconRoot(tree, stree, gene2species)
        
        parsimony(tree, nseqs, seqs);
    }
    
    
    // init likelihood score
    ExtendArray<int> recon(nnodes);
    ExtendArray<int> events(nnodes);
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);
    
    logl = treelk(tree, stree,
                  recon, events, params,
                  -1, 0,
                  predupprob, dupprob, errorlogl);
    toplogl = logl;
    toptree = tree->copy();
    
    float speed = 0;
    
    // MCMC loop
    for (int i=0; i<niter; i++) {
        printf("iter %d\n", i);
        Node *node1, *node2, *node3=NULL, *node4=NULL;
        int change1=0, change2=0;
    
        // propose new tree
        proposeRandomNni(tree, &node1, &node2, &change1);
        proposeNni(tree, node1, node2, change1);
        
        if (frand() < .5) {
            proposeRandomNni(tree, &node3, &node4, &change2);        
            proposeNni(tree, node3, node4, change2);
        }
        
        // TODO: need random reroot or recon root.
        Node *oldroot = tree->root->children[0];
        //int choice = int((rand() / float(RAND_MAX)) * tree->nnodes);
        //tree->reroot(tree->nodes[choice]);
        int choice1 = int(frand() * 2);
        int choice2 = int(frand() * 2);
        if (tree->root->children[choice1]->nchildren == 2)
            tree->reroot(tree->root->children[choice1]->children[choice2]);
        else
            tree->reroot(tree->root->children[!choice1]->children[choice2]);
        
        assert(tree->assertTree());
        
        parsimony(tree, nseqs, seqs);
        
        //findMLBranchLengthsHky(tree, nseqs, seqs, 
        //                       bgfreq, ratio, tree->nnodes);
        
        // calc new likelihood
        reconcile(tree, stree, gene2species, recon);
        labelEvents(tree, recon, events);
        
        float nextlogl = treelk(tree, stree,
                                recon, events, params,
                                -1, 0,
                                predupprob, dupprob, errorlogl);
        
        // acceptance rule
        //printf("%f %f\n", nextlogl - logl, log(rand() / float(RAND_MAX)));
        if (nextlogl > logl ||
            nextlogl - logl + speed > log(rand() / float(RAND_MAX)))
        {
            printf("accept %f  %f\n", nextlogl, logl);
            // accept
            logl = nextlogl;
            speed /= 2.0;

            // keep track of toptree            
            if (logl > toplogl) {
                delete toptree;
                speed = 0.0;
                toptree = tree->copy();
                toplogl = logl;
            }
        } else {
            printf("reject\n");
            //speed = (speed + 1) * 1.3;
            
            // reject, undo topology change
            tree->reroot(oldroot);
            //printf("NNI %d %d %d %d\n", node1->name, node1->parent->name, 
            //       node2->name, node2->nchildren);
            if (node3)
                proposeNni(tree, node3, node3->parent, change2);
            proposeNni(tree, node1, node1->parent, change1);
        }
    }
    
    
    return toptree;
}
