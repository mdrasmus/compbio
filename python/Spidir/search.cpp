//=============================================================================
//  SPIDIR - tree search


#include "common.h"
#include "Tree.h"
#include "Matrix.h"
#include "parsimony.h"


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
        assert(node2->nchildren >= 2);
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



void proposeRandomNni(Tree *tree)
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
    
    proposeNni(tree, tree->nodes[choice], tree->nodes[choice]->parent, 
               int(rand() / float(RAND_MAX) * 2));
}



void searchMCMC(Tree *initTree, SpeciesTree *stree,
                SpidirParams *params, int *gene2species,
                int nseqs, int seqlen, char **seqs,
                int niter=500)
{
    Tree *toptree = NULL;
    float toplogl = -1e-10;
    int iter = 0;
    Tree *tree = NULL;
    
    // init with NJ    
    if (initTree != NULL) {
        // need to do a copy
        tree = initTree;
    } else {
        int nnodes = nseqs * 2 - 1;
        int *ptree = new int [nnodes];
        float *dists = new float [nnodes];
        Matrix<float> distmat(nseqs, nseqs);
        
        calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
        neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);
        tree = new Tree(nnodes);
        ptree2tree(nnodes, ptree, tree);
        
        // reconroot        
        // tree = phylo.reconRoot(tree, stree, gene2species)
        
        parsimony(tree, nseqs, seqs);
    }
    
    /*
    float predupprob=.001, dupprob=1.0, errorlogl=0;
    
    // init likelihood score
    toplogl = treelk(tree, stree,
                     recon, events, params,
                     -1, 0,
                     predupprob, dupprob, errorlogl);
    toptree = tree->copy();
    
    
    for (int i=0; i<niter; i++) {
        proposeRandomNni(tree);
        parsimony(tree, nseqs, seqs);
        float logl = treelk(tree, stree,
                            recon, events, params,
                            -1, 0,
                            predupprob, dupprob, errorlogl);
        
    }
    */
    
    /*    
    // proposal function
    def propose(chain, tree):
        tree2 = proposeFunc(conf, tree,  distmat, labels, 
                            stree, gene2species, params, visited)
        
        # check visited dict
        thash = phylo.hashTree(tree2)
        if thash in visited:
            logl, tree2, count = visited[thash]
            this.nold += 1
        else:
            Spidir.setTreeDistances(conf, tree2, distmat, labels)
            logl = Spidir.treeLogLikelihood(conf, tree2, stree, gene2species, params)
            
            #tree2, logl = Spidir.treeLogLikelihoodAllRoots(conf, tree2, stree, 
            #                                        gene2species, params)
            this.nold = 0
        
        addVisited(conf, visited, tree2, gene2species, thash)
        
        # best yet tree
        if logl > this.toplogl:
            printMCMC(conf, "%d:%d" % (chain.name, this.iter), 
                      tree2, stree, gene2species, visited)
            this.toplogl = logl
            this.toptree = tree2.copy()
            
            # move some other chains to best state
            #chains2 = sorted(chains, key=lambda x: x.logl)
            #for chain in chains2[:1]:
            #    chain.state = this.toptree.copy()
            #    chain.logl = this.toplogl
        
        # alter logl to influence search only
        chain.relax = conf["speedup"] * this.nold
               
        return tree2, logl
        
    # init chains    
    chains = []
    for i in range(conf["nchains"]):
        chains.append(McmcChain(i, tree.copy(), this.toplogl, propose))
    
    
    # run chains
    for i in xrange(1, conf["iters"]):
        this.iter += 1
        
        for chain in chains:
            chain.step()   
   

    return this.toptree, this.toplogl
    */
}
