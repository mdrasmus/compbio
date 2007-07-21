//=============================================================================
// SPIDIR Tree datastructure

#include <assert.h>
#include <stdio.h>


#include "Tree.h"
#include "Matrix.h"




// return a copy of the tree
Tree *Tree::copy()
{
    Tree *tree2 = new Tree(nnodes);
    Node **nodes2 = tree2->nodes;
    
    for (int i=0; i<nnodes; i++) {
        nodes2[i]->setChildren(nodes[i]->nchildren);
        nodes2[i]->name = i;
        nodes2[i]->dist = nodes[i]->dist;
    }
    
    for (int i=0; i<nnodes; i++) {
        for (int j=0; j<nodes[i]->nchildren; j++) {
            Node *child = nodes[i]->children[j];
            if (child)
                nodes2[i]->children[j] = nodes2[child->name];
            else
                nodes2[i]->children[j] = NULL;
        }
        Node *parent = nodes[i]->parent;
        if (parent)
            nodes2[i]->parent = nodes2[parent->name];
        else
            nodes2[i]->parent = NULL;
    }
    
    tree2->root = nodes2[root->name];
    
    return tree2;
}


// root tree by a new branch/node 
void Tree::reroot(Node *newroot, bool onBranch)
{
    // handle trivial case, newroot is root
    if (root == newroot ||
        (onBranch &&
         root->nchildren == 2 &&
         (root->children[0] == newroot ||
          root->children[1] == newroot)))
        return;
            
    
    // determine where to stop ascending
    Node *oldroot = root;
    Node *stop1=NULL, *stop2=NULL;
    
    if (isRooted()) {
        stop1 = root->children[0];
        stop2 = root->children[1];
    } else {
        stop1 = root;
    }

    // start the reversal
    Node *ptr1 = NULL, *ptr2 = NULL;
    float nextDist = 0;
    float rootdist;
    
    if (onBranch) {
        if (isRooted()) {
            // just need to stick current root somewhere else
            Node *other = newroot->parent;            
            rootdist = stop1->dist + stop2->dist;
            
            oldroot->children[0] = newroot;
            oldroot->children[1] = other;            
            newroot->parent = oldroot;
            newroot->dist /= 2.0;
            
            ptr1 = other;

            int oldchild = findval(ptr1->children, ptr1->nchildren, newroot);
            assert(oldchild != -1);
            
            // prepare for reversing loop
            ptr1->children[oldchild] = oldroot;
            ptr2 = oldroot;
            nextDist = newroot->dist;
        } else {
            // need to add a new node to be root
            assert(0);
        }
    } else {
        if (isRooted()) {
            // need to remove the root node, and make tribranch
            assert(0);
        } else {
            // just need to swap node positions
            assert(0);
        }
    }
    
    
    // reverse parent child relationships
    while (ptr1 != stop1 && ptr1 != stop2) {
        int oldchild = findval(ptr1->children, ptr1->nchildren, ptr2);
        assert(oldchild != -1);
        
        Node *next = ptr1->parent;
        
        // ptr1 is now fixed
        ptr1->children[oldchild] = next;
        ptr1->parent = ptr2;
        
        // swap distances
        float tmpdist = ptr1->dist;
        ptr1->dist = nextDist;
        nextDist = tmpdist;
        
        // move pointers
        ptr2 = ptr1;
        ptr1 = next;
    }
    
    
    // handle last two nodes
    if (stop2 != NULL) {
        // make stop1 parent of stop2
        if (stop2 == ptr1) {        
            Node *tmp = stop1;
            stop1 = ptr1;
            stop2 = tmp;
        }
        assert(ptr1 == stop1);
        
        int oldchild = findval(stop1->children, stop1->nchildren, ptr2);        
        stop1->children[oldchild] = stop2;
        stop1->parent = ptr2;
        stop1->dist = nextDist;
        stop2->parent = stop1;
        stop2->dist = rootdist;
    } else {
        assert(0);
    }
    
    
    // renumber nodes
    // - all leaves don't change numbers
    assert(root->name = nnodes-1);
}


void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node)
{
    if (!node)
        node = tree->root;

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreePostOrder(tree, nodes, node->children[i]);
    
    // record post-process
    nodes->append(node);
}


void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node)
{
    if (!node)
        node = tree->root;

    // record pre-process
    nodes->append(node);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreePostOrder(tree, nodes, node->children[i]);
}



//=============================================================================
// Tree Input/Output


Node *Tree::readNode(FILE *infile, Node *parent, int &depth)
{
    char chr, char1;
    Node *node;
    string token;

    // read first character
    if (!(char1  = readChar(infile, depth))) {
        printError("unexpected end of file");
        return NULL;
    }
    
    
    if (char1 == '(') {
        // read internal node
    
        int depth2 = depth;
        node = addNode(new Node());
        if (parent)
            parent->addChild(node);
        
        // read all child nodes at this depth
        while (depth == depth2) {
            Node *child = readNode(infile, node, depth);
            if (!child)
                return NULL;
        }
        
        // read distance for this node
        char chr = readUntil(infile, token, "):,;", depth);
        if (chr == ':') {
            node->dist = readDist(infile, depth);
            if (!(chr = readUntil(infile, token, "):,", depth)))
                return NULL;
        }
        //if (chr == ';' && depth == 0)
        //    return node;
        
        return node;
    } else {
        // read leaf
        
        node = addNode(new Node());
        if (parent)
            parent->addChild(node);
        
        // read name
        if (!(chr = readUntil(infile, token, ":),", depth)))
            return NULL;
        token = char1 + trim(token.c_str());
        node->leafname = token;
        
        // read distance for this node
        if (chr == ':') {
            node->dist = readDist(infile, depth);
            if (!(chr = readUntil(infile, token, ":),", depth)))
                return NULL;
        }
        
        return node;
    }
}


int nodeNameCmp(const void *_a, const void *_b)
{
    Node *a = *((Node**) _a);
    Node *b = *((Node**) _b);
    
    if (a->nchildren == 0) {
        if (b->nchildren == 0)
            return 0;
        else
            return -1;
    } else {
        if (b->nchildren == 0)
            return 1;
        else
            return 0;
    }
}


bool Tree::readNewick(FILE *infile)
{
    int depth = 0;
    root = readNode(infile, NULL, depth);
    
    // renumber nodes in a valid order
    // 1. leaves come first
    // 2. root is last
    
    if (root != NULL) {
        nodes[root->name] = nodes[nnodes-1];    
        nodes[nnodes-1] = root;
        nodes[root->name]->name = root->name;
        root->name = nnodes - 1;
        
        qsort((void*) nodes.get(), nodes.size(), sizeof(Node*),
               nodeNameCmp);
        
        // update names
        for (int i=0; i<nnodes; i++)
            nodes[i]->name = i;
        
        return true;
    } else
        return false;
}


bool Tree::readNewick(const char *filename)
{
    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'\n", filename);
        return false;
    }

    bool ret = readNewick(infile);
    
    if (!ret) {
        printError("format error in '%s'\n", filename);
    }
    
    fclose(infile);
    
    return ret;
}


// write out the newick notation of a tree
void Tree::writeNewick(FILE *out, Node *node, int depth)
{
    if (node == NULL) {
        assert(root != NULL);
        writeNewick(out, root, 0);
        fprintf(out, ";\n");
    } else {
        if (node->nchildren == 0) {
            for (int i=0; i<depth; i++) fprintf(out, "  ");
            fprintf(out, "%s:%f", node->leafname.c_str(), node->dist);
        } else {
            // indent
            for (int i=0; i<depth; i++) fprintf(out, "  ");
            fprintf(out, "(\n");
            
            for (int i=0; i<node->nchildren - 1; i++) {
                writeNewick(out, node->children[i], depth+1);
                fprintf(out, ",\n");
            }
            
            writeNewick(out, node->children[node->nchildren-1], depth+1);
            fprintf(out, "\n");
            
            for (int i=0; i<depth; i++) fprintf(out, "  ");
            fprintf(out, ")");
            
            if (depth > 0)
                fprintf(out, ":%f", node->dist);
        }
    }
}


bool Tree::writeNewick(const char *filename)
{
    FILE *out = NULL;
    
    if ((out = fopen(filename, "w")) == NULL) {
        printError("cannot write file '%s'\n", filename);
        return false;
    }

    writeNewick(out);
    return true;
}


//=============================================================================
// Visualization


void drawLine(Matrix<char> &matrix, char chr, int x1, int y1, int x2, int y2)
{
    float stepx, stepy;
    int steps;
    
    steps = max(abs(x2-x1), abs(y2-y1));
    
    stepx = float(x2 - x1) / steps;
    stepy = float(y2 - y1) / steps;
    
    float x = x1;
    float y = y1;
    for (int i=0; i<steps+1; i++) {
        matrix[int(y)][int(x)] = chr;
        x += stepx;
        y += stepy;
    }
}



// TODO: finish display
void displayNode(Matrix<char> &matrix, Node *node, int *xpt, int* ypt)
{
    //drawLine(matrix, '-', x, y, x2, y);
    int x = xpt[node->name];
    int y = ypt[node->name];
    

    
    for (int i=0; i<node->nchildren; i++) {
        displayNode(matrix, node->children[i], xpt, ypt);
    }
    
    //if (x >= 0 && y >= 0 && x < matrix.numCols() && y < matrix.numRows()) {
        // horizontal line
        if (node->parent) {
            drawLine(matrix, '-', x, y, 
                     xpt[node->parent->name], y);
        }
        
        // vertical line
        if (!node->isLeaf()) {
            int l = node->nchildren - 1;
            drawLine(matrix, '|', x,
                                  ypt[node->children[0]->name],
                                  x,
                                  ypt[node->children[l]->name]);
            matrix[ypt[node->children[0]->name]][x] = '/';
            matrix[ypt[node->children[l]->name]][x] = '\\';

            matrix[y][x] = '+';                    
        }
    //}
}


int treeLayout(Node *node, ExtendArray<int> &xpt, ExtendArray<int> &ypt,
               float xscale, int yscale, int y=0)
{
    if (node->parent == NULL) {
        xpt[node->name] = int(xscale * node->dist);
        ypt[node->name] = 0;
    } else {
        xpt[node->name] = xpt[node->parent->name] + int(xscale * node->dist + 1);
        ypt[node->name] = y;
    }
    
    if (node->isLeaf()) {
        y += yscale;
    } else {
        for (int i=0; i<node->nchildren; i++) {
            y = treeLayout(node->children[i], xpt, ypt, xscale, yscale, y);
        }
        int l = node->children[node->nchildren-1]->name;
        ypt[node->name] = (ypt[l] + ypt[node->children[0]->name]) / 2;
    }
    
    return y;
}


// TODO: finish display
void displayTree(Tree *tree, FILE *outfile, float xscale, int yscale)
{
    ExtendArray<int> xpt(tree->nnodes);
    ExtendArray<int> ypt(tree->nnodes);
    treeLayout(tree->root, xpt, ypt, xscale, yscale);

    int width = 0;    
    for (int i=0; i<tree->nnodes; i++) {
        if (xpt[i] > width) 
            width = xpt[i];
    }
    
    Matrix<char> matrix(tree->nnodes+1, width+1);
    matrix.setAll(' ');
    
    displayNode(matrix, tree->root, xpt, ypt);
    
    // write out matrix
    for (int i=0; i<matrix.numRows(); i++) {
        for (int j=0; j<matrix.numCols(); j++)
            fprintf(outfile, "%c", matrix[i][j]);
        fprintf(outfile, "\n");
    }
}


//=============================================================================
// assert that the tree datastructure is self-consistent
bool Tree::assertTree()
{
    if (root == NULL) return false;
    if (nnodes != nodes.size()) return false;
    if (root->parent != NULL) return false;
    if (root->name != nnodes - 1) return false;
    
    bool leaves = true;
    for (int i=0; i<nnodes; i++) {
        //printf("assert %d\n", i);
        if (nodes[i] == NULL) return false;
        
        // names are correct
        if (nodes[i]->name != i) return false;
        
        // do leaves come first 
        if (nodes[i]->isLeaf()) {
            if (!leaves)
                return false;
        } else
            leaves = false;
        
        // check parent child pointers
        for (int j=0; j<nodes[i]->nchildren; j++) {
            //printf("assert %d %d\n", i, j);
            if (nodes[i]->children[j] == NULL) return false;
            //printf("assert %d %d parent\n", i, j);
            if (nodes[i]->children[j]->parent != nodes[i]) return false;
        }
    }
    
    //printf("done\n");
    
    return true;
}


//=============================================================================
// phylogeny functions


// Neighbor-joining algorithm
void neighborjoin(int ngenes, float **distmat, int *ptree, float *branches)
{
    Matrix<float> dists(ngenes*2-1, ngenes*2-1);
    float *restdists = new float [ngenes*2-1];
    int *leaves = new int [ngenes];
    int nleaves = ngenes;
    int newnode = ngenes;
    
    // initialize distances
    for (int i=0; i<ngenes; i++) {
        float r = 0.0;
        for (int j=0; j<ngenes; j++) {
            dists[i][j] = distmat[i][j];
            r += distmat[i][j];
        }
        restdists[i] = r / (ngenes - 2);
    }
    
    // initialize leaves
    for (int i=0; i<ngenes; i++)
        leaves[i] = i;
    
    
    // join loop
    while (nleaves > 2) {
        // search for closest genes
        float low = MAX_FLOAT;
        int lowi = -1, lowj = -1;
        
        for (int i=0; i<nleaves; i++) {
            for (int j=i+1; j<nleaves; j++) {
                int gene1 = leaves[i];
                int gene2 = leaves[j];
                float dist = dists[gene1][gene2] - restdists[gene1] 
                                                 - restdists[gene2];
                if (dist < low) {
                    low = dist;
                    lowi = i;
                    lowj = j;
                }
            }
        }
        
        // join gene1 and gene2
        int lowgene1 = leaves[lowi];
        int lowgene2 = leaves[lowj];
        int parent = newnode++;
        ptree[lowgene1] = parent;
        ptree[lowgene2] = parent;
        
        // set distances
        branches[lowgene1] = (dists[lowgene1][lowgene2] + 
                              restdists[lowgene1] - 
                              restdists[lowgene2]) / 2.0;
        branches[lowgene2] = dists[lowgene1][lowgene2] - branches[lowgene1];
        
        // gene1 and gene2 are no longer leaves, remove them from leaf set
        leaves[lowi] = parent;
        leaves[lowj] = leaves[nleaves-1];
        nleaves--;
        
        float r = 0;
        for (int i=0; i<nleaves; i++) {
            int gene = leaves[i];
            if (gene != parent) {
                float v = (dists[lowgene1][gene] + 
                           dists[lowgene2][gene] -
                           dists[lowgene1][lowgene2]) / 2.0;
                dists[parent][gene] = v;
                dists[gene][parent] = v;
                r += v;
            }
        }
        
        if (nleaves > 2)
            restdists[parent] = r / (nleaves - 2);
    }
    
    // join the last two genes, split the remaining dist evenly
    int gene1 = leaves[0];
    int gene2 = leaves[1];
    int parent = newnode++;
    
    ptree[gene1] = parent;
    ptree[gene2] = parent;
    ptree[parent] = -1;
    branches[gene1] = dists[gene1][gene2] / 2.0;
    branches[gene2] = dists[gene1][gene2] / 2.0;
    branches[parent] = 0.0;
    
    assert(parent == ngenes*2-2);
    
    delete [] restdists;
    delete [] leaves;
}



void reconRoot(Tree *tree, SpeciesTree *stree, int *gene2species)
{
    
}


/*
def reconRoot(gtree, stree, gene2species = gene2species, 
               rootby = "duploss", newCopy=True):
    # make a consistent unrooted copy of gene tree
    if newCopy:
        gtree = gtree.copy()
    treelib.unroot(gtree, newCopy=False)
    treelib.reroot(gtree, 
                   gtree.nodes[util.sort(gtree.leafNames())[0]].parent.name, 
                   onBranch=False, newCopy=False)
    
    
    # make recon root consistent for rerooting tree of the same names
    # TODO: there is the possibility of ties, they are currently broken
    # arbitrarily.  In order to make comparison of reconRooted trees with 
    # same gene names accurate, hashOrdering must be done, for now.
    hashOrderTree(gtree, gene2species)
    
    # get list of edges to root on
    edges = []
    def walk(node):
        edges.append((node, node.parent))
        if not node.isLeaf():
            node.recurse(walk)
            edges.append((node, node.parent))
    for child in gtree.root.children:
        walk(child)
    
    
    # try initial root and recon    
    treelib.reroot(gtree, edges[0][0].name, newCopy=False)
    recon = reconcile(gtree, stree, gene2species)
    events = labelEvents(gtree, recon)     
    
    # find reconciliation that minimizes loss
    minroot = edges[0]
    rootedge = sorted(edges[0])
    if rootby == "dup": 
        cost = countDup(gtree, events)
    elif rootby == "loss":
        cost = len(findLoss(gtree, stree, recon))
    elif rootby == "duploss":
        cost = countDupLoss(gtree, stree, recon, events)
    else:
        raise "unknown rootby value '%s'"  % rootby
    mincost = cost
    
    
    # try rooting on everything
    for edge in edges[1:-1]:
        if sorted(edge) == rootedge:
            continue
        rootedge = sorted(edge)
        
        node1, node2 = edge
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2, "%s %s" % (node1.name, node2.name)
        
        # uncount cost
        if rootby in ["dup", "duploss"]:
            if events[gtree.root] == "dup":
                cost -= 1
            if events[node2] == "dup":
                cost -= 1
        if rootby in ["loss", "duploss"]:
            cost -= len(findLossNode(gtree.root, recon))
            cost -= len(findLossNode(node2, recon))
        
        # new root and recon
        treelib.reroot(gtree, node1.name, newCopy=False)        
        
        recon[node2] = reconcileNode(node2, stree, recon)
        recon[gtree.root] = reconcileNode(gtree.root, stree, recon)
        events[node2] = labelEventsNode(node2, recon)
        events[gtree.root] = labelEventsNode(gtree.root, recon)
        
        if rootby in ["dup", "duploss"]:
            if events[node2] ==  "dup":
                cost += 1
            if events[gtree.root] ==  "dup":
                cost += 1
        if rootby in ["loss", "duploss"]:
            cost += len(findLossNode(gtree.root, recon))
            cost += len(findLossNode(node2, recon))
        
        #print edge[0].name, edge[1].name, cost
        
        # keep track of min cost
        if cost < mincost:
            mincost = cost
            minroot = edge
    
    # root tree by minroot
    if edge != minroot:
        node1, node2 = minroot
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2
        treelib.reroot(gtree, node1.name, newCopy=False)
    
    return gtree
*/



// Find Last Common Ancestor
Node *treeLca(SpeciesTree *stree, Node *node1, Node *node2)
{
    int depth1 = stree->depths[node1->name];
    int depth2 = stree->depths[node2->name];
        
    // get nodes to same depth
    if (node1 != node2) {
        while (depth1 > depth2) {
            node1 = node1->parent;
            depth1 = stree->depths[node1->name];
        }
        
        while (depth2 > depth1) {
            node2 = node2->parent;
            depth2 = stree->depths[node2->name];
        }
    }
    
    // walk up both nodes until they meet
    while (node1 != node2) {
        node1 = node1->parent;
        node2 = node2->parent;
    }
    
    return node1;
}


// NOTE: assumes binary species tree
void reconcile_helper(Tree *tree, Node *node, SpeciesTree *stree, int *recon)
{
    // recurse
    for (int i=0; i<node->nchildren; i++)
        reconcile_helper(tree, node->children[i], stree, recon);
    
    if (node->nchildren > 0) {
        int sname1 = recon[node->children[0]->name];
        int sname2 = recon[node->children[1]->name];
    
        // this node's species is lca of children species
        recon[node->name] = treeLca(stree, 
                                    stree->nodes[sname1], 
                                    stree->nodes[sname2])->name;
    }
}


// reconcile a gene tree with a species tree
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon)
{  
    // label gene leaves with their species
    for (int i=0; i<tree->nnodes; i++)
        if (tree->nodes[i]->nchildren == 0)
            recon[i] = gene2species[i];
    
    reconcile_helper(tree, tree->root, stree, recon);    
}


// label events for each node in tree
// NOTE: assumes binary gene tree
void labelEvents(Tree *tree, int *recon, int *events)
{
    Node **nodes = tree->nodes;

    for (int i=0; i<tree->nnodes; i++) {
        if (nodes[i]->nchildren == 0)
            events[i] = EVENT_GENE;
        else 
        if (recon[i] == recon[nodes[i]->children[0]->name] ||
            recon[i] == recon[nodes[i]->children[1]->name])
            events[i] = EVENT_DUP;
        else
            events[i] = EVENT_SPEC;
    }
}



//=============================================================================
// primitive tree format conversion functions

// creates a forward tree from a parent tree
void makeFtree(int nnodes, int *ptree, int ***ftree)
{
    *ftree = new int* [nnodes];
    int **ftree2 = *ftree;
    
    // initialize
    for (int i=0; i<nnodes; i++) {
        ftree2[i] = new int [2];
        ftree2[i][0] = -1;
        ftree2[i][1] = -1;
    }
    
    // populate
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];
        
        if (parent != -1) {
            if (ftree2[parent][0] == -1)
                ftree2[parent][0] = i;
            else
                ftree2[parent][1] = i;
        }
    }
}


void freeFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++)
        delete [] ftree[i];
    delete [] ftree;
}


// create a tree object from a parent tree array
void ptree2tree(int nnodes, int *ptree, Tree *tree)
{
    Node **nodes = tree->nodes;
    
    // allocate children
    for (int i=0; i<nnodes; i++) {
        nodes[i]->allocChildren(2);
        nodes[i]->name = i;
        nodes[i]->nchildren = 0;
    }
    
    // store parent and child pointers
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];
        
        if (parent != -1) {
            Node *parentnode = nodes[parent];            
            parentnode->children[parentnode->nchildren++] = nodes[i];
            nodes[i]->parent = parentnode;
        } else {
            nodes[i]->parent = NULL;
        }
    }
    
    // set root
    tree->root = nodes[nnodes - 1];
    assert(tree->assertTree());
}


// create a tree object from a parent tree array
void tree2ptree(Tree *tree, int *ptree)
{
    Node **nodes = tree->nodes;
    int nnodes = tree->nnodes;
    
    for (int i=0; i<nnodes; i++) {
        if (nodes[i]->parent)
            ptree[i] = nodes[i]->parent->name;
        else
            ptree[i] = -1;
    }
}


//=============================================================================
// primitive input/output


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}


// write out the names of internal nodes
void printTree(Tree *tree, Node *node, int depth)
{
    if (node == NULL) {
        if (tree->root != NULL) {
            printTree(tree, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (node->nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=%s:%f", node->name, node->leafname.c_str(), node->dist);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=(\n", node->name);
            
            for (int i=0; i<node->nchildren - 1; i++) {
                printTree(tree, node->children[i], depth+1);
                printf(",\n");
            }
            
            printTree(tree, node->children[node->nchildren-1], depth+1);
            printf("\n");
            
            for (int i=0; i<depth; i++) printf("  ");
            printf(")");
            
            if (depth > 0)
                printf(":%f", node->dist);
        }
    }
}
