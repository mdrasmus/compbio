//=============================================================================
// phylogeny functions



#include "phylogeny.h"
#include "Matrix.h"
#include <utility>

namespace spidir {


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
        float low = INFINITY;
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


//=============================================================================
// reconciliation functions


typedef pair<Node*, Node*> Edge;
void getReconRootOrder(Node *node, ExtendArray<Edge> *edges)
{
    edges->append(Edge(node, node->parent));
    
    if (!node->isLeaf()) {
        for (int i=0; i<node->nchildren; i++)
            getReconRootOrder(node->children[i], edges);
        edges->append(Edge(node, node->parent));
    }
}



// NOTE: assumes binary tree
void reconRoot(Tree *tree, SpeciesTree *stree, int *gene2species)
{
    
    // determine rooting order
    ExtendArray<Edge> edges(0, tree->nnodes);
    edges.append(Edge(tree->root->children[0],
                      tree->root->children[1]));
        
    for (int i=0; i<tree->root->nchildren; i++) {
        Node *node = tree->root->children[i];
        for (int j=0; j<node->nchildren; j++) {
            getReconRootOrder(node->children[j], &edges);
        }
        edges.append(Edge(tree->root->children[0],
                          tree->root->children[1]));
    }
    
    // try initial root and recon
    tree->reroot(edges[0].first, edges[0].second);
    ExtendArray<int> recon(tree->nnodes);
    ExtendArray<int> events(tree->nnodes);
    reconcile(tree, stree, gene2species, recon);
    labelEvents(tree, recon, events);

    int minroot = 0;
    int mincost = countDuplications(events.size(), events);
    int cost = mincost;
    
    // try other roots
    for (int i=0; i<edges.size(); i++) {
        // get new edge
        Edge edge = edges[i];
        if (edge.first->parent != edge.second)
            swap(edge.first, edge.second);
    
        // uncount cost
        if (events[tree->root->name] == EVENT_DUP)
            cost--;
        if (events[edge.second->name] == EVENT_DUP)
            cost--;
        
        // reroot
        tree->reroot(edge.first, edge.second);
        
        // Recompute recon and events
        recon[edge.second->name] = reconcileNode(edge.second, stree, recon);
        recon[tree->root->name] = reconcileNode(tree->root, stree, recon);
        events[edge.second->name] = labelEventsNode(edge.second, recon);
        events[tree->root->name] = labelEventsNode(tree->root, recon);

        // count any new duplications
        if (events[tree->root->name] ==  EVENT_DUP)
            cost ++;        
        if (events[edge.second->name] ==  EVENT_DUP)
            cost++;
        
        // record mincost root
        if (cost < mincost) {
            mincost = cost;
            minroot = i;
        }
    }
    
    // root tree by minroot
    tree->reroot(edges[minroot].first, edges[minroot].second);
}



// Find Last Common Ancestor
Node *treeLca(SpeciesTree *stree, Node *node1, Node *node2)
{
    int index1 = stree->preorder[node1->name];
    int index2 = stree->preorder[node2->name];
    
    while (index1 != index2) {
        if (index1 > index2) {
            node1 = node1->parent;
            index1 = stree->preorder[node1->name];
        } else {
            node2 = node2->parent;
            index2 = stree->preorder[node2->name];
        }
    }
    
    return node1;
}


// NOTE: assumes binary species tree
void reconcile_recurse(Tree *tree, Node *node, SpeciesTree *stree, int *recon)
{
    // recurse
    for (int i=0; i<node->nchildren; i++)
        reconcile_recurse(tree, node->children[i], stree, recon);
    
    if (node->nchildren > 0) {
        int sname1 = recon[node->children[0]->name];
        int sname2 = recon[node->children[1]->name];
    
        // this node's species is lca of children species
        recon[node->name] = treeLca(stree, 
                                    stree->nodes[sname1], 
                                    stree->nodes[sname2])->name;
    }
}


// TODO: implement more efficiently with post order traversal
// reconcile a gene tree with a species tree
void reconcile(Tree *tree, SpeciesTree *stree,
               int *gene2species, int *recon)
{  
    // label gene leaves with their species
    for (int i=0; i<tree->nnodes; i++)
        if (tree->nodes[i]->isLeaf())
            recon[i] = gene2species[i];
    
    reconcile_recurse(tree, tree->root, stree, recon);    
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


int countLoss_recurse(Node *node, SpeciesTree *stree, int *recon)
{
    int loss = countLossNode(node, stree, recon);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        loss += countLoss_recurse(node->children[i], stree, recon);

    return loss;
}


int countLoss(Tree *tree, SpeciesTree *stree, int *recon)
{
    return countLoss_recurse(tree->root, stree, recon);
}


// assumes binary tree
int countLossNode(Node *node, SpeciesTree *stree, int *recon)
{
    int loss = 0;
    
    // if not parent, then no losses
    if (node->parent == NULL)
        return 0;
    
    // determine starting and ending species
    Node *sstart = stree->nodes[recon[node->name]];
    Node *send = stree->nodes[recon[node->parent->name]];
    
    // the species path is too short to have losses
    if (sstart == send)
        return 0;
    
    // determine species path of this gene branch (node, node->parent)
    Node *ptr = sstart;
    while (ptr != send) {
        // process ptr
        if (ptr != sstart)
            loss += ptr->nchildren - 1;
        
        // go up species tree
        ptr = ptr->parent;
    }
    
    // determine whether node->parent is a dup
    // if so, send (species end) is part of species path
    if (send->name == recon[node->parent->children[0]->name] ||
        send->name == recon[node->parent->children[1]->name])
         loss += send->nchildren - 1;
        
    return loss;
}

//=============================================================================
// Gene2species

const string Gene2species::NULL_SPECIES;

bool Gene2species::read(const char *filename)
{
    BufferedReader reader;
    if (!reader.open(filename, "r"))
        return false;
    
    char *line;
    string expr, species;
    char *ptr;
    
    // process each line of the file    
    while ((line = reader.readLine())) {
        // skip blank lines
        if (strlen(line) < 2)
            continue;
    
        char *e = strtok_r(line, "\t", &ptr);
        char *s = strtok_r(NULL, "\n", &ptr);
        
        // if bad format, quit
        if (e == NULL || s == NULL)
            return false;
        
        expr = e;
        species = s;
        
        if (expr.size() == 0) {
            // bad gene name expression
            return false;
        } else if (expr[0] == '*') {
            // suffix
            m_rules.append(Gene2speciesRule(Gene2speciesRule::SUFFIX,
                                            expr.substr(1, expr.size()-1), 
                                            species));
        } else if (expr[expr.size() - 1] == '*') {
            // prefix
            m_rules.append(Gene2speciesRule(Gene2speciesRule::PREFIX,
                                            expr.substr(0, expr.size()-1), 
                                            species));
        } else {
            // exact match
            m_exactLookup[expr] = species;
        }
    }

    return true;
}

string Gene2species::getSpecies(string gene)
{
    for (int i=0; i<m_rules.size(); i++) {
        switch (m_rules[i].rule) {
            case Gene2speciesRule::PREFIX:
                if (gene.find(m_rules[i].expr, 0) == 0)
                    return m_rules[i].species;
                break;

            case Gene2speciesRule::SUFFIX:
                if (gene.rfind(m_rules[i].expr, gene.size()-1) == 
                    gene.size() - m_rules[i].expr.size())
                    return m_rules[i].species;
                break;
        }
    }

    // try to find gene in exact match hashtable
    return m_exactLookup[gene];
}

// Returns the gene to species mapping in the 'map' output array
// Only leaf nodes will have defined values for gene2species
bool Gene2species::getMap(string *genes, int ngenes, 
                          string *species, int nspecies, int *map)
{
    // process each gene
    for (int i=0; i<ngenes; i++) {
        // get the species name for the gene
        string sp = getSpecies(genes[i]);

        if (sp.size() == 0) {
            // no species mapping
            // map[i] = -1;
            return false;
        } else {
            // find the index of the species name
            map[i] = -1;
            for (int j=0; j<nspecies; j++) {
                if (sp == species[j]) {
                    map[i] = j;
                    break;
                }
            }
            
            // no species found
            if (map[i] == -1)
                return false;
        }
    }

    return true;
}



} // namespace spidir
