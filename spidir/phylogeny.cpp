//=============================================================================
// phylogeny functions



#include "phylogeny.h"
#include "Matrix.h"


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
    while ((line = reader.readLine())) {
        //chomp(line);

        expr = strtok_r(line, "\t", &ptr);
        species = strtok_r(NULL, "\n", &ptr);

        if (expr[0] == '*') {
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
            //assert(0);
        }
    }

    return false;
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

            case Gene2speciesRule::EXACT:
                break;
        }
    }

    return NULL_SPECIES;
}

bool Gene2species::getMap(string *genes, int ngenes, 
                          string *species, int nspecies, int *map)
{
    for (int i=0; i<ngenes; i++) {
        string sp = getSpecies(genes[i]);

        if (sp.size() == 0) {
            map[i] = -1;
        } else {
            map[i] = -1;
            for (int j=0; j<nspecies; j++) {
                if (sp == species[j])
                    map[i] = j;
            }
        }
    }

    return true;
}





/*

def makeGene2species(maps):
    # find exact matches and expressions
    exacts = {}
    exps = []
    for mapping in maps:
        if "*" not in mapping[0]:
            exacts[mapping[0]] = mapping[1]
        else:
            exps.append(mapping)
    
    # create mapping function
    def gene2species(gene):
        # eval expressions first in order of appearance
        for exp, species in exps:
            if exp[-1] == "*":
                if gene.startswith(exp[:-1]):
                    return species
            elif exp[0] == "*":
                if gene.endswith(exp[1:]):
                    return species
        
        if gene in exacts:
            return exacts[gene]
        
        raise Exception("Cannot map gene '%s' to any species" % gene)
    return gene2species


*/

} // namespace spidir
