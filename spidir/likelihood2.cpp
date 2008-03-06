#include <math.h>
#include "Tree.h"
#include "phylogeny.h"


namespace spidir
{

float birthDeathDensity(float *times, int ntimes, float maxtime, 
                        float birthRate, float deathRate)
{
    const float l = birthRate;
    const float u = deathRate;
    const float r = l - u;
    const float a = u / l;
    const float T = maxtime;

    const float uT = (1.0 - exp(-r*T)) / (1.0 - a * exp(-r*T));
    const float p0 = a*uT;
    
    
    
    if (ntimes == 0) {
        return p0;
    } else if (ntimes == 1) {
        const float p1 = (1.0 - a*uT) * (1.0 - uT);
        return p1;
    } else {
        float p = (1.0 - p0);
        
        // (n - 1)! r^(n-2) (1 - a)^n
        for (int i=1; i<=ntimes-2; i++)
            p *= i * r * (1.0 - a);
        p *= (ntimes - 1) * (1.0 - a) * (1.0 - a);
        
        //printf("p_1 %f\n", p);
        
        // exp( r * sum_{i=3}^n x_i)
        float sum = 0.0;
        for (int i=2; i<ntimes; i++)
            sum += T - times[i];
        p *= exp(r * sum);
        
        // prod_{i=2}^n 1/(exp(r x_i) - a)^2
        float prod = 1.0;
        for (int i=1; i<ntimes; i++) {
            float tmp = exp(r*(T-times[i])) - a;
            prod *= (tmp*tmp);
        }
        p *= 1.0 / prod;

        // f(t_2 - t_1, 1)
        p *= (r * exp(-r * (times[1] - times[0]))) / 
             (1.0 - a * exp(-r * (T - times[0])));
        
        //printFloatArray(times, ntimes);
        //printf("p %f; T %f\n", p, T);
        return p;
    }
}


void birthDeathTreePrior_recurse(Node *node, float time, int *events, 
                                 ExtendArray<float> &times)
{
    if (events[node->name] == EVENT_DUP) {
        times.append(time += node->dist);
        
        for (int j=0; j<node->nchildren; j++) {
            birthDeathTreePrior_recurse(node->children[j], time + node->dist,
                                        events, times);
        }
    }
}


// NOTE: assumes binary species tree
float birthDeathTreePrior(Tree *tree, SpeciesTree *stree, int *recon, 
                          int *events, float birthRate, float deathRate)
{
    float prob = 0.0;

    // prepare gene tree with added implicit speciations
    ExtendArray<int> recon2(0, tree->nnodes);
    ExtendArray<int> events2(0, tree->nnodes);
    recon2.extend(recon, tree->nnodes);
    events2.extend(events, tree->nnodes);
    //Tree *tree2 = tree->copy();
    
    int addedNodes = -tree->nnodes;
    addImpliedSpecNodes(tree, stree, recon2, events2);
    addedNodes += tree->nnodes;
    
    // I will add code to handle this case later
    assert(events[tree->root->name] != EVENT_DUP);
 
    float _times[tree->nnodes];
    ExtendArray<float> times(0, tree->nnodes, _times);
    
    // loop through speciation nodes in tree2
    for (int i=0; i<tree->nnodes; i++) {
        const Node *node = tree->nodes[i];
        
        if (events2[node->name] == EVENT_SPEC) {
            if (node->nchildren == 1) {
                // calc loss prob
                const Node *snode = stree->nodes[recon2[node->name]];
                float maxtime;
                
                if (snode->children[0]->name == 
                    recon2[node->children[0]->name])
                {
                    maxtime = snode->children[1]->dist;
                } else {
                    maxtime = snode->children[0]->dist;
                }
                prob += log(birthDeathDensity(NULL, 0, maxtime, 
                                              birthRate, deathRate));
                
            }
            
            for (int j=0; j<node->nchildren; j++) {
                const float maxtime =
                    stree->nodes[recon2[node->children[j]->name]]->dist;
            
                if (events2[node->children[j]->name] != EVENT_DUP) {
                    // no apparent dup/loss  1->1
                    prob += log(birthDeathDensity(NULL, 1, maxtime, 
                                                  birthRate, deathRate));
                } else {
                    // duplication
                    times.clear();
                    times.append(0.0);
                    
                    birthDeathTreePrior_recurse(node->children[j], 0.0, events2, 
                                                times);
                    
                    //printf("maxtime %f\n", maxtime);
                    
                    prob += log(birthDeathDensity(times, times.size(),
                                               maxtime, birthRate, deathRate));
                }
            }
        }
    }
    
    removeImpliedSpecNodes(tree, addedNodes);
    
    
    // clean up
    times.detach();
    
    return prob;
}

} // namespace spidir
