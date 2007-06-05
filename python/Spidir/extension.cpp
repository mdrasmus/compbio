#include "Python.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#define PI 3.1415926

// events
enum {
    EVENT_GENE = 0,
    EVENT_SPEC = 1,
    EVENT_DUP = 2
};

// fractional branches
enum {
    FRAC_NONE,
    FRAC_DIFF,
    FRAC_PARENT,
    FRAC_NODE
};


// prototypes
void makeFtree(int nnodes, int *ptree, int ***ftree);
void freeFtree(int nnodes, int **ftree);
float normallog(float x, float u, float s);
void printIntArray(int *array, int size);
void printFloatArray(float *array, int size);



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
    ReconParams(int nnodes) :
        nnodes(nnodes)
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
};


// right now works for pre-order traversal
class TreeWalker
{
public:
    TreeWalker(int nnodes, int **ftree, int start) :
        nnodes(nnodes),
        node(start),        
        ftree(ftree)
    {
        stack = new int [nnodes];
        stacki = 0;
        stack[0] = start;
    }
    
    ~TreeWalker()
    {
        delete [] stack;
    }
    
    bool recurse(int node)
    {
        // descend tree
        if (ftree[node][0] != -1) {
            stack[++stacki] = ftree[node][1];
            stack[++stacki] = ftree[node][0];
            return true;
        } else
            return false;
    }
    
    int next()
    {
        // done walking
        if (stacki < 0)
            return -1;
        
        // get next node
        node = stack[stacki];
        
        // pop node of stack
        stacki--;
        
        return node;
    }
    
    
    int nnodes;
    int node;
    int **ftree;
    
    int *stack;
    int stacki;
};


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
    
    if (totvar < 1e-10) {
        printf("%f %f\n", totmean, totvar);
        printf("%d %d\n", startfrac, endfrac);
        assert(0);
    }
    
    return normallog(dist, totmean, sqrt(totvar));

    /*
    totmean = 0.0
    totvar  = 0.0
    
    #print k[node], startfrac, midparams, endfrac, endparams
    
    if startfrac == FRAC_DIFF:
        totmean += (k[node] - k[node.parent]) * startparams[0]
        totvar  += (k[node] - k[node.parent]) * startparams[1] ** 2
    elif startfrac == FRAC_PARENT:
        totmean += (1.0 - k[node.parent]) * startparams[0]
        totvar  += (1.0 - k[node.parent]) * startparams[1] ** 2
    #else startfrac == FRAC_NONE:
    #    pass
    
    if midparams != None:
        totmean += midparams[0]
        totvar  += midparams[1] ** 2
    
    if endfrac == FRAC_PARENT:
        totmean += (1.0 - k[node.parent]) * endparams[0]
        totvar  += (1.0 - k[node.parent]) * endparams[1] ** 2
    elif endfrac == FRAC_NODE:
        totmean += k[node] * endparams[0]
        totvar  += k[node] * endparams[1] ** 2
    #else endfrac == FRAC_NONE:
    #    pass
    
    if totvar <= 0.0:
        print "!!!!"
        print k[node], k[node.parent]
        print startfrac, startparams, midparams, endfrac, endparams
    
    return log(stats.normalPdf(dist, [totmean, math.sqrt(totvar)]))
    */
}


// Calculate the likelihood of a subtree
float subtreelk(int nnodes, int *ptree, int **ftree, float *dists, int root,
                int nsnodes, int *pstree, 
                int *recon, int *events,
                float *mu, float *sigma, float generate,
                int nsamples=100)
{
    float logl = 0.0;
    int sroot = nsnodes - 1;
    
    ReconParams reconparams = ReconParams(nnodes);
    
    if (events[root] != EVENT_DUP) {
        // single branch, no integration needed
        
        if (recon[root] != sroot) {
            // no midpoints, no integration needed
            //logl = normallog(dists[root] / generate,
            //                  mu[recon[root]], sigma[recon[root]]);
            
            reconBranch(root, ptree, pstree, recon, events, 
                        mu, sigma, &reconparams);
            logl = branchlk(dists[root] / generate, 
                            root, ptree, &reconparams);
            //if (logl < -10) {
            //    printf("debug: %f %f %f %f\n", dists[root], generate, 
            //                             mu[recon[root]], sigma[recon[root]]);
            //}
        }
        
        //printf("logl.s: %f\n", logl);
    } else {
        // multiple branches, integrate
        
        // set reconparams by traversing subtree
        TreeWalker walk = TreeWalker(nnodes, ftree, root);
        int node;
        int *nodes =  new int [nnodes];
        int nodesi = 0;
        
        while ((node = walk.next()) != -1) {
            //printf("node> %d %d\n", node, ptree[node]);
            reconBranch(node, ptree, pstree, 
                        recon, events, mu, sigma, &reconparams);
            nodes[nodesi++] = node;
            
            if (events[node] == EVENT_DUP) {
                walk.recurse(node);
            }
        }
        
        // choose number of samples based on number of nodes to integrate over
        nsamples = int(100*log(nodesi)) + 50;
        
        // perform integration by sampling
        float prob = 0.0;
        for (int i=0; i<nsamples; i++) {
            float sampleLogl = 0.0;
            
            // propose a setting of midpoints
            setRandomMidpoints(root, ptree, nodes, nodesi,
                               recon, events, &reconparams);
            
            // loop through all branches in subtree
            for (int j=0; j<nodesi; j++) {
                int node = nodes[j];
                
                if (recon[node] != sroot) {
                    sampleLogl += branchlk(dists[node] / generate, 
                                           node, ptree, &reconparams);
                }
            }
            
            prob += expf(sampleLogl) / nsamples;
        }
        
        logl = log(prob);
        //printf("logl.i: %f\n", logl);
        
        // cleanup
        delete [] nodes;
    }
    
    

    /*
    logl = 0.0
    midpoints[root.parent] = 1.0
    
    # integration is only needed if child is dup
    if events[child] != "dup":
        setMidpoints(child, events, recon, midpoints, [])
        logl += calcSubtreeLikelihood(child, recon, events, stree, params, 
                  midpoints, extraBranches, baserate)
    else:
        startparams = {}
        startfrac = {}
        midparams = {}
        endparams = {}
        endfrac = {}
        kdepend = {}

        # recon subtree
        nodes = []
        def walk(node):
            nodes.append(node)
            startparams[node], startfrac[node], midparams[node], \
                endparams[node], endfrac[node], kdepend[node] = \
                reconBranch2(node, recon, events, params)

            if events[node] == "dup":
                for child in node.children:
                    walk(child)
        walk(child)

        print "sample"
        samples = 100
        val = 0.0            

        for i in xrange(samples):
            setMidpointsRandom2(child, events, recon, midpoints)                
            
            val2 = 0.0
            for node in nodes:
                if recon[node] != stree.root:
                    val2 += branchLikelihood2(node.dist / baserate, 
                                  node, midpoints, 
                                  startparams[node], startfrac[node],
                                  midparams[node], endparams[node], 
                                  endfrac[node])
            val += math.exp(val2)
        logl += log(val / float(samples)) 
    */

    return logl;
}


// Calculate the likelihood of a tree
float treelk(int nnodes, int *ptree, float *dists,
             int nsnodes, int *pstree, 
             int *recon, int *events,
             float *mu, float *sigma, float generate)
{
    float logl = 0.0;
    int root = nnodes - 1;

    // make forward tree
    int **ftree;
    makeFtree(nnodes, ptree, &ftree);
    
    
    // loop through independent subtrees
    for (int i=0; i<nnodes; i++) {
        if (events[i] == EVENT_SPEC || i == root) {
            for (int j=0; j<2; j++) {
                logl += subtreelk(nnodes, ptree, ftree, dists, ftree[i][j],
                                  nsnodes, pstree, 
                                  recon, events,
                                  mu, sigma, generate);
            }
        }
    }
    
    // clean up
    freeFtree(nnodes, ftree);
    
    return logl;
}



// computes the log(normalPdf(x | u, s^2))
float normallog(float x, float u, float s)
{
    assert(s > 1e-10);
    
    //return 1.0/(s * sqrt(2.0*PI)) * exp(- (x-u)*(x-u) / (2.0 * s*s));
    return - logf(s * sqrt(2.0*PI)) - (x-u)*(x-u) / (2.0*s*s);
}


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


//=============================================================================
// debug

void printIntArray(int *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%d ", array[i]);
    printf("\n");
}

void printFloatArray(float *array, int size)
{
    for (int i=0; i<size; i++)
        printf("%f ", array[i]);
    printf("\n");
}


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}


//=============================================================================
// Python interface
extern "C" {


bool makeIntArray(PyObject *obj, int **array, int *size)
{
    *size = PyList_GET_SIZE(obj);
    *array = new int[*size];
    
    for (int i=0; i<*size; i++) {
        PyObject *item = PyList_GET_ITEM(obj, i);
        if (!PyInt_Check(item)) {
            delete [] array;
            return false;
        }
        (*array)[i] = PyInt_AS_LONG(item);
    }
    
    return true;
}

bool makeFloatArray(PyObject *obj, float **array, int *size)
{
    *size = PyList_GET_SIZE(obj);
    *array = new float[*size];
    
    for (int i=0; i<*size; i++) {
        PyObject *item = PyList_GET_ITEM(obj, i);
        if (!PyFloat_Check(item))
            return false;
        (*array)[i] = PyFloat_AS_DOUBLE(item);
    }
    
    return true;
}



// Calculate the likelihood of a tree
static PyObject *
extension_treelk(PyObject *self, PyObject *args)
{
    
    // check number of args
    if (PyTuple_GET_SIZE(args) < 8) {
        printf("wrong number of args\n");
        return NULL;
    }
    
    // parse args
    PyObject *pyptree = PyTuple_GET_ITEM(args, 0);
    PyObject *pydists = PyTuple_GET_ITEM(args, 1);
    PyObject *pypstree = PyTuple_GET_ITEM(args, 2);
    PyObject *pyrecon = PyTuple_GET_ITEM(args, 3);
    PyObject *pyevents = PyTuple_GET_ITEM(args, 4);
    PyObject *pymu = PyTuple_GET_ITEM(args, 5);
    PyObject *pysigma = PyTuple_GET_ITEM(args, 6);
    PyObject *pygenerate = PyTuple_GET_ITEM(args, 7);
    
    // check arg types
    if (!PyList_Check(pyptree) || 
        !PyList_Check(pydists) ||
        !PyList_Check(pypstree) ||
        !PyList_Check(pyrecon) ||
        !PyList_Check(pyevents) ||
        !PyList_Check(pymu) ||
        !PyList_Check(pysigma) ||
        !PyFloat_Check(pygenerate))
    {
        printf("wrong argument types\n");
        return NULL;
    }
    
    
    // gene tree
    int nnodes;
    int *ptree = NULL;
    float *dists = NULL;
    
    // species tree
    int nsnodes;
    int *pstree = NULL;
    
    // reconciliation
    int *recon = NULL;
    int *events = NULL;
    
    // params
    float *mu = NULL;
    float *sigma = NULL;
    float generate = PyFloat_AS_DOUBLE(pygenerate);
    
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        goto error;
    }

    if (!makeFloatArray(pydists, &dists, &nnodes)) {
        printf("bad dists\n");
        goto error;
    }
    
    if (!makeIntArray(pypstree, &pstree, &nsnodes)) {
        printf("bad pstree\n");
        goto error;
    }

    if (!makeIntArray(pyrecon, &recon, &nnodes)) {
        printf("bad recon\n");
        goto error;
    }
    
    if (!makeIntArray(pyevents, &events, &nnodes)) {
        printf("bad events\n");
        goto error;
    }
    
    if (!makeFloatArray(pymu, &mu, &nsnodes)) {
        printf("bad mu\n");
        goto error;
    }
    
    if (!makeFloatArray(pysigma, &sigma, &nsnodes)) {
        printf("bad sigma\n");
        goto error;
    }
    
    // error handling, cleanup
    if (0) {
        error:
            if (ptree) delete [] ptree;
            if (dists) delete [] dists;
            if (pstree) delete [] pstree;
            if (recon) delete [] recon;
            if (events) delete [] events;
            if (mu) delete [] mu;
            if (sigma) delete [] sigma;

            return NULL;
    }
    
    /*
    // display all information
    printIntArray(ptree, nnodes);
    printFloatArray(dists, nnodes);
    printIntArray(pstree, nsnodes);
    printIntArray(recon, nnodes);
    printIntArray(events, nnodes);
    printFloatArray(mu, nnodes);
    printFloatArray(sigma, nnodes);
    */
    
    // calculate likelihood
    float logl = treelk(nnodes, ptree, dists,
                        nsnodes, pstree, 
                        recon, events,
                        mu, sigma, generate);
    
    
    return Py_BuildValue("f", logl);    
}




PyMODINIT_FUNC
initextension(void)
{
    static PyMethodDef methods[] = {
        {"treelk",  extension_treelk, METH_VARARGS,
         "Tree likelihood"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
    
    PyObject *m = Py_InitModule("extension", methods);
}


} // extern C
