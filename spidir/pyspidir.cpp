// python headers
#include <Python.h>

// c++ headers
#include <string.h>

// spidir headers
#include "common.h"
#include "ExtendArray.h"
#include "mldist.h"
#include "likelihood.h"
#include "phylogeny.h"
#include "parsimony.h"


using namespace spidir;


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


// returns a new reference
PyObject *makeFloatListPy(float *array, int size)
{
    PyObject *list = PyList_New(size);
    
    for (int i=0; i<size; i++) {
        PyObject *item = PyFloat_FromDouble(array[i]);
        PyList_SET_ITEM(list, i, item);
    }
    
    return list;
}


bool makeStringArray(PyObject *obj, char ***array, int *size)
{
    *size = PyList_GET_SIZE(obj);
    *array = new char*[*size];
    
    for (int i=0; i<*size; i++) {
        PyObject *item = PyList_GET_ITEM(obj, i);
        if (!PyString_Check(item))
            return false;
            
        char *str = PyString_AS_STRING(item);
        int len = strlen(str) + 1;
        (*array)[i] = new char [len];
        strncpy((*array)[i], str, len);
    }
    
    return true;    
}


void freeStringArray(char **array, int size)
{
    for (int i=0; i<size; i++)
        delete [] array[i];
    delete [] array;
}



// Calculate the likelihood of a tree
static PyObject *
pyspidir_treelk(PyObject *self, PyObject *args)
{
    
    // check number of args
    if (PyTuple_GET_SIZE(args) < 11) {
        printf("wrong number of args\n");
        return NULL;
    }
    
    // parse args
    PyObject *pyptree = PyTuple_GET_ITEM(args, 0);
    PyObject *pydists = PyTuple_GET_ITEM(args, 1);
    PyObject *pypstree = PyTuple_GET_ITEM(args, 2);
    PyObject *pygene2species = PyTuple_GET_ITEM(args, 3);
    PyObject *pymu = PyTuple_GET_ITEM(args, 4);
    PyObject *pysigma = PyTuple_GET_ITEM(args, 5);
    PyObject *pyalpha = PyTuple_GET_ITEM(args, 6);
    PyObject *pybeta = PyTuple_GET_ITEM(args, 7);
    PyObject *pygenerate = PyTuple_GET_ITEM(args, 8);
    PyObject *pypredupprob = PyTuple_GET_ITEM(args, 9);
    PyObject *pydupprob = PyTuple_GET_ITEM(args, 10);
    
    // check arg types
    if (!PyList_Check(pyptree) || 
        !PyList_Check(pydists) ||
        !PyList_Check(pypstree) ||
        !PyList_Check(pygene2species) ||
        !PyList_Check(pymu) ||
        !PyList_Check(pysigma) ||
        !PyFloat_Check(pyalpha) ||
        !PyFloat_Check(pybeta) ||        
        !PyFloat_Check(pygenerate) ||
        !PyFloat_Check(pypredupprob) ||
        !PyFloat_Check(pydupprob)
        )
    {
        printf("wrong argument types\n");
        return NULL;
    }
    
    
    // gene tree
    int nnodes;
    StackArray<int> ptree;
    StackArray<float> dists;
    
    // species tree
    int nsnodes;
    StackArray<int> pstree;
    
    // reconciliation
    StackArray<int> gene2species;
    
    // params
    StackArray<float> mu;
    StackArray<float> sigma;
    float alpha = PyFloat_AS_DOUBLE(pyalpha);
    float beta = PyFloat_AS_DOUBLE(pybeta);    
    float generate = PyFloat_AS_DOUBLE(pygenerate);
    float predupprob = PyFloat_AS_DOUBLE(pypredupprob);
    float dupprob = PyFloat_AS_DOUBLE(pydupprob);
    
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        return NULL;
    }

    if (!makeFloatArray(pydists, &dists, &nnodes)) {
        printf("bad dists\n");
        return NULL;
    }
    
    if (!makeIntArray(pypstree, &pstree, &nsnodes)) {
        printf("bad pstree\n");
        return NULL;
    }
    
    if (!makeIntArray(pygene2species, &gene2species, &nnodes)) {
        printf("bad gene2species\n");
        return NULL;
    }
    
    if (!makeFloatArray(pymu, &mu, &nsnodes)) {
        printf("bad mu\n");
        return NULL;
    }
    
    if (!makeFloatArray(pysigma, &sigma, &nsnodes)) {
        printf("bad sigma\n");
        return NULL;
    }
    
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);

    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();

    // reconcile gene tree to species tree
    ExtendArray<int> recon(nnodes);
    ExtendArray<int> events(nnodes);

    reconcile(&tree, &stree, gene2species, recon);
    labelEvents(&tree, recon, events);

    // calculate likelihood
    float logl = treelk(nnodes, ptree, dists,
                  nsnodes, pstree, 
                  recon, events,
                  mu, sigma, generate, 
                  predupprob, dupprob, alpha, beta);
    
    return Py_BuildValue("f", logl);
}




// Calculate the likelihood of a tree
static PyObject *
pyspidir_est_generate(PyObject *self, PyObject *args)
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
    PyObject *pygene2species = PyTuple_GET_ITEM(args, 3);
    PyObject *pymu = PyTuple_GET_ITEM(args, 4);
    PyObject *pysigma = PyTuple_GET_ITEM(args, 5);
    PyObject *pyalpha = PyTuple_GET_ITEM(args, 6);
    PyObject *pybeta = PyTuple_GET_ITEM(args, 7);
    
    // check arg types
    if (!PyList_Check(pyptree) || 
        !PyList_Check(pydists) ||
        !PyList_Check(pypstree) ||
        !PyList_Check(pygene2species) ||
        !PyList_Check(pymu) ||
        !PyList_Check(pysigma) ||
        !PyFloat_Check(pyalpha) ||
        !PyFloat_Check(pybeta)
        )
    {
        printf("wrong argument types\n");
        return NULL;
    }
    
    
    // gene tree
    int nnodes;
    StackArray<int> ptree;
    StackArray<float> dists;
    
    // species tree
    int nsnodes;
    StackArray<int> pstree;
    
    // reconciliation
    StackArray<int> gene2species;
    
    // params
    StackArray<float> mu;
    StackArray<float> sigma;
    float alpha = PyFloat_AS_DOUBLE(pyalpha);
    float beta = PyFloat_AS_DOUBLE(pybeta);    
    
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        return NULL;
    }

    if (!makeFloatArray(pydists, &dists, &nnodes)) {
        printf("bad dists\n");
        return NULL;
    }
    
    if (!makeIntArray(pypstree, &pstree, &nsnodes)) {
        printf("bad pstree\n");
        return NULL;
    }
    
    if (!makeIntArray(pygene2species, &gene2species, &nnodes)) {
        printf("bad gene2species\n");
        return NULL;
    }
    
    if (!makeFloatArray(pymu, &mu, &nsnodes)) {
        printf("bad mu\n");
        return NULL;
    }
    
    if (!makeFloatArray(pysigma, &sigma, &nsnodes)) {
        printf("bad sigma\n");
        return NULL;
    }
    
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setDists(dists);

    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();

    // reconcile gene tree to species tree
    ExtendArray<int> recon(nnodes);
    ExtendArray<int> events(nnodes);

    reconcile(&tree, &stree, gene2species, recon);
    labelEvents(&tree, recon, events);

    SpidirParams params(nsnodes, NULL, mu, sigma, alpha, beta);

    // find max a posteriori gene rate
    float generate = maxPosteriorGeneRate(&tree, &stree, recon, events, &params);
    
    return Py_BuildValue("f", generate);
}


static PyObject *
pyspidir_parsimony(PyObject *self, PyObject *args)
{
    // check number of args
    if (PyTuple_GET_SIZE(args) < 2) {
        printf("wrong number of args\n");
        return NULL;
    }
    
    // parse args
    PyObject *pyptree = PyTuple_GET_ITEM(args, 0);
    PyObject *pyseqs = PyTuple_GET_ITEM(args, 1);
    
    
    // check arg types
    if (!PyList_Check(pyptree) || 
        !PyList_Check(pyseqs)
        )
    {
        printf("wrong argument types\n");
        return NULL;
    }
    
    
    // gene tree
    int nnodes;
    StackArray<int> ptree;
    char** seqs;
    //int *ptree = NULL;
    //char **seqs = NULL;
    int nseqs;
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        //error = true;
        //goto cleanup;
        return NULL;
    }

    if (!makeStringArray(pyseqs, &seqs, &nseqs)) {
        printf("bad seqs\n");
        freeStringArray(seqs, nseqs);
        return NULL;
        //error = true;
        //goto cleanup;
    }
    


    ExtendArray<float> dists(nnodes); // = new float [nnodes];
    for (int i=0; i<nnodes; i++)
        dists[i] = 0;
    parsimony(nnodes, ptree, nseqs, seqs, dists);
    PyObject *ret = makeFloatListPy(dists, nnodes);
    
    freeStringArray(seqs, nseqs);
    
    return ret;
}


static PyObject *
pyspidir_mlhkydist(PyObject *self, PyObject *args)
{
    // check number of args
    if (PyTuple_GET_SIZE(args) < 5) {
        printf("wrong number of args\n");
        return NULL;
    }
    
    // parse args
    PyObject *pyptree = PyTuple_GET_ITEM(args, 0);
    PyObject *pyseqs = PyTuple_GET_ITEM(args, 1);
    PyObject *pybgfreq = PyTuple_GET_ITEM(args, 2);
    PyObject *pyratio = PyTuple_GET_ITEM(args, 3);
    PyObject *pymaxiter = PyTuple_GET_ITEM(args, 4);
    
    
    // check arg types
    if (!PyList_Check(pyptree) || 
        !PyList_Check(pyseqs) ||
        !PyList_Check(pybgfreq) ||
        !PyFloat_Check(pyratio) ||
        !PyInt_Check(pymaxiter))
    {
        printf("wrong argument types\n");
        return NULL;
    }
    
    
    // gene tree
    int nnodes;
    StackArray<int> ptree; // *ptree = NULL;
    char **seqs = NULL;
    StackArray<float> bgfreq; // = NULL;
    int nseqs;
    int nbases;
    float ratio = PyFloat_AS_DOUBLE(pyratio);
    int maxiter = PyInt_AS_LONG(pymaxiter);
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        return NULL;
    }

    if (!makeStringArray(pyseqs, &seqs, &nseqs)) {
        printf("bad seqs\n");
        freeStringArray(seqs, nseqs);
        return NULL;
    }
    
    if (!makeFloatArray(pybgfreq, &bgfreq, &nbases)) {
        printf("bad bgfreq\n");
        return NULL;
    }
    
    
    // call C code    
    ExtendArray<float> dists(nnodes); // = new float [nnodes];
    for (int i=0; i<nnodes; i++) 
        dists[i] = 0.0;
    float logl = findMLBranchLengthsHky(nnodes, ptree, nseqs, seqs, 
                           dists, bgfreq, ratio, maxiter,
                           true);
    PyObject *ret = makeFloatListPy(dists, nnodes);
    
    freeStringArray(seqs, nseqs);
    
    return Py_BuildValue("Nf", ret, logl);
}


static PyObject *
pyspidir_set_log(PyObject *self, PyObject *args)
{
    // check number of args
    if (PyTuple_GET_SIZE(args) < 2) {
        printf("wrong number of args\n");
        return NULL;
    }
    
    // parse args
    PyObject *pylevel = PyTuple_GET_ITEM(args, 0);
    PyObject *pyfile = PyTuple_GET_ITEM(args, 1);    
    
    // check arg types
    if (!PyInt_Check(pylevel) || 
        !PyString_Check(pyfile))
    {
        printf("wrong argument types\n");
        return NULL;
    }
    

    int level = PyInt_AS_LONG(pylevel);
    char *filename = PyString_AS_STRING(pyfile);
    

    if (!strcmp(filename, "")) {
        openLogFile(stdout);
    } else{
        if (!openLogFile(filename)) {
            printf("could not open file\n");
            return NULL;
        }
    }
    
    setLogLevel(level);
    Py_RETURN_NONE;
}


static PyObject *
pyspidir_close_log(PyObject *self, PyObject *args)
{
    closeLogFile();
    Py_RETURN_NONE;
}




PyMODINIT_FUNC
initpyspidir(void)
{
    srand(time(NULL));
    
    static PyMethodDef methods[] = {
        {"treelk",  pyspidir_treelk, METH_VARARGS,
         "Tree likelihood"},
        {"est_generate", pyspidir_est_generate, METH_VARARGS,
         "Estimates max a posteriori gene rate"},
        {"parsimony",  pyspidir_parsimony, METH_VARARGS,
         "Parsimony method"},
        {"mlhkydist", pyspidir_mlhkydist, METH_VARARGS,
         "ML estimates of branch lengths by HKY"},
        {"set_log", pyspidir_set_log, METH_VARARGS,
         "Sets the log level and output file"},
        {"close_log", pyspidir_close_log, METH_VARARGS,
         "Closes a log file"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
    
    PyObject *m = Py_InitModule("pyspidir", methods);
}


} // extern C
