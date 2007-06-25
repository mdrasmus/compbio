#include <Python.h>
#include <string.h>

#include "spidir.h"
#include "common.h"
#include "Tree.h"


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
    if (PyTuple_GET_SIZE(args) < 13) {
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
    PyObject *pydisterror = PyTuple_GET_ITEM(args, 9);
    PyObject *pypredupprob = PyTuple_GET_ITEM(args, 10);
    PyObject *pydupprob = PyTuple_GET_ITEM(args, 11);
    PyObject *pyerrorprob = PyTuple_GET_ITEM(args, 12);
    
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
        !PyFloat_Check(pydisterror) ||
        !PyFloat_Check(pypredupprob) ||
        !PyFloat_Check(pydupprob) ||
        !PyFloat_Check(pyerrorprob)
        )
    {
        printf("wrong argument types\n");
        return NULL;
    }
    
    bool error = false;
    
    // return value
    float logl = 0;
    
    // gene tree
    int nnodes;
    int *ptree = NULL;
    float *dists = NULL;
    
    // species tree
    int nsnodes;
    int *pstree = NULL;
    
    // reconciliation
    int *gene2species = NULL;
    
    // params
    float *mu = NULL;
    float *sigma = NULL;
    float alpha = PyFloat_AS_DOUBLE(pyalpha);
    float beta = PyFloat_AS_DOUBLE(pybeta);    
    float generate = PyFloat_AS_DOUBLE(pygenerate);
    float disterror = PyFloat_AS_DOUBLE(pydisterror);
    float predupprob = PyFloat_AS_DOUBLE(pypredupprob);
    float dupprob = PyFloat_AS_DOUBLE(pydupprob);
    float errorprob = PyFloat_AS_DOUBLE(pyerrorprob);
    
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        error = true;
        goto cleanup;
    }

    if (!makeFloatArray(pydists, &dists, &nnodes)) {
        printf("bad dists\n");
        error = true;
        goto cleanup;
    }
    
    if (!makeIntArray(pypstree, &pstree, &nsnodes)) {
        printf("bad pstree\n");
        error = true;
        goto cleanup;
    }
    
    if (!makeIntArray(pygene2species, &gene2species, &nnodes)) {
        printf("bad gene2species\n");
        error = true;
        goto cleanup;
    }
    
    if (!makeFloatArray(pymu, &mu, &nsnodes)) {
        printf("bad mu\n");
        error = true;
        goto cleanup;
    }
    
    if (!makeFloatArray(pysigma, &sigma, &nsnodes)) {
        printf("bad sigma\n");
        error = true;
        goto cleanup;
    }
    
    
    /*
    // display all information
    printIntArray(ptree, nnodes);
    printFloatArray(dists, nnodes);
    printIntArray(pstree, nsnodes);
    printIntArray(recon, nnodes);
    printIntArray(events, nnodes);
    printIntArray(gene2species, nnodes);
    printFloatArray(mu, nnodes);
    printFloatArray(sigma, nnodes);
    */
    
    {
        // make tree object
        Tree tree(nnodes);
        ptree2tree(nnodes, ptree, &tree);
        tree.setDists(dists);

        SpeciesTree stree(nsnodes);
        ptree2tree(nsnodes, pstree, &stree);
        stree.setDepths();

        // reconcile gene tree to species tree
        int *recon = new int [nnodes];
        int *events = new int [nnodes];

        reconcile(&tree, &stree, gene2species, recon);
        labelEvents(&tree, recon, events);


        /*
        printf("events\n");
        printIntArray(events, nnodes);
        printIntArray(events2, nnodes);
        printf("\n");
        */

        //printTree(&tree);
        ///printTree(&stree);
        //tree2ptree(tree, ptree);


        // calculate likelihood
        logl = treelk(nnodes, ptree, dists,
                      nsnodes, pstree, 
                      recon, events,
                      mu, sigma, generate, disterror,
                      predupprob, dupprob, errorprob, alpha, beta);



        delete recon;
        delete events;
    }
    
    cleanup:
    if (ptree) delete [] ptree;
    if (dists) delete [] dists;
    if (pstree) delete [] pstree;
    if (gene2species) delete [] gene2species;
    if (mu) delete [] mu;
    if (sigma) delete [] sigma;
        
    
    if (error)
        return NULL;
    else
        return Py_BuildValue("f", logl);
}


static PyObject *
pyspidir_parsimony(PyObject *self, PyObject *args)
{
    PyObject *ret = NULL;
    bool error = false;
    
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
    int *ptree = NULL;
    char **seqs = NULL;
    int nseqs;
    
    // convert data
    if (!makeIntArray(pyptree, &ptree, &nnodes)) {
        printf("bad ptree\n");
        error = true;
        goto cleanup;
    }

    if (!makeStringArray(pyseqs, &seqs, &nseqs)) {
        printf("bad seqs\n");
        error = true;
        goto cleanup;
    }
    
    
    // call C code
    {
        //for (int i=0; i<nseqs; i++)
        //    printf("%s\n", seqs[i]);

    
        float *dists = new float [nnodes];
        for (int i=0; i<nnodes; i++) dists[i] = 0;
        parsimony(nnodes, ptree, nseqs, seqs, dists);
        ret = makeFloatListPy(dists, nnodes);
        delete [] dists;    
    }
    
    cleanup:
        if (ptree) delete [] ptree;
        if (seqs) freeStringArray(seqs, nseqs);

    
    return ret;
}


PyMODINIT_FUNC
initpyspidir(void)
{
    srand(time(NULL));
    
    static PyMethodDef methods[] = {
        {"treelk",  pyspidir_treelk, METH_VARARGS,
         "Tree likelihood"},
        {"parsimony",  pyspidir_parsimony, METH_VARARGS,
         "Parsimony method"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
    
    PyObject *m = Py_InitModule("pyspidir", methods);
}


} // extern C
