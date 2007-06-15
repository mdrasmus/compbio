#include <Python.h>
#include "spidir.h"
#include "common.h"


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
    
    if (!makeIntArray(pygene2species, &gene2species, &nnodes)) {
        printf("bad gene2species\n");
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
            if (gene2species) delete [] gene2species;
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
    printIntArray(gene2species, nnodes);
    printFloatArray(mu, nnodes);
    printFloatArray(sigma, nnodes);
    */
    
    
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
    printf("recons\n");
    printIntArray(recon, nnodes);
    printIntArray(recon2, nnodes);
    printf("\n");
    
    printf("events\n");
    printIntArray(events, nnodes);
    printIntArray(events2, nnodes);
    printf("\n");
    */
    
    //printTree(&tree);
    ///printTree(&stree);
    //tree2ptree(tree, ptree);
    
    
    // calculate likelihood
    float logl = treelk(nnodes, ptree, dists,
                        nsnodes, pstree, 
                        recon, events,
                        mu, sigma, generate, disterror,
                        predupprob, dupprob, errorprob, alpha, beta);
    

    
    delete recon;
    delete events;
    
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
