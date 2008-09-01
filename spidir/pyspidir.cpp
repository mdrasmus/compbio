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



#define def_py2array(func_name, T, checkType, convertType) \
bool func_name(PyObject *obj, T **array, int *size) \
{ \
    *size = PyList_GET_SIZE(obj); \
    *array = new T[*size]; \
    \
    for (int i=0; i<*size; i++) { \
        PyObject *item = PyList_GET_ITEM(obj, i); \
        if (!checkType(item)) { \
            delete [] array; \
            return false; \
        } \
        (*array)[i] = convertType(item); \
    } \
    \
    return true; \
}

def_py2array(makeIntArray, int, PyInt_Check, PyInt_AS_LONG);
def_py2array(makeFloatArray, float, PyFloat_Check, PyFloat_AS_DOUBLE);
def_py2array(py2array, int, PyInt_Check, PyInt_AS_LONG);
def_py2array(py2array, float, PyFloat_Check, PyFloat_AS_DOUBLE);

bool py2array(PyObject *obj, char ***array, int *size)
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

#define makeStringArray py2array



// returns a new reference
#define def_array2py(func_name, T, convertType) \
PyObject *func_name(T *array, int size) \
{ \
    PyObject *list = PyList_New(size); \
    for (int i=0; i<size; i++) { \
        PyObject *item = convertType(array[i]); \
        PyList_SET_ITEM(list, i, item); \
    } \
    return list; \
}

def_array2py(makeFloatListPy, float, PyFloat_FromDouble);
def_array2py(array2py, float, PyFloat_FromDouble);
def_array2py(array2py, int, PyInt_FromLong);



void freeStringArray(char **array, int size)
{
    for (int i=0; i<size; i++)
        delete [] array[i];
    delete [] array;
}


bool ParsePy(PyObject *args, const char *fmt, ...)
{
    va_list ap;   
    
    bool status = true;
    int *d;
    int **dl;
    float *f;
    float **fl;
    char ***cl;
    void **v;
    string *str;
    bool *b;
    
    va_start(ap, fmt);

    if (PyTuple_GET_SIZE(args) < (int) strlen(fmt))
    {
        printf("wrong number of arguments\n");
        return false;
    }
    
    // loop through format string
    int i=0;
    for (const char *argtype = fmt; *argtype && status; argtype++, i++) {        
        // get next arg in lst
        PyObject *arg = PyTuple_GET_ITEM(args, i);
        
        switch (*argtype) {
            case 'd':
                if (!PyInt_Check(arg)) {
                    printf("expected integer for argument %d", i);
                    status = false;
                    break;
                }
                
                d = va_arg(ap, int *);
                *d = PyInt_AS_LONG(arg);
                break;
            case 'f':
                if (!PyFloat_Check(arg)) {
                    printf("expected float for argument %d", i);
                    status = false;
                    break;
                }
                
                f = va_arg(ap, float *);
                *f = float(PyFloat_AS_DOUBLE(arg));
                break;
            case 'D':
                if (!PyList_Check(arg)) {
                    printf("expected list for argument %d", i);
                    status = false;
                    break;
                }
                
                dl = va_arg(ap, int **); 
                d = va_arg(ap, int *); 
                if (!makeIntArray(arg, dl, d)) {
                    printf("error parsing int list\n");
                    status = false;
                }
                break;
            case 'F':
                if (!PyList_Check(arg)) {
                    printf("expected list for argument %d", i);
                    status = false;
                    break;
                }
                
                fl = va_arg(ap, float **); 
                d = va_arg(ap, int *); 
                if (!makeFloatArray(arg, fl, d)) {
                    printf("error parsing float list\n");
                    status = false;
                }
                break;
            case 's':
                if (!PyString_Check(arg)) {
                    printf("expected string for argument %d", i);
                    status = false;
                    break;
                }
                
                str = va_arg(ap, string *);
                *str = PyString_AS_STRING(arg);
                break;
            case 'C':
                if (!PyList_Check(arg)) {
                    printf("expected list for argument %d", i);
		    status = false;
		    break;
		}
		
		cl = va_arg(ap, char ***);
		d = va_arg(ap, int *);
                if (!makeStringArray(arg, cl, d)) {
                    printf("error parsing char matrix");
		    status = false;
                }
		break;
            case 'b':
                b = va_arg(ap, bool *);
                *b = (arg != Py_False);
                break;
            case 'v':
                v = va_arg(ap, void **);
                *v = (void*) arg;
                break;
        }
    }
    
    va_end(ap);
    
    return status;
}


//=============================================================================
// Python interface
extern "C" {


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

    float lossprob = dupprob;

    // calculate likelihood
    float logl = treelk(nnodes, ptree, dists,
                  nsnodes, pstree, 
                  recon, events,
                  mu, sigma, generate, 
                  predupprob, dupprob, lossprob, alpha, beta);
    
    return Py_BuildValue("f", logl);
}


// Calculate the likelihood of a tree
static PyObject *
pyspidir_genbranches(PyObject *self, PyObject *args)
{
    // gene tree
    int nnodes;
    StackArray<int> ptree;

    // species tree
    int nsnodes;
    StackArray<int> pstree;
    
    // reconciliation
    StackArray<int> gene2species;
    
    // params
    StackArray<float> mu;
    StackArray<float> sigma;
    float alpha;
    float beta;
        

    if (!ParsePy(args, "DDDFFff", 
                 &ptree, &nnodes,
                 &pstree, &nsnodes,
                 &gene2species, &nnodes,
                 &mu, &nsnodes,
                 &sigma, &nsnodes,
                 &alpha, &beta))
        return NULL;    
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);

    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();

    // reconcile gene tree to species tree
    ExtendArray<int> recon(nnodes);
    ExtendArray<int> events(nnodes);

    reconcile(&tree, &stree, gene2species, recon);
    labelEvents(&tree, recon, events);

    // generate branch lengths
    ExtendArray<float> dists(nnodes);
    for (int i=0; i<nnodes; i++) 
        dists[i] = 0.0;
    generateBranchLengths(nnodes, ptree, 
                          nsnodes, pstree,
                          recon, events,
                          mu, sigma,
                          alpha, beta,
                          dists);
    
    PyObject *ret = makeFloatListPy(dists, nnodes);
    return ret;
}


// Calculate the likelihood of a tree
static PyObject *
pyspidir_est_generate(PyObject *self, PyObject *args)
{
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
    float alpha;
    float beta;
    
    
    if (!ParsePy(args, "DFDDFFff", 
                 &ptree, &nnodes,
                 &dists, &nnodes,
                 &pstree, &nsnodes,
                 &gene2species, &nnodes,
                 &mu, &nsnodes,
                 &sigma, &nsnodes,
                 &alpha, &beta))
        return NULL;   
    
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


void samples_posterior_gene_rate(float generate, Tree *tree, void *userdata)
{
    PyObject *callback = (PyObject*) userdata;

    PyObject *args = PyTuple_New(1);
    PyObject *pygenerate = PyFloat_FromDouble(generate);
    PyTuple_SET_ITEM(args, 0, pygenerate);
    
    Py_DECREF(PyObject_CallObject(callback, args));
    Py_DECREF(args);
}


// Calculate the likelihood of a tree
static PyObject *
pyspidir_sample_gene_rate(PyObject *self, PyObject *args)
{
    int nsamples;

    // gene tree
    int nnodes;
    StackArray<int> ptree;
    
    // species tree
    int nsnodes;
    StackArray<int> pstree;
    
    // reconciliation
    StackArray<int> gene2species;
    
    // params
    StackArray<float> mu;
    StackArray<float> sigma;
    float alpha;
    float beta;

    // sequence model params
    int nbgfreq;
    StackArray<float> bgfreq;
    float tsvratio;

    // sequence matrix
    int nseqs;
    char **seqs;

    PyObject *pycallback;
    
    
    if (!ParsePy(args, "dDDDFFffFfCv", 
                 &nsamples,
                 &ptree, &nnodes,
                 &pstree, &nsnodes,
                 &gene2species, &nnodes,
                 &mu, &nsnodes,
                 &sigma, &nsnodes,
                 &alpha, &beta,
                 &bgfreq, &nbgfreq,
                 &tsvratio,
		 &seqs, &nseqs,
                 &pycallback))
        return NULL;


    Py_INCREF(pycallback);
    
    // make tree object
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);

    // assert binary tree
    for (int i=0; i<tree.nnodes; i++) {
	assert(tree.nodes[i]->nchildren <= 2);
    }

    SpeciesTree stree(nsnodes);
    ptree2tree(nsnodes, pstree, &stree);
    stree.setDepths();

    SpidirParams params(nsnodes, NULL, mu, sigma, alpha, beta);


    samplePosteriorGeneRate(&tree, 
                            nseqs, seqs, 
                            bgfreq, tsvratio,
                            &stree, 
                            gene2species,
                            &params,
                            nsamples,
                            samples_posterior_gene_rate,
                            (void*) pycallback);
    
    Py_DECREF(pycallback);

    Py_RETURN_NONE;
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


// Calculate the likelihood of a tree
static PyObject *
pyspidir_hkymatrix(PyObject *self, PyObject *args)
{
    // args
    int nbgfreq;
    StackArray<float> bgfreq;
    float tsvratio;
    float time;
    
    // parse args
    if (!ParsePy(args, "Fff", 
                 &bgfreq, &nbgfreq,
                 &tsvratio, &time))
        return NULL;    

    if (nbgfreq != 4) {
        printError("expected 4 base frequencies");
    }
    
    // compute matrix
    float matrix[16];
    makeHkyMatrix(bgfreq, tsvratio, time, matrix);
    
    // create matrix of floats
    PyObject *pymatrix = PyList_New(4);
    
    for (int i=0; i<4; i++) {
        PyObject *row = PyList_New(4);
        PyList_SET_ITEM(pymatrix, i, row);
    
        for (int j=-0; j<4; j++) {            
            PyObject *val = PyFloat_FromDouble(matrix[matind(4, i, j)]);
            PyList_SET_ITEM(row, j, val);
        }
    }
    
    return pymatrix;
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
        {"genbranches",  pyspidir_genbranches, METH_VARARGS,
         "Generate branch lengths"},
        {"est_generate", pyspidir_est_generate, METH_VARARGS,
         "Estimates max a posteriori gene rate"},
        {"sample_gene_rate", pyspidir_sample_gene_rate, METH_VARARGS,
         "Samples gene rates from the posterior distribution"},
        {"parsimony",  pyspidir_parsimony, METH_VARARGS,
         "Parsimony method"},
        {"mlhkydist", pyspidir_mlhkydist, METH_VARARGS,
         "ML estimates of branch lengths by HKY"},
        {"hkymatrix", pyspidir_hkymatrix, METH_VARARGS,
         "compute the transition matrix for the HKY model"},
        {"set_log", pyspidir_set_log, METH_VARARGS,
         "Sets the log level and output file"},
        {"close_log", pyspidir_close_log, METH_VARARGS,
         "Closes a log file"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
    };
    
    PyObject *m = Py_InitModule("pyspidir", methods);
}


} // extern C
