#ifndef CLUTO_TREE_H
#define CLUTO_TREE_H

// creates a permutation vector 'perm' of a matrix for best display
void ClutoReorder(
    int nrows, int ncols, int *rowptr, int *rowind, float *rowval,
    int treetype, int nparts, int *part, int *perm);


void ClutoReorderTree(
    int nrows, int ncols, int *rowptr, int *rowind, float *rowval,
    int *ptree, int *perm);

// build a parent tree for a matrix
void BuildParentTree(
    int nrows, int ncols, int *rowptr, int *rowind, float *rowval,
    int treetype, int nparts, int *part, int *ptree, int *root);

// this is used by BuildParentTree if treetype = TOP
void BuildPseudoParentTree(
    int nrows, int ncols, int *part, int *ptree, int *ptree_top, int *root);

// build a forward tree from a parent tree
void BuildForwardTree(int nrows, int *ptree, int **ftree);

// build a permutation vector from a forward tree
void ForwardTreeToPerm(int **ftree, int n, int *perm, int &index);




#endif

