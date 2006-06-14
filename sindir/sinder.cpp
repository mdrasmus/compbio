/*******************************************************************************
 SINDER - Matt Rasmussen
 sinder.cpp
 11/19/05

*******************************************************************************/

#include "Tree.h"

using namespace Sinder;


int main(int argc, char **argv)
{
    Tree tree;
    FILE *infile = fopen("test/test.tree", "r");
    
    if (tree.Read(infile)) {
        tree.Write(stdout);
    }
    
    fclose(infile);
    

    return 0;
}


