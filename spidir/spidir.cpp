/*=============================================================================

    Test SPIDIR functions

=============================================================================*/

#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>


#include "common.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "Tree.h"
#include "ConfigParam.h"


using namespace std;



int main(int argc, char **argv)
{
    
    // parameters
    string alignfile;    
    string smapfile;
    string streefile;
    string paramsfile;
    int niter = 0;
    string lenfitter;
    float tsvratio;
    string bgfreqstr;
    bool help = false;
    
    
    // parse arguments
    ConfigParser config;
    config.add(new ConfigParam<string>(
        "-a", "--align", "<alignment fasta>", &alignfile, 
        "sequence alignment in fasta format"));
    config.add(new ConfigParam<string>(
        "-S", "--smap", "<species map>", &smapfile, 
        "gene to species map"));
    config.add(new ConfigParam<string>(
        "-s", "--stree", "<species tree>", &streefile, 
        "species tree file in newick format"));
    config.add(new ConfigParam<string>(
        "-p", "--param", "<spidir params file>", &paramsfile, 
        "SPIDIR branch length parameters file"));
    
    
    config.add(new ConfigParamComment("Miscellaneous"));
    config.add(new ConfigParam<int>(
        "-i", "--niter", "<# iterations>", &niter, 100, 
        "number of iterations"));
    config.add(new ConfigParam<string>(
        "-l", "--lengths", "(hky|parsimony)", &lenfitter, "hky",
        "algorithm for determining branch lengths"));
    config.add(new ConfigParam<float>(
        "-r", "--tsvratio", "<transition/transversion ratio>", &tsvratio, 0.5,
        "used for HKY model"));
    config.add(new ConfigParam<string>(
        "-f", "--bgfreq", "<A freq>,<C ferq>,<G freq>,<T freq>", 
        &bgfreqstr, ".25,.25,.25,.25",
        "background frequencies"));
    config.add(new ConfigSwitch(
        "-h", "--help", &help, "display help information"));

    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }
    
    if (help) {
        config.printHelp();
        return 0;
    }

    
    //============================================================
    // read species tree
    SpeciesTree stree;
    stree.readNewick(streefile.c_str());
    stree.setDepths();
    
    // read gene2species map
    Gene2species g;
    g.read(smapfile.c_str());
    

    // read sequences
    Sequences *aln;
    
    if ((aln = readAlignFasta(alignfile.c_str())) == NULL ||
        !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
        return 1;
    }

    // read SPIDIR parameters
    SpidirParams *params;
    if ((params = readSpidirParams(paramsfile.c_str())) == NULL)
    {
        printError("bad parameters file");
        return 1;
    }
    
    if (!params->order(&stree)) {
        printError("parameters do not correspond to the given species tree");
        return 1;
    }
    
    // determine background base frequency
    float bgfreq[4];
    vector<string> tokens = split(bgfreqstr.c_str(), ",");
    if (tokens.size() != 4) {
        printError("bgfreq requires four base frequencies e.g .25,.25,.25,.25");
        return 1;
    }
    for (unsigned int i=0; i<tokens.size(); i++)
        bgfreq[i] = atof(tokens[i].c_str());
    
    
    int nnodes = aln->nseqs * 2 - 1;

    // produce mapping array
    ExtendArray<string> genes(0, nnodes);
    ExtendArray<string> species(stree.nnodes);
    ExtendArray<int> gene2species(nnodes);
    genes.extend(aln->names, aln->nseqs);
    for (int i=aln->nseqs; i<nnodes; i++)
        genes.append("");
    stree.getLeafNames(species);   
    
    g.getMap(genes, nnodes, species, stree.nnodes, gene2species);
    
    
    // determine branch length algorithm
    BranchLengthFitter *fitter = NULL;
    if (lenfitter == "parsimony") {
        fitter = new ParsimonyFitter(aln->nseqs, aln->seqlen, aln->seqs);
    }
    else if (lenfitter == "hky") {
        int maxiter = 2*nnodes;
        fitter = new HkyFitter(aln->nseqs, aln->seqlen, aln->seqs, 
                               bgfreq, tsvratio, maxiter);
    } else {
        printError("unknown branch length fitting algorithm: '%s'", 
                   lenfitter.c_str());
        return 1;
    }
        
    
    // search
    Tree *toptree = searchMCMC(NULL, &stree,
                               params, gene2species,
                               aln->nseqs, aln->seqlen, aln->seqs,
                               niter, &nniProposer,
                               fitter);
    
    toptree->setLeafNames(genes);
    toptree->writeNewick();
    
    delete toptree;
}


int test_gene2species(int argc, char **argv)
{
    Gene2species g;
    Tree tree, stree;
    
    tree.readNewick(argv[1]);
    stree.readNewick(argv[2]);
    g.read(argv[3]);
    
    ExtendArray<string> genes(tree.nnodes);
    tree.getLeafNames(genes);

    ExtendArray<string> species(stree.nnodes);
    stree.getLeafNames(species);

    
    for (int i=0; i<tree.nnodes; i++) {
        printf("'%s' -> '%s'\n", genes[i].c_str(), 
               g.getSpecies(genes[i]).c_str());
    }
    
    ExtendArray<int> map(tree.nnodes);
    g.getMap(genes, tree.nnodes, species, stree.nnodes, map);
    
    printIntArray(map, tree.nnodes);
    
    return 0;
}
