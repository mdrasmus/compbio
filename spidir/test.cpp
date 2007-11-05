/*=============================================================================

    Test SPIDIR functions

=============================================================================*/

// c++ headers
#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

// spidir headers
#include "common.h"
#include "mldist.h"
#include "likelihood.h"
#include "phylogeny.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "ConfigParam.h"
#include "Sequences.h"


using namespace std;
using namespace spidir;




int test_reconstruct(int argc, char **argv)
{    
    // parameters
    string alignfile;    
    string smapfile;
    string streefile;
    string paramsfile;
    int niter = 0;
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
    
    
    config.add(new ConfigParamComment("Miscellaneous options"));
    config.add(new ConfigParam<int>(
        "-i", "--niter", "<# iterations>", &niter, 100, 
        "number of iterations"));
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
    
    // search
    Tree *toptree = searchMCMC(NULL, &stree,
                               params, gene2species,
                               genes, aln->nseqs, aln->seqlen, aln->seqs,
                               niter, 0.01, 1.0);
    
    toptree->setLeafNames(genes);
    toptree->writeNewick();
    
    delete toptree;
    
    return 0;
}


int test_mledist(int argc, char **argv)
{
    // parameters
    string alignfile;
    
    // parse arguments
    ConfigParser config;
    config.add(new ConfigParam<string>(
        "-a", "--align", "<alignment fasta>", &alignfile, 
        "sequence alignment in fasta format"));
    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }
    
    
    // read sequences
    Sequences *aln;
    
    if ((aln = readAlignFasta(alignfile.c_str())) == NULL ||
        !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
        return 1;
    }
    
    int nnodes = aln->nseqs * 2 - 1;
    int nseqs = aln->nseqs;
    int seqlen = aln->seqlen;
    char **seqs = aln->seqs;
    
    ExtendArray<int> ptree(nnodes);
    ExtendArray<float> dists(nnodes);
    Matrix<float> distmat(nseqs, nseqs);

    calcDistMatrix(nseqs, seqlen, seqs, distmat.getMatrix());
    neighborjoin(nseqs, distmat.getMatrix(), ptree, dists);
    Tree tree(nnodes);
    ptree2tree(nnodes, ptree, &tree);
    tree.setLeafNames(aln->names);
    
    parsimony(&tree, nseqs, seqs);
    displayTree(&tree, stdout, 100);
    
    float bgfreq[4] = {.25, .25, .25, .25};
    float ratio = .5;
    int maxiter = 20;
    
    
    
    findMLBranchLengthsHky(&tree, nseqs, seqs, bgfreq, ratio, maxiter);
    displayTree(&tree, stdout, 100);
    
    return 0;
}


int test_hky(int argc, char **argv)
{
    float tsvratio, time;
    float bgfreq[4];
    string bgfreqstr;

    ConfigParser config;
    config.add(new ConfigParam<float>(
        "-r", "--tsvratio", "<transition/transversion ratio>", &tsvratio, 0.5,
        "used for HKY model (default=0.5)"));
    config.add(new ConfigParam<float>(
        "-t", "--time", "<time units>", &time, 1.0,
        "used for HKY model (default=1.0)"));
    config.add(new ConfigParam<string>(
        "-f", "--bgfreq", "<A freq>,<C ferq>,<G freq>,<T freq>", 
        &bgfreqstr, ".25,.25,.25,.25",
        "background frequencies (default=0.25,0.25,0.25,0.25"));
    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }
    
    
    // determine background base frequency
    vector<string> tokens = split(bgfreqstr.c_str(), ",");
    if (tokens.size() != 4) {
        printError("bgfreq requires four base frequencies e.g .25,.25,.25,.25");
        return 1;
    }
    for (unsigned int i=0; i<tokens.size(); i++)
        bgfreq[i] = atof(tokens[i].c_str());

    float matrix[16];
    makeHkyMatrix(bgfreq, tsvratio, time, matrix);
    
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            printf("%f ", matrix[4*i+j]);
        }
        printf("\n");
    }
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




int test_reroot(int argc, char **argv)
{
    Tree tree;
    
    tree.readNewick(argv[1]);
    
    char filename[100];
    
    for (int i=0; i<tree.nnodes; i++) {
        tree.reroot(tree.nodes[i]);
        snprintf(filename, 100, "%d.reroot.tree", i);
        tree.writeNewick(filename);
    }
}


int test_reconroot(int argc, char **argv)
{
    Tree tree;
    tree.readNewick(argv[1]);
    
    SpeciesTree stree;
    stree.readNewick(argv[2]);
    stree.setDepths();
    
    // read gene2species map
    Gene2species g;
    g.read(argv[3]);
    
    // produce mapping array
    ExtendArray<string> genes(tree.nnodes);
    ExtendArray<string> species(stree.nnodes);
    tree.getLeafNames(genes);
    stree.getLeafNames(species);
    
    ExtendArray<int> gene2species(tree.nnodes);
    g.getMap(genes, tree.nnodes, species, stree.nnodes, gene2species);
    
    reconRoot(&tree, &stree, gene2species);
    displayTree(&tree, stdout, 100);
}


int test_genbranches(int argc, char **argv)
{
    string treefile;
    string streefile;
    string smapfile;
    string paramsfile;

    // parse arguments
    ConfigParser config;
    config.add(new ConfigParam<string>(
        "-t", "--tree", "<tree topology file>", &treefile, 
        "topology to simulate"));
    config.add(new ConfigParam<string>(
        "-S", "--smap", "<species map>", &smapfile, 
        "gene to species map"));
    config.add(new ConfigParam<string>(
        "-s", "--stree", "<species tree>", &streefile, 
        "species tree file in newick format"));
    config.add(new ConfigParam<string>(
        "-p", "--param", "<spidir params file>", &paramsfile, 
        "SPIDIR branch length parameters file"));
    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }


    Tree tree;
    tree.readNewick(treefile.c_str());
    
    SpeciesTree stree;
    stree.readNewick(streefile.c_str());
    stree.setDepths();
    
    // read gene2species map
    Gene2species g;
    g.read(smapfile.c_str());
    
    // read SPIDIR parameters
    SpidirParams* params;
    if ((params = readSpidirParams(paramsfile.c_str())) == NULL)
    {
        printError("bad parameters file");
        return 1;
    }
    
    if (!params->order(&stree)) {
        printError("parameters do not correspond to the given species tree");
        return 1;
    }    
    
    
    // produce mapping array
    ExtendArray<string> genes(tree.nnodes);
    ExtendArray<string> species(stree.nnodes);
    tree.getLeafNames(genes);
    stree.getLeafNames(species);
    
    ExtendArray<int> gene2species(tree.nnodes);
    g.getMap(genes, tree.nnodes, species, stree.nnodes, gene2species);
    
    
    // reconcile gene tree to species tree
    ExtendArray<int> recon(tree.nnodes);
    ExtendArray<int> events(tree.nnodes);

    reconcile(&tree, &stree, gene2species, recon);
    labelEvents(&tree, recon, events);
     
    // generate branch lengths
    generateBranchLengths(&tree,
                          &stree,
                          recon, events,
                          params);
    
    displayTree(&tree, stdout, 100);
    tree.writeNewick(stdout);
    
    
    //for (int i=0; i<10000; i++)
    //    printf("%f\n", gammavariate(10, 2));
        //printf("%f\n", normalvariate(10, 1));
}



int main(int argc, char **argv)
{
    srand(time(NULL));

    if (argc < 2) {
        printf("choose a test:\n"
               "  reconstruct\n"
               "  gene2species\n"
               "  mledist\n"
               "  reroot\n"
               "  reconroot\n"
               "  genbranches\n");
        return 1;
    }

    string testname = argv[1];
    
    if (testname == "reconstruct") {
        test_reconstruct(argc-1, &argv[1]);
        
    } else if (testname == "gene2species") {
        test_gene2species(argc-1, &argv[1]);
        
    } else if (testname == "mledist") {
        test_mledist(argc-1, &argv[1]);
    } else if (testname == "hky") {
        test_hky(argc-1, &argv[1]);
    } else if (testname == "reroot") {
        test_reroot(argc-1, &argv[1]);
    } else if (testname == "reconroot") {
        test_reconroot(argc-1, &argv[1]);
    } else if (testname == "genbranches") {
        test_genbranches(argc-1, &argv[1]);
    } else {
        printf("unknown test\n");
        return 1;
    }

    
    return 0;
}
