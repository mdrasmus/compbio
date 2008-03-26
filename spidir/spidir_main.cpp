/*=============================================================================

    Test SPIDIR functions

=============================================================================*/

// c++ headers
#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>

// spidir headers
#include "common.h"
#include "phylogeny.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "ConfigParam.h"
#include "Sequences.h"
#include "spidir.h"


#define VERSION_INFO  "\
   ___    SPIDIR v0.9 (beta) 2008 \n\
  /0 0\\   SPecies Informed DIstanced-base Reconstruction \n\
  \\___/   Matt Rasmussen \n\
 /// \\\\\\  CSAIL, MIT \n\
"


using namespace std;
using namespace spidir;



int main(int argc, char **argv)
{
    // seed random number generator
    srand(time(NULL));
    
    // parameters
    string alignfile;    
    string smapfile;
    string streefile;
    string paramsfile;
    string outprefix;
    string search;
    string correctFile;
    string lkfuncopt;
    int niter = 0;
    string lenfitter;
    float tsvratio;
    string bgfreqstr;
    float predupprob = 0.01;
    float dupprob = 1.0;
    float lossprob = 1.0;
    string logfile;
    int verbose = LOG_QUIET;
    bool help = false;
    bool version = false;
    bool estGenerate = false;
    int bootiter = 1;
    
    
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
    config.add(new ConfigParam<string>(
        "-o", "--output", "<output filename prefix>", &outprefix, "spidir",
        "prefix for all output filenames"));
    
    
    config.add(new ConfigParamComment("Sequence model evolution"));
    config.add(new ConfigParam<string>(
        "-l", "--lengths", "(hky|spidir|hky_spidir|parsimony|birthdeath)", &lenfitter, "hky",
        "algorithm for determining branch lengths (default: hky)"));    
    config.add(new ConfigParam<float>(
        "-r", "--tsvratio", "<transition/transversion ratio>", &tsvratio, 1.0,
        "used for HKY model (default=1.0)"));
    config.add(new ConfigParam<string>(
        "-f", "--bgfreq", "<A freq>,<C ferq>,<G freq>,<T freq>", 
        &bgfreqstr, ".25,.25,.25,.25",
        "background frequencies (default=0.25,0.25,0.25,0.25)"));

    
    config.add(new ConfigParamComment("Miscellaneous"));
    config.add(new ConfigParam<string>(
        "", "--search", "mcmc|climb", &search, "mcmc", 
        "search algorithm"));
    config.add(new ConfigParam<int>(
        "-i", "--niter", "<# iterations>", &niter, 100, 
        "number of iterations"));
    config.add(new ConfigParam<float>(
        "-D", "--dupprob", "<duplication probability>", &dupprob, -1.0,
        "probability of a node being a duplication (default=-1.0, do not use)"));
    config.add(new ConfigParam<float>(
        "-L", "--lossprob", "<loss probability>", &lossprob, -1.0,
        "probability of loss (default=-1.0, do not use)"));
    config.add(new ConfigParam<float>(
        "-P", "--predupprob", "<pre-duplication probability>", &predupprob, 0.01,
        "probability of a node being a pre-duplication (default=0.01)"));
    config.add(new ConfigSwitch(
        "-g", "--generate", &estGenerate, "estimate generate"));
    config.add(new ConfigParam<string>(
        "-c", "--correct", "<correct tree file>", &correctFile, ""
        "check if correct tree is visited in search"));
    config.add(new ConfigParam<string>(
        "", "--lkfunc", "hky|spidir|duploss|birthdeath|none", &lkfuncopt, "spidir",
        "function for branch length likelihood"));
    config.add(new ConfigParam<int>(
        "-b", "--boot", "<# bootstraps>", &bootiter,
        "number of bootstraps to perform (default: 1)"));

    config.add(new ConfigParam<int>(
        "-V", "--verbose", "<verbosity level>", &verbose, LOG_LOW, 
        "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
    config.add(new ConfigParam<string>(
        "", "--log", "<log filename>", &logfile, "", 
        "log filename.  Use '-' to display on stdout."));
    config.add(new ConfigSwitch(
        "-v", "--version", &version, "display version information"));
    config.add(new ConfigSwitch(
        "-h", "--help", &help, "display help information"));

    
    
    if (!config.parse(argc, (const char**) argv)) {
        if (argc < 2)
            config.printHelp();
        return 1;
    }
    
    // display help
    if (help) {
        config.printHelp();
        return 0;
    }
    
    // display version info
    if (version) {
        printf(VERSION_INFO);
        return 0;
    }
    
    
    //============================================================
    // output filenames
    string outtreeFilename = outprefix  + ".tree";
    
    // use default log filename
    if (logfile == "")
        logfile = outprefix + ".log";
    
    if (logfile == "-") {
        // use standard out
        openLogFile(stdout);
    } else {
        if (!openLogFile(logfile.c_str())) {
            printError("cannot open log file '%s'.", logfile.c_str());
            return 1;
        }
    }
    
    setLogLevel(verbose);
    
    if (isLogLevel(LOG_LOW)) {
        printLog(LOG_LOW, "SPIDIR executed with the following arguments:\n");
        for (int i=0; i<argc; i++) {
            printLog(LOG_LOW, "%s ", argv[i]);
        }
        printLog(LOG_LOW, "\n\n");
    }
    
    
    //============================================================
    // read species tree
    SpeciesTree stree;
    if (!stree.readNewick(streefile.c_str())) {
        printError("error reading species tree '%s'", streefile.c_str());
        return 1;
    }
    stree.setDepths();
    
    
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
        printError("error reading parameters file '%s'", paramsfile.c_str());
        return 1;
    }
    
    if (!params->order(&stree)) {
        printError("parameters do not correspond to the given species tree");
        return 1;
    }
    
    // initialize species tree branch lengths
    float generate = params->alpha / params->beta;
    for (int i=0; i<stree.nnodes; i++)
        stree.nodes[i]->dist = params->mu[i] * generate;
    
    
    // determine background base frequency
    float bgfreq[4];
    vector<string> tokens = split(bgfreqstr.c_str(), ",");
    if (tokens.size() != 4) {
        printError("bgfreq requires four base frequencies e.g .25,.25,.25,.25");
        return 1;
    }
    for (unsigned int i=0; i<tokens.size(); i++) {
        if (sscanf(tokens[i].c_str(), "%f", &bgfreq[i]) != 1) {
            printError("bgfreq must be floats");
            return 1;
        }
    }
    
    
    int nnodes = aln->nseqs * 2 - 1;

    // read gene2species map
    Gene2species mapping;
    if (!mapping.read(smapfile.c_str())) {
        printError("error reading gene2species mapping '%s'", smapfile.c_str());
        return 1;
    }

    // produce mapping array
    ExtendArray<string> genes(0, aln->nseqs);
    genes.extend(aln->names, aln->nseqs);    
    
    ExtendArray<string> species(stree.nnodes);
    stree.getLeafNames(species);
        
    ExtendArray<int> gene2species(nnodes);
    mapping.getMap(genes, aln->nseqs, species, stree.nnodes, gene2species);
    
    
    //=====================================================
    // init likelihood function
    BranchLikelihoodFunc *lkfunc;
    
    if (lkfuncopt == "none")
        lkfunc = new BranchLikelihoodFunc();
    else if (lkfuncopt == "spidir")
        lkfunc = new SpidirBranchLikelihoodFunc(nnodes, &stree, params, 
                                                gene2species,
                                                predupprob, dupprob, lossprob,
                                                estGenerate);
    else if (lkfuncopt == "duploss")
        lkfunc = new SpidirBranchLikelihoodFunc(nnodes, &stree, params, 
                                                gene2species,
                                                predupprob, dupprob, lossprob, 
                                                estGenerate,
                                                true);
    else if (lkfuncopt == "hky") 
        lkfunc = new HkyBranchLikelihoodFunc(aln->nseqs, aln->seqlen, aln->seqs, 
                                             bgfreq, tsvratio);
    else if (lkfuncopt == "birthdeath") 
        lkfunc = new BranchLikelihoodFunc();
    else {
        printError("unknown lkfunc '%s'", lkfuncopt.c_str());
        return 1;
    }
    

    //========================================================
    // branch lengths
    
    // determine branch length algorithm
    BranchLengthFitter *fitter = NULL;
    if (lenfitter == "parsimony") {
        fitter = new ParsimonyFitter(aln->nseqs, aln->seqlen, aln->seqs);
    }
    else if (lenfitter == "hky") {
        const int maxiter = 5;
        fitter = new HkyFitter(aln->nseqs, aln->seqlen, aln->seqs, 
                               bgfreq, tsvratio, maxiter);
    } else if (lenfitter == "spidir") {
        fitter = new SpidirSample(&stree, params, gene2species);
    } else if (lenfitter == "hky_spidir") {
        const int maxiter = 1000;
        fitter = new HkySpidirSample(&stree, params, gene2species,
                                     aln->nseqs, aln->seqlen, aln->seqs, 
                                     bgfreq, tsvratio, maxiter);
    } else if (lenfitter == "birthdeath") {
        fitter = new BirthDeathFitter(aln->nseqs, aln->seqlen, aln->seqs, 
                                      bgfreq, tsvratio, 
                                      &stree, gene2species,
                                      dupprob, lossprob);
    } else {
        printError("unknown branch length fitting algorithm: '%s'", 
                   lenfitter.c_str());
        return 1;
    }
    
    
    //========================================================
    // initialize search
    
    // init topology proposer
    const int quickiter = 50;
    const float sprRatio = .3;
    SprNniProposer proposer2(&stree, gene2species, niter, sprRatio);
    DupLossProposer proposer(&proposer2, quickiter, niter);
    
    // load correct tree
    Tree correctTree;    
    if (correctFile != "") {
        if (!correctTree.readNewick(correctFile.c_str())) {
            printError("cannot read correct tree '%s'", correctFile.c_str());
            return 1;
        }
        // TODO: aborts if leaves mismatch, should catch error
        correctTree.reorderLeaves(genes);
        proposer.setCorrect(&correctTree);
    }    
    
    Tree *tree = getInitialTree(genes, aln->nseqs, aln->seqlen, aln->seqs,
                                &stree, gene2species);
    
    
    //=======================================================
    // search
    time_t startTime = time(NULL);
    Tree *toptree;
    if (search == "mcmc") {
        string mcmcfilename = outprefix + ".mcmc";
        FILE *mcmcfile = NULL;
        
        if (!(mcmcfile = fopen(mcmcfilename.c_str(), "w"))) {
            printError("cannot open mcmc file '%s'", mcmcfilename.c_str());
            return 1;
        }
        SampleFunc samples(mcmcfile);
    
        toptree = searchMCMC(tree, 
                             genes, aln->nseqs, aln->seqlen, aln->seqs,
                             &samples,
                             lkfunc,
                             &proposer,
                             fitter);
    } else if (search == "climb") {
        toptree = searchClimb(tree, 
                              genes, aln->nseqs, aln->seqlen, aln->seqs,
                              lkfunc,
                              &proposer,
                              fitter);
    } else {
        printError("unknown search '%s'", search.c_str());
        return 1;
    }
    
    //========================================================
    // output final tree
    toptree->setLeafNames(genes);
    toptree->writeNewick(outtreeFilename.c_str());
    
    
    // log tree correctness
    if (correctFile != "") {
        if (proposer.seenCorrect()) {
            printLog(LOG_LOW, "SEARCH: correct visited\n");
        } else {
            printLog(LOG_LOW, "SEARCH: correct NEVER visited\n");
        }
        
        if (toptree->sameTopology(&correctTree)) {
            printLog(LOG_LOW, "RESULT: correct\n");
        } else {
            printLog(LOG_LOW, "RESULT: wrong\n");
        }
    }
    
    
    // clean up
    delete toptree;
    delete params;
    delete fitter;
    delete lkfunc;
    
    
    // log runtime
    time_t runtime = time(NULL) - startTime;
    printLog(LOG_LOW, "runtime seconds: %d\n", runtime);
    printLog(LOG_LOW, "runtime minutes: %.1f\n", float(runtime / 60.0));
    printLog(LOG_LOW, "runtime hours: %.1f\n", float(runtime / 3600.0));
    closeLogFile();
}



