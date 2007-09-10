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
#include "phylogeny.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "ConfigParam.h"
#include "Sequences.h"
#include "spidir.h"



using namespace std;
using namespace spidir;

#define VERSION_INFO  "\
   ___    SPIDIR v0.8 (beta) July 2007 \n\
  /0 0\\   SPecies Informed DIstanced-base Reconstruction \n\
  \\___/   Matt Rasmussen \n\
 /// \\\\\\  CSAIL, MIT \n\
"


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
    int niter = 0;
    string lenfitter;
    float tsvratio;
    string bgfreqstr;
    float predupprob=0.01;
    float dupprob=1.0;
    string logfile;
    int verbose = LOG_QUIET;
    bool help = false;
    bool version = false;
    
    
    // parse arguments
    ConfigParser config;
    config.add(new ConfigParam<string>(
        "-a", "--align", "<alignment fasta>", &alignfile, 
        "sequence alignment in fasta format"));
    config.add(new ConfigParam<string>(
        "-o", "--output", "<output filename prefix>", &outprefix, "maxml",
        "prefix for all output filenames"));
    
    
    config.add(new ConfigParamComment("Sequence model evolution"));
    config.add(new ConfigParam<string>(
        "-l", "--lengths", "(hky|parsimony)", &lenfitter, "hky",
        "algorithm for determining branch lengths"));    
    config.add(new ConfigParam<float>(
        "-r", "--tsvratio", "<transition/transversion ratio>", &tsvratio, 0.5,
        "used for HKY model (default=0.5)"));
    config.add(new ConfigParam<string>(
        "-f", "--bgfreq", "<A freq>,<C ferq>,<G freq>,<T freq>", 
        &bgfreqstr, ".25,.25,.25,.25",
        "background frequencies (default=0.25,0.25,0.25,0.25"));
    
    
    config.add(new ConfigParamComment("Miscellaneous"));
    config.add(new ConfigParam<int>(
        "-i", "--niter", "<# iterations>", &niter, 100, 
        "number of iterations"));

    config.add(new ConfigParam<int>(
        "-V", "--verbose", "<verbosity level>", &verbose, LOG_QUIET, 
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
    // read sequences
    Sequences *aln;
    
    if ((aln = readAlignFasta(alignfile.c_str())) == NULL ||
        !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
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
   
    
    

    //=====================================================
    // init likelihood function
    BranchLikelihoodFunc lkfunc;
    
    // init topology proposer
    NniProposer nniProposer(NULL, NULL, niter);
    
    
    // determine branch length algorithm
    BranchLengthFitter *fitter = NULL;
    if (lenfitter == "parsimony") {
        fitter = new ParsimonyFitter(aln->nseqs, aln->seqlen, aln->seqs);
    }
    else if (lenfitter == "hky") {
        int maxiter = 3*nnodes;
        fitter = new HkyFitter(aln->nseqs, aln->seqlen, aln->seqs, 
                               bgfreq, tsvratio, maxiter);
    } else {
        printError("unknown branch length fitting algorithm: '%s'", 
                   lenfitter.c_str());
        return 1;
    }
    
    
    // search
    Tree *toptree = searchMCMC(NULL, 
                               aln->names, aln->nseqs, aln->seqlen, aln->seqs,
                               &lkfunc,
                               &nniProposer,
                               fitter);
    
    toptree->setLeafNames(aln->names);
    toptree->writeNewick(outtreeFilename.c_str());
    
    delete toptree;
    delete fitter;
    
    closeLogFile();
}


