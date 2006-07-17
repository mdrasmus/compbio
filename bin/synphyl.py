#!/usr/bin/env python

import sys, time, os
from rasmus import synteny, synphylweb, bionj, ensembl, phyloutil, env
from rasmus import fasta, util, genomeio, genomeutil
from rasmus import muscle, phylip, algorithms, clustalw

options = [
    ["c:", "config=", "config", "<config file>"],
    ["s", "synteny", "synteny", "", 
        {"help": "cluster genes by synteny"}],
    ["R", "norefine", "norefine", "", 
        {"help": "don't refine syneny"}],
    ["p", "pull", "pull", "", 
        {"help": "finish clustering"}],
    ["P", "nopull", "nopull", "", 
        {"help": "don't use pulling partition"}],
    ["a", "align", "align", "", 
        {"help": "build alignments from clustering"}],
    ["t", "tree", "tree", "", 
        {"help": "build trees from alignments"}],
    ["l", "altsubtrees", "altsubtrees", "", 
        {"help": "count alternative species subtrees"}],
    ["L", "alttrees", "alttrees", "", 
        {"help": "count alternative species trees"}],
    ["h", "homology", "homology", "", 
        {"help": "find homology in trees"}],
    ["e", "estimate", "estimate", ""],
    ["o", "orthsets", "orthsets", "", 
        {"help": "build ortholog sets"}],
    ["v", "vissyn", "vissyn", "", 
        {"help": "build synteny visualizations"}],
    ["w", "web", "web", "", 
        {"help": "generate synphyl website"}],
    ["f", "showconfig", "showconfig", "", 
        {"help": "show current configuration"}],
    ["q", "quick", "quick", "", 
        {"help": "do quick processing"}],
    ["r:", "report=", "report", "<report type>", 
        {"help": "print a report"}],
    ["P:",  "paths=", "paths", "<data path>"]
]



param = util.parseOptions(sys.argv, options, quit=True)



# parse arguments
if "config" in param:
    configfile = param["config"][-1]
else:
    configfile = "config.py"


requiredConfig = [
    "root_name",
    "genomes"
]


defaultConfig = {
    "paths": [],
    
    "root_name": "all",
    "output_dir": ".",
    
    # matches defaults
    "match_cols": [0,1,11],
    "match_ext": "blastp",
    
    # synteny defaults
    "synteny_resume": False,
    "synteny_refine": False,
    "synteny_sigscore": 500,
    
    "usepull": True,
    "quick": False,
    "quick_minpartsize": 20,
    "save_parttrees": True,
    "boot_iters": 1,
    
    # alignment defaults
    "align_start": 0,
    "align_end": None,
    "align_fast": True,
    "align_prefix": "",
    "align_suffix": ".align",
    "align_bigtree": 300,
    
    # tree building defaults
    "tree_start": 0,
    "tree_end": None,
    "tree_boot_suffix": ".boot.tree",
    
    # vissynteny
    "vissynteny_context": 50e3,
    "vissynteny_image_size": (800, 400)    
}



def setDefaultConfig(config, param):   
    env.addPaths(* config["paths"])
    
    # genome files
    # automatically find genome files
    if "genome_files" not in config:
        config["genome_files"] = {}
        for genome in config["genomes"]:
            config["genome_files"][genome] = env.findFile(genome + ".coord")
    
    # species map
    if "smap" in config:
        config["gene2species"] = \
                    genomeutil.readGene2species(env.findFile(config["smap"]))

    # default loading of protein sequences
    #if "protein_seqs" not in config:
    #    seqs = []
    #    for genome in config["genomes"]:
    #        seqs.append(genomeio.simplepepfile(genome))
    #    config["protein_seqs"] = seqs

    config["output_dir"] = ensureProperPath(config["output_dir"])
    
    # setup output paths
    outputDirs = [
        ["align_dir", "alignments/"],
        ["synteny_dir", "synteny/"],
        ["report_coarse_dir", "report_coarse/"],
        ["report_fine_dir", "report_fine/"],
        ["vissynteny_dir", "vissynteny/"],
        ["homology_dir", "homology/"]
    ]
    
    
    # calculate full paths
    for key, path in outputDirs:    
        config.setdefault(key, config["output_dir"] + path)
        config[key] = ensureProperPath(config[key])
        ensurePathExist(config[key])
    
    
    # ensure web paths are proper
    if "web" in param:
        config["local_www_dir"] = ensureProperPath(config["local_www_dir"])
        config["remote_www_dir"] = ensureProperPath(config["remote_www_dir"])
    
        # assert that www dir exists and make symbolic link to it from output dir
        config.setdefault("wwwbase_dir", config["output_dir"] + "wwwbase/")
        assert os.path.isdir(config["local_www_dir"])
        if not os.path.isdir(config["wwwbase_dir"]):
            cwd = os.getcwd()
            os.chdir(config["output_dir"])
            os.symlink(config["local_www_dir"], "wwwbase")
            os.chdir(cwd)
        
    
    




def generateConf(config):
    assert "genconf" in config, "ERROR"
    out = file(config["genconf"], "w")
    
    print >>out, "genomes = []"
    
    
    out.close()


def showConfig(config):
    keys = config.keys()
    keys.sort()
    
    keys.remove("__builtins__")
    
    util.log("configuration")
    util.log("-------------")
    for key in keys:
        util.log("  %s = %s" % (key, str(config[key])))
    util.log()


#
# filename formats
#

def local2remote(config, path):
    return os.path.realpath(path).replace(config["local_www_dir"], 
                                          config["remote_www_dir"])


def matchFile(config, genome1, genome2):
    return env.findFile("%s_%s.%s" % (genome1, genome2, config["match_ext"]))

def synPartFile(config):
    if config["synteny_refine"]:
        return "%s.refine.syncomps" % (config["synteny_dir"] + config["root_name"])
    else:
        return "%s.syncomps" % (config["synteny_dir"] + config["root_name"])

def pulledPartFile(config):
    if config["usepull"]:
        return config["synteny_dir"] + config["root_name"] + ".pull.part"
    else:
        return synPartFile(config)

def syntenyFile(config):
    if config["synteny_refine"]:
        return config["synteny_dir"] + config["root_name"] + ".refine.synteny"
    else:
        return config["synteny_dir"] + config["root_name"] + ".synteny"

def syntenyStatsFile(config):
    return config["synteny_dir"] + config["root_name"] + ".synstats"

def alignFile(config, partid):
    return "%s%s%d%s" % (config["align_dir"], config["align_prefix"], partid,
                         config["align_suffix"])

def bootTreeFile(config, partid):
    return "%s%s%d%s" % (config["align_dir"], config["align_prefix"], partid,
                         config["tree_boot_suffix"])

def alignPublicFile(config, partid):
    return local2remote(config, alignFile(config, partid))
    
def bootTreePublicFile(config, partid):
    return  local2remote(config, bootTreeFile(config, partid))

def speciesTreeFile(config):
    return config["species_tree"]

def homologyFile(config):
    return config["homology_dir"] + config["root_name"] + ".homology"

def orthologSetFile(config):
    return config["homology_dir"] + config["root_name"] + ".homology.part"

def visSyntenyIndexFile(config):
    return config["vissynteny_dir"] + "vissynteny.index"

def altTreesFile(config):
    return config["align_dir"] + config["root_name"] + ".alttrees.html"

def altSubtreesFile(config):
    return config["align_dir"] + config["root_name"] + ".altsubtrees.html"

def fastaFile(config, genome):
    return env.findFile("%s.fasta" % genome)

def logFile(config):
    return config["output_dir"] + "synphyl.log"



def ensureProperPath(path):
    if not path.endswith("/"):
        path += "/"
    return path

def ensurePathExist(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise "cannot create %s" % path


def readProteinSeqs(config):
    seqs = {}
    for f in config["protein_seqs"]:
        seqs.update(fasta.readFasta(f, valuefunc=fasta.removestar))
    return seqs


def readGeneDescriptions(config):
    # load descriptions
    desc = {}
    for genome in config["genomes"]:
        if os.path.isfile(genomeio.descfile(genome)):
            desc.update(genomeio.readGeneDesc(genomeio.descfile(genome)))
    return desc


def readGenomes(config):
    matching = genomeutil.Matching()
    for genome, genomefile in config["genome_files"].items():
        util.tic("reading '%s'" % genome)
        matching.readGenomes(genomefile, config["gene2species"])
        util.toc()
    matching.autoconf()
    return matching
    

def setupSyntenyInputs(config):    
    # setup input data
    data = []
    
    if "match_pairs" in config:    
        for genome1, genome2 in config["match_pairs"]:
            data.append([genome1, genome2, 
                         config["genome_files"][genome1], 
                         config["genome_files"][genome2],
                         "",    # fasta file not needed
                         "",    # fasta file not needed
                         matchFile(config, genome1, genome2)])
    else:
        for genome1 in config["genomes"]:
            for genome2 in config["genomes"]:
                if genome1 == genome2: continue
                try:
                    data.append([genome1, genome2, 
                         config["genome_files"][genome1], 
                         config["genome_files"][genome2],
                         "",    # fasta file not needed
                         "",    # fasta file not needed
                         matchFile(config, genome1, genome2)])
                except env.PathError:
                    pass

    
    
    # setup configuration
    conf = synteny.initConf()

    mapping = [
        ["cols", "match_cols"],
        ["minscore", "synteny_minscore"],
        ["minBlockSize", "synteny_minblocksize"],
        ["justBbh", "synteny_just_bbh"],
        ["sigScore", "synteny_sigscore"],
        ["sigBbhScore", "synteny_sigbbhscore"],
        ["nearRange", "synteny_near"],
        ["rangeGenes", "synteny_range"],
        ["pullFrac", "synteny_pullfrac"],
        ["minPullLen", "synteny_minpulllen"],
        ["minGeneSize", "synteny_mingenesize"],
        ["relScore", "synteny_relscore"]
    ]
    
    for key in config:
        conf[key] = config[key]
    
    for name1, name2 in mapping:
        if name2 in config:
            conf[name1] = config[name2]
    
    return data, conf


def runSynteny(config):
    util.tic("SYNPHYL SYNTENY")
    data, conf = setupSyntenyInputs(config)
    
    # find synteny
    if not config["synteny_refine"]:
        comps = synteny.findMultiSynteny(
            data, conf, config["synteny_dir"], config["root_name"])
    else:
        comps = util.readDelim(config["synteny_refine"])
    
    # do final statistics
    os.system("syntenystats.py -c %s -s %s -g %s -S %s -P %s > %s" % 
              (synPartFile(config),
               syntenyFile(config),
               ",".join(conf["genomes"]),
               conf["smap"],
               ":".join(conf["paths"]),
               syntenyStatsFile(config)))
    
    util.toc()    

"""
def runGenePulling(config):
    util.tic("SYNPHYL GENE PULLING")
    data, conf = setupSyntenyInputs(config)
    matching = readGenomes(config)
    genes = matching.getGenes()
    comps = util.readDelim(synPartFile(config))
    
    synteny.pullGenes(genes, comps, data, conf, config["logfile"])

    # save partitions with pulled genes
    util.writeDelim(pulledPartFile(config), comps)
    
    util.toc()
"""


def runAlignments(config):
    util.tic("SYNPHYL ALIGNMENTS")

    fastas = ""
    for genome in config["genomes"]:
        fastas += "-f %s.fasta " % genome

    os.system("partalign.py -p %s -o %s -t %s" % 
              (synPartFile(config), config["align_dir"], fastas))
    
    util.toc()
    
    
"""
def runAlignments2(config):
    util.tic("SYNPHYL ALIGNMENTS")
    
    # read inputs
    seqs = readProteinSeqs(config)
    parts = util.readDelim(pulledPartFile(config))
    
    print pulledPartFile(config)
    
    if "align_inds" not in config:
        if config["align_end"] == None:
            partinds = xrange(config["align_start"], len(parts))
        else:
            partinds = xrange(config["align_start"], config["align_end"])
    else:
        partinds = config["align_inds"]
    
    
    # estimation
    if "estimate" in config:
        n = 0
        for i in xrange(len(parts)):
            part = parts[i]
            # decide if tree partitioning is need
            if config["quick"] and not treePartNeeded(config, part):
                continue
            n += 1
        print "will create %d alignments" % n
        util.toc()
        return
    
    
    for i in partinds:
        part = parts[i]
        util.tic("part %d of %d, size %d" % (i, len(parts), len(part)))
        seqs2 = util.subdict(seqs, part)

        if len(seqs2) < 2:
            util.toc()
            continue
        
        # quick processing
        if config["quick"] and not treePartNeeded(config, part):
            util.toc()
            continue
        
        if len(part) > config["align_bigtree"]:
            tree = muscle.buildBigTree(seqs2)
            tree.writeNewick(bootTreeFile(config, i))
            util.toc()
            continue
        
        try:
            if config["align_fast"]:
                names, seqs3 = muscle.muscleFastOrdered(seqs2)
                aln = fasta.array2dict(names, seqs3)
            else:
                names, seqs3 = muscle.muscleOrdered(seqs2)
                aln = fasta.array2dict(names, seqs3)
        except:
            # sometimes muscle runs out of memory
            # if so, use clustalw instead
            aln = clustalw.clustalw(seqs2)
            names = aln.keys()

        fasta.writeFasta(file(alignFile(config, i), "w"), aln, order=names)
        util.toc()
    
    util.toc()



def runTreeBuilding(config):
    util.tic("SYNPHYL TREES")
    files = util.listFiles(config["align_dir"], config["align_suffix"])
    
    
    # estimation
    if "estimate" in config:
        print "will build %d trees" % len(files[config["tree_start"]:])
        util.toc()
        return
    

    for align in files[config["tree_start"]:config["tree_end"]]:
        util.log(align)
        aln = fasta.readFasta(align)
        
        if len(aln) > 2:
            if config["boot_iters"] > 1:
                trees = phylip.bootNeighbor(aln, config["bootiters"])
                phylip.writeBootTrees(align.replace(config["align_suffix"], 
                                      config["tree_boot_suffix"]), trees)
            else:
                tree = bionj.bionj(aln)
                tree.writeNewick(file(align.replace(config["align_suffix"], 
                                      config["tree_boot_suffix"]), "w"))
    
    util.toc()



def treePartNeeded(config, part):
    genomes = map(ensembl.id2genome, part)
    counts = util.histDict(genomes)
    return len(filter(util.gefunc(2), counts.values())) > 2 and \
           len(part) >= config["quick_minpartsize"]
    
    

def runOrthologySets(config):
    homology = phyloutil.Homology()
    homology.read(homologyFile(config))
    parts = phyloutil.homology2orthologSets(homology)
    
    util.writeDelim(orthologSetFile(config))
    
    


def runTreeProcessing(config):
    parts = util.readDelim(pulledPartFile(config))
    
    stree = algorithms.Tree()
    stree.readNewick(speciesTreeFile(config))
    homology = phyloutil.Homology()
    
    # estimation
    if "estimate" in config:
        n = 0    
        for i in xrange(len(parts)):
            part = parts[i]
            # decide if tree partitioning is need
            if treePartNeeded(config, part) and \
               os.path.isfile(bootTreeFile(config, i)):
                n += 1
        print "will process %d trees" % n
        return

    
    # process each partition
    for i in xrange(len(parts)):
        part = parts[i]
        treefile = bootTreeFile(config, i)
        
        # decide if tree partitioning is need
        if treePartNeeded(config, part) and \
           os.path.isfile(treefile):                    
            util.tic("process %s" % treefile)
            
            infile = file(treefile)
            
            # process bootstrapped trees
            for i in xrange(config["bootiters"]):
                util.tic("iter %d" % i)

                tree = algorithms.Tree()
                tree.readNewick(infile)
                tree = phyloutil.reconRoot(tree, stree, conf["gene2species"])
                trees = phyloutil.findBootTreeHomology(tree, stree, homology, 
                    config["bootiters"])
                
                if config["save_parttrees"]:
                    for i in range(len(trees)):
                        trees[i].writeNewick(treefile.replace(".tree", 
                                                              ".part%d.tree"))
                
                util.toc()
        else:
            util.tic("quick process %s" % treefile)
            phyloutil.findSimplePartHomology(parts[i], homology)
            util.toc()

    homology.write(file(homologyFile(config), "w"))
    oparts = phyloutil.homology2orthologSets(homology)
    
    util.writeDelim(orthologSetFile(config), oparts)



def countSubtreeTopologies(config):
    parts = util.readDelim(pulledPartFile(config))
    
    # count subtrees
    counts = {}
    for i in xrange(len(parts)):
        util.log("processing tree %d of %d" % (i, len(parts)))
        
        treefile = bootTreeFile(config, i)
        if os.path.isfile(treefile):
            tree = algorithms.Tree()
            tree.readNewick(treefile)
            
            trees = phyloutil.findRootedSubtrees(tree)
            
            for tree2 in trees:
                hashkey = algorithms.hashTree(tree2, ensembl.id2genome)
                if hashkey not in counts:
                    counts[hashkey] = 1
                else:
                    counts[hashkey] += 1
    
    return counts


def countTreeTopologies(config):
    parts = util.readDelim(pulledPartFile(config))
    stree = algorithms.Tree()
    stree.readNewick(speciesTreeFile(config))
    
    # count subtrees
    counts = {}
    for i in xrange(len(parts)):
        util.log("processing tree %d of %d" % (i, len(parts)))
        
        treefile = bootTreeFile(config, i)
        if os.path.isfile(treefile):
            tree = algorithms.Tree()
            tree.readNewick(treefile)
            tree = phyloutil.reconRoot(tree, stree)
            hashkey = algorithms.hashTree(tree, ensembl.id2genome)
            if hashkey not in counts:
                counts[hashkey] = 1
            else:
                counts[hashkey] += 1
    
    return counts


def runAltSpeciesSubtrees(config):
    counts = countSubtreeTopologies(config)
    
    # sort by counts
    keys = counts.keys()
    keys.sort(lambda a,b: cmp(counts[b], counts[a]))
    
    out = file(altSubtreesFile(config), "w")
    print >>out, "<html><body><table>"
    for key in keys[:500]:
        print >>out, "<tr><td>%s</td><td>%d</td></tr>" % (key, counts[key])
    print >>out, "</table></body></html>"
    out.close()


def runAltSpeciesTrees(config):
    counts = countTreeTopologies(config)
    
    # sort by counts
    keys = counts.keys()
    keys.sort(lambda a,b: cmp(counts[b], counts[a]))
    
    out = file(altTreesFile(config), "w")
    print >>out, "<html><body><table>"
    for key in keys[:500]:
        print >>out, "<tr><td>%s</td><td>%d</td></tr>" % (key, counts[key])
    print >>out, "</table></body></html>"
    out.close()
"""
    

def runVisSynteny(config):
    util.tic("SYNPHYL SYNTENY VISUALIZATION")
    
    util.tic("make svg")
    cmd = ("vissynteny.py "+
        " -S " + config["smap"] +
        " -g " + ",".join(config["genomes"]) + 
        " -c " + synPartFile(config) +
        " -s " + syntenyFile(config) +
        " -l -x " + str(config["vissynteny_context"]) +
        " -r " + config["refgenome"] + 
        (" -w %dx%d " % config["vissynteny_image_size"]) +
        " -v " + config["vissynteny_dir"] +
        " -W " + visSyntenyIndexFile(config))
    
    util.log("executing command: " + cmd)
    os.system(cmd)
    util.toc()

    
    # make png files
    util.tic("make png")
    os.system("for x in %s/*.svg; do echo $x; convert $x ${x/.svg/.png}; done"
              % config["vissynteny_dir"])
    util.toc()
    
    # make pdf files
    util.tic("make pdf")
    os.system("svg2pdf.py %s/*.svg" % config["vissynteny_dir"])
    util.toc()
    
    util.toc()
    
    

def runWebGeneration(config):
    util.tic("SYNPHYL WEB GENERATION")
    matching = genomeutil.Matching()
    genomeio.readGenomes(matching, config["genomes"])
        
    desc = readGeneDescriptions(config)
    
    if os.path.isfile(visSyntenyIndexFile(config)):
        config["vissynteny_index_file"] = visSyntenyIndexFile(config)
        util.log("using synteny visualization")
    
    if os.path.isfile(pulledPartFile(config)):
        util.tic("generate coarse report")
        parts = util.readDelim(pulledPartFile(config))
        config["report_dir"] = config["report_coarse_dir"]
        config["align_file_func"] = lambda partid: alignPublicFile(config, partid)
        config["tree_file_func"] = lambda partid: bootTreePublicFile(config, partid)
        synphylweb.writeAll(config, matching, parts, desc)
        util.toc()
    
    """
    if os.path.isfile(orthologSetFile(config)):
        util.tic("generate fine report")
        parts = util.readDelim(orthologSetFile(config))
        config["report_dir"] = config["report_fine_dir"]
        del config["align_file_func"]
        del config["tree_file_func"]
        synphylweb.writeAll(config, matching, parts, desc)
        util.toc()
    """
    
    util.toc()


def main():
    env.addEnvPaths("DATAPATH")

    # default config
    config = {}
    config.update(defaultConfig)
    
    # generate a configuration template
    if "genconf" in param:
        config["genconf"] = param["genconf"][-1]
        generateConf(config)
        return
    
    # read config file
    if (os.path.isfile(configfile)):
        execfile(configfile, config)
    setDefaultConfig(config, param)
    
    
    # setup logging
    util.globalTimer().addStream(file(logFile(config), "a"))
    util.log("=" * 78)
    util.log("SYNPHYL %s " % time.asctime())
    util.log("=" * 78)
    
    # show options
    util.log("Command line options are:")
    keys = param.keys()
    keys.sort()
    for key in keys:
        util.log("  %s = %s" % (key, str(param[key])))
    util.log()
    
    
    # parse parameters
    if "norefine" in param:
        config["synteny_refine"] = False
    
    if "nopull" in param:
        config["usepull"] = False
    
    if "quick" in param:
        config["quick"] = True
    
    if "estimate" in param:
        config["estimate"] = True
    
    if "showconfig" in param:
        showConfig(config)
            
    #if "paths" in param:        
    #    env.addPaths(* param["paths"])
    
    # stages
    if "synteny" in param:
        runSynteny(config)
    
    if "pull" in param:
        runGenePulling(config)
    
    if "align" in param:
        runAlignments(config)
    
    if "tree" in param:
        runTreeBuilding(config)
    
    if "altsubtrees" in param:
        runAltSpeciesSubtrees(config)
    
    if "alttrees" in param:
        runAltSpeciesTrees(config)
    
    if "homology" in param:
        runTreeProcessing(config)
    
    if "vissyn" in param:
        runVisSynteny(config)
    
    if "web" in param:
        runWebGeneration(config)
    
    if "report" in param:
        for report in param["report"]:
            if report == "synteny":
                os.system("partstats.py -p %s -g %s -S %s -c -s -m -z" % 
                    (synPartFile(config), ",".join(config["genomes"]),
                     config["smap"]))
            if report == "coarse":
                os.system("partstats.py -p %s -g %s -d %s -c -s -m -z" % 
                    (pulledPartFile(config), ",".join(config["genomes"]),
                     config["desc"]))
            if report == "fine":
                os.system("partstats.py -p %s -g %s -d %s -c -s -m -z" %
                    (orthologSetFile(config), ",".join(config["genomes"]),
                    config["desc"]))
                
    util.log("=" * 78)
    util.log("SYNPHYL COMPLETED %s " % time.asctime())
    util.log("=" * 78)

    
main()

