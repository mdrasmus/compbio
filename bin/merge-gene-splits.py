#!/usr/bin/env python

# python libs
import os
import sys


# rasmus libs
from rasmus import env
from rasmus import fasta
from rasmus import genecall
from rasmus import genomeutil
from rasmus import gff
from rasmus import regionlib
from rasmus import util


options = [
    ["S:", "smap=", "smap", "<gene2species mapping>",
        {"single": True}],
    ["r:", "gff=", "gff", "<gff file>"],
    ["p:", "pep=", "pep", "<peptide fasta file>"],
    ["n:", "nt=", "nt", "<nucleotide fasta file>"],
    ["P:", "outpep=", "outpep", "<output peptide fasta file>",
        {"single": True}],
    ["N:", "outnt=", "outnt", "<output nucleotide fasta file>",
        {"single": True}],
    ["A:", "alignext=", "alignext", "<alignment file extension>",
        {"single": True}],
    ["B:", "alignoldext=", "alignoldext", "<old alignment file extension>",
        {"single": True}],        
    ["F:", "fastaext=", "fastaext", "<new fasta file extension>",
        {"single": True}],
    ["G:", "fastaoldext=", "fastaoldext", "<old fasta file extension>",
        {"single": True}]
]


conf = util.parseOptions(sys.argv, options, resthelp="<alignments> ...")

def main(conf):
    env.addEnvPaths("DATAPATH")
    
    # read gene2species
    gene2species = genomeutil.readGene2species(env.findFile(conf["smap"]))

    # read gene regions
    util.tic("read gene regions")
    regions = []
    for gffFile in conf["gff"]:
        util.tic("read " + gffFile)
        lst = gff.readGff(env.findFile(gffFile), format=gff.GFF3,
                          regionFilter = lambda x: x.feature == "gene")
        for gene in lst:
            gene.species = gene2species(gene.data["ID"])
        regions.extend(lst)
        util.toc()
    regiondb = regionlib.RegionDb(regions)
    util.toc()
    
    
    # read sequence data
    util.tic("read sequence data")
    pepseqs = fasta.FastaDict()
    ntseqs = fasta.FastaDict()
    
    for pepfile in conf["pep"]:
        pepseqs.read(env.findFile(pepfile), useIndex=True)
    
    for ntfile in conf["nt"]:
        ntseqs.read(env.findFile(ntfile), useIndex=True)
    util.toc()
    
    
    # look for split genes in alignments
    outpepseqs = fasta.FastaDict()
    outntseqs = fasta.FastaDict()
    
    util.tic("process alignments")
    for alignfile in conf["REST"]:
        aln = fasta.readFasta(alignfile)
        
        merges = genecall.findFragments(regiondb, aln, overlapCutoff=.5)
        
        if len(merges) > 0:
            util.tic(alignfile)
            util.logger("merges", merges)
            
            oldnames = util.remove(aln.keys(), * util.flatten(merges))
            newfasta = pepseqs.get(oldnames)
            
            for merge in merges:
                newname = "+".join(merge)            
                peps = []
                nts = []
                
                for name in merge:
                    if name not in pepseqs:
                        util.logger(name, "is an unknown peptide sequence!")
                    else:
                        peps.append(pepseqs[name])
                    
                    if name not in ntseqs:
                        util.logger(name, "is an unknown nucleotide sequence!")
                    else:
                        nts.append(ntseqs[name])
                
                
                if len(peps) == len(merge):
                    outpepseqs.add(newname, "".join(peps))
                    newfasta.add(newname, "".join(peps))
                    util.logger("created '%s' peptide sequence" % newname)
                else:
                    util.logger("unable to create necessary peptide sequence '%s'!" % newname)
                    break
                
                if len(nts) == len(merge):
                    outntseqs.add(newname, "".join(nts))
                    util.logger("created '%s' nucleotide sequence" % newname)
                else:
                    util.logger("unable to create nucleotide sequence '%s'!" % newname)
                                    
            else:
                # create new fasta with new gene names
                fafile = util.replace_ext(alignfile, conf["alignext"], 
                                         conf["fastaext"])
                faoldfile = util.replace_ext(alignfile, conf["alignext"], 
                                         conf["fastaoldext"])
                alignoldfile = util.replace_ext(alignfile, conf["alignext"], 
                                         conf["alignoldext"])

                os.rename(fafile, faoldfile)
                newfasta.write(fafile)
                os.rename(alignfile, alignoldfile)
                util.logger("wrote '%s'" % fafile)
                util.logger("backup '%s'" % faoldfile)
                util.logger("backup '%s'" % alignoldfile)
            util.toc()            


    # write out new sequence files
    outpepseqs.write(conf["outpep"])
    outntseqs.write(conf["outnt"])
    util.toc()

main(conf)
