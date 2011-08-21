#!/bin/bash
# Example script of analyzing fungi gene families
# Fri Aug 19 15:30:29 EDT 2011


# External dependencies
#
# BLAST programs:
#   formatdb
#   blastall
#
#   you can install these on Debian with 
#     apt-get install blast2
#
# MUSCLE sequence alignment program:
#   muscle
#
#   you can install this on Debian with 
#     apt-get install muscle
#
# RAXML phylogenetic reconstruction
#   raxmlHPC
#

#
# In practice, one would rewrite these commands in order to run programs
# in parallel on a compute cluster.
#


# untar fungi example dataset
tar zxvf data/fungi.tar.gz -C data



# how many sequences are in this example?
grep -c '>' data/fungi/fasta-pep/*.fa
:<<EOF
data/fungi/fasta-pep/agos.pep.fa:59
data/fungi/fasta-pep/calb.pep.fa:73
data/fungi/fasta-pep/cgla.pep.fa:73
data/fungi/fasta-pep/cgui.pep.fa:71
data/fungi/fasta-pep/clus.pep.fa:74
data/fungi/fasta-pep/cpar.pep.fa:72
data/fungi/fasta-pep/ctro.pep.fa:78
data/fungi/fasta-pep/dhan.pep.fa:71
data/fungi/fasta-pep/klac.pep.fa:66
data/fungi/fasta-pep/kwal.pep.fa:73
data/fungi/fasta-pep/lelo.pep.fa:72
data/fungi/fasta-pep/sbay.pep.fa:97
data/fungi/fasta-pep/scas.pep.fa:66
data/fungi/fasta-pep/scer.pep.fa:81
data/fungi/fasta-pep/smik.pep.fa:97
data/fungi/fasta-pep/spar.pep.fa:95
EOF


# perform pairwise blast
# keep only the best hit between any pair of proteins
mkdir -p output/fungi/blast
for pep1 in data/fungi/fasta-pep/*.fa; do
    formatdb -i $pep1 -o T

    for pep2 in data/fungi/fasta-pep/*.fa; do
	if [[ $pep1 < $pep2 ]]; then
	    sp1=$(basename ${pep1/.pep.fa})
	    sp2=$(basename ${pep2/.pep.fa})

	    echo blasting $sp1 $sp2...
	    blastall -p blastp -d $pep1 -i $pep2 -m8 | \
		blastbest > output/fungi/blast/$sp1-$sp2.blast
	fi
    done
done


# how many blast hits do we have?
wc -l output/fungi/blast/*

# filter hits for >40 bits and find connected components
awk -vOFS="\t" '$12 > 40 {print $1, $2}' output/fungi/blast/*.blast \
    | connected-comp > output/fungi/fungi.part.txt

# what is family size distribution?
awk '{print NF}' output/fungi/fungi.part.txt | sort -g | uniq -c


# partition fasta sequences into family directories
partfasta -p output/fungi/fungi.part.txt \
    -o output/fungi/fams \
    -F .pep.fa \
    data/fungi/fasta-pep/*.fa

# list all family fasta files
phylofiles -e output/fungi/fams .pep.fa


# make peptide alignments with MUSCLE
phylofiles -e output/fungi/fams .pep.fa | (while read fa; do
    muscle < $fa > ${fa/.fa/.align}
done)

# view an alignment
viewalign output/fungi/fams/0/0.pep.align

# ensure all alignments completed
phylofiles -e output/fungi/fams .pep.align | wc -l


# reverse translate peptide alignments into dna alignments
phylofiles -e output/fungi/fams .pep.align | xargs revtrans \
    $(for x in data/fungi/fasta-nt/*.fa;do echo -d $x; done) \
    -o .pep.align \
    -n .nt.align

# ensure all nucleotide alignments are complete
phylofiles -e output/fungi/fams .nt.align | wc -l


# make raxml trees
phylofiles -n output/fungi/fams .nt.raxml.tree | xargs run-raxml \
    -A .nt.align \
    -T .nt.raxml.tree \
    -O .nt.raxml.output \
    --seqtype=dna

# view a raxml tree
viewtree output/fungi/fams/0/0.nt.raxml.tree


# ensure all trees are complete
phylofiles -e output/fungi/fams .nt.raxml.tree | wc -l


# reroot gene trees
phylofiles -e output/fungi/fams .nt.raxml.tree | xargs reconroot \
    -s data/config/fungi.stree \
    -S data/config/fungi.smap 
    

# evaluate the gene trees using the duplication consistency score
phylofiles -e output/fungi/fams .nt.raxml.tree | dupconsistency \
    -s data/config/fungi.stree \
    -S data/config/fungi.smap \
    > output/fungi/fungi.nt.raxml.dupcons.txt

cut -f4 output/fungi/fungi.nt.raxml.dupcons.txt | hist
:<<EOF
item            count  percent  
0.0               210  79.8479  
1.0                16   6.0837  
0.333333333333      8   3.0418  
0.25                7   2.6616  
0.5                 6   2.2814  
0.111111111111      4   1.5209  
0.666666666667      3   1.1407  
0.166666666667      2   0.7605  
0.142857142857      2   0.7605  
0.75                2   0.7605  
0.125               1   0.3802  
0.0625              1   0.3802  
0.444444444444      1   0.3802  
EOF


#=============================================================================
# clean up

# remove BLAST index files
rm -f data/fungi/fasta-pep/*.{phr,pin,psd,psq,psi}
rm -f formatdb.log

# remove output for fungi
rm -rf output/fungi


#=============================================================================
# MISC code

# make initial files

phylofiles data/fungi-fams/ .nt.align | xargs cat | python -c '
from rasmus.common import *
from compbio import fasta

gene2species = phylo.read_gene2species("data/config/fungi.smap")

db = defaultdict(lambda: fasta.FastaDict())

for key, seq in fasta.iter_fasta(sys.stdin):
  seq = seq.replace("-", "")

  seq = seq.replace("CYT", "CCT")
  seq = seq.replace("GTW", "GTA")

  db[gene2species(key)][key] = seq

for sp in db:
  db[sp].write("data/fungi/fasta-nt/%s.nt.fa" % sp)

  pep = fasta.FastaDict()

  for key in db[sp]:
    pep[key] = seqlib.translate(db[sp][key])
  
  pep.write("data/fungi/fasta-pep/%s.pep.fa" % sp)
'

