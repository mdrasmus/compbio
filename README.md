compbio
-------

A collection of python libraries and utilities for computational biology.

The package contains implementations of algorithms related to several
areas of genomics, phylogenetics, and population genetics. Some of the
highlights include:

- reading, writing, and manipulating phylogenetic trees
- reconciling gene-trees with species-trees
- inferring gene duplications, losses, and horizontal transfers
- methods for coalescent processes, incomplete lineage sorting
- methods for ancestral recombination graphs (ARGs)
- finding syntenic regions (i.e. co-linear orthology) between genomes
- processing common file formats: FASTA, PHYLIP, newick, nexus, etc.

In addition to computational biology-specific methods this package also
contains general utilities for working with scientific data:

- sparse matrix file formats
- reading, writing, and manipulating tables of data
- working with intervals (e.g. intersection, union, etc)
- plotting (Gnuplot, Rpy)
- statistics
- general data-structures and algorithms:
  quad trees, Union-Find, HHMs, clustering

These libraries were built up over the course of the author's Ph.D.
(Matthew D. Rasmussen <rasmus@alum.mit.edu>). Many of the methods here were
utilized in several published projects including: SPIDIR, SPIMAP, DLCoal, and
ARGweaver.

- ARGweaver: Rasmussen, Siepel. Genome-wide inference of ancestral
recombination graphs. ArXiv. 2013/
- DLCoal: Rasmussen, Kellis. Unified modeling of gene duplication,
loss, and coalescence using a locus tree.  Genome Research. 2012.
- SPIMAP: Rasmussen, Kellis. A Bayesian approach for fast and accurate
gene tree reconstruction. Molecular Biology and Evolution. 2010.
- SPIDIR: Rasmussen, Kellis. Accurate gene-tree reconstruction by
learning gene- and species-specific substitution rates across multiple
complete genomes. Genome Research. 2007.
