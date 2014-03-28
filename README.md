compbio
=======

Python libraries and utilities for computational biology.

## About

This package contains algorithms related to several areas of genomics,
phylogenetics, and population genetics. Some of the highlights
include:

- reading, writing, and manipulating phylogenetic trees
- reconciling gene-trees with species-trees
- inferring gene duplications, losses, and horizontal transfers
- methods for coalescent processes, incomplete lineage sorting
- methods for ancestral recombination graphs (ARGs)
- finding syntenic regions (i.e. co-linear orthology) between genomes
- processing common file formats: FASTA, PHYLIP, newick, nexus, etc.

In addition to computational biology-specific methods, this package also
contains general utilities for working with scientific data:

- sparse matrix file formats
- reading, writing, and manipulating tables of data
- working with intervals (e.g. intersection, union, etc)
- plotting (Gnuplot, Rpy)
- statistics
- general data-structures and algorithms:
  quad trees, Union-Find, HHMs, clustering

## Download

The `compbio` package is available for download from several sources:

- https://github.com/mdrasmus/compbio
- https://pypi.python.org/pypi/compbio

## Requirements

Most modules in this package can be used without any additional dependencies.

For plotting modules, the dependencies include:

- [gnuplot](http://www.gnuplot.info/)
- [R](http://www.r-project.org/)

For some scientific methods, the dependencies include:

- [scipy](http://www.scipy.org/)

For development of the `compbio` package itself, dependencies can be installed
with `pip`:

```
pip install -r requirements-dev.txt
```

## INSTALL

The `compbio` package is available on pypi, and can be installed using pip:
```
pip install compbio
```

These packages can be installed from the source directory using:
```
python setup.py install
```

Optionally, the libraries can be used directly from the source directory by
configuring one's environment variables as follows (assuming bash shell):

```
export PATH=$PATH:path/to/compbio/bin
export PYTHONPATH=$PYTHONPATH:path/to/compbio
```

## Author

These libraries were built up over the course of the Ph.D. of the
author, Matthew D. Rasmussen (http://mattrasmus.com,
<rasmus@alum.mit.edu>). Many of the methods here were utilized in
several published software projects including:

- [ARGweaver](http://mdrasmus.github.io/argweaver/):
  Rasmussen, Siepel. Genome-wide inference of ancestral
  recombination graphs. ArXiv. 2013.
- [DLCoal](http://compbio.mit.edu/dlcoal/):
  Rasmussen, Kellis. Unified modeling of gene duplication,
  loss, and coalescence using a locus tree.  Genome Research. 2012.
- [SPIMAP](http://compbio.mit.edu/spimap/):
  Rasmussen, Kellis. A Bayesian approach for fast and accurate
  gene tree reconstruction. Molecular Biology and Evolution. 2010.
- [SPIDIR](http://compbio.mit.edu/spidir/):
  Rasmussen, Kellis. Accurate gene-tree reconstruction by
  learning gene- and species-specific substitution rates across multiple
  complete genomes. Genome Research. 2007.

*Minor note:* Although the libraries of this package supports each of
these software packages, this package is not a required
dependency. Instead each software package contains its own private
copy of modules taken from this package.
