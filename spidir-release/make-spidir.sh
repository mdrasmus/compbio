#!/bin/sh

SPIDIR=spidir-0.7

mkdir -p ${SPIDIR}
mkdir -p ${SPIDIR}/lib
mkdir -p ${SPIDIR}/lib/Spidir
mkdir -p ${SPIDIR}/lib/rasmus
mkdir -p ${SPIDIR}/lib/rasmus/vis
mkdir -p ${SPIDIR}/lib/rasmus/bio
mkdir -p ${SPIDIR}/bin


LIBS="__init__ \
    algorithms \
    cluster \
    env \
    graph \
    matrix \
    options \
    plotting \
    progress \
    regionlib \
    stats \
    svg \
    tablelib \
    textdraw \
    timer \
    treelib \
    util \
    vector \
    bio/__init__ \
    bio/bionj \
    bio/blast \
    bio/fasta \
    bio/genomeutil \
    bio/gff \
    bio/phylip \
    bio/phylo \
    bio/seqlib \
    vis/__init__ \
    vis/treesvg \
    "


cp -r ../python/Spidir/*.py ${SPIDIR}/lib/Spidir

if [ ! -d ${SPIDIR}/spidir-examples ]; then
    cp -r spidir-examples ${SPIDIR}
fi

for LIB in $LIBS; do
    cp ../python/rasmus/${LIB}.py ${SPIDIR}/lib/rasmus/${LIB}.py
done


cp ../bin/spidir.py ${SPIDIR}/bin

tar zcvf ${SPIDIR}.tar.gz ${SPIDIR}
