#!/bin/sh


../../spidir \
    -p ../data/flies.nt.param \
    -S ../data/flies.smap \
    -s ../data/flies.stree \
    -a ../data/0.nt.align \
    -o ../output/fly \
    --search climb \
    -V 2 \
    --correct ../data/0.nt.tree $*
