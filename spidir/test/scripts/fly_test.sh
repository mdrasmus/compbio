#!/bin/sh


../../spidir \
    -p ../data/flies.nt.param \
    -S ../data/flies.smap \
    -s ../data/flies.stree \
    -a ../data/0.nt.align \
    -o ../output/fly \
    --search climb \
    -D 0.2 \
    -L 0.19 \
    -V 2 \
    -i 1000 \
    --correct ../data/0.nt.tree $*
