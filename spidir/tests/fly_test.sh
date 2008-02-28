#!/bin/sh


../spidir \
    -p ../test/flies.nt.param \
    -S ../test/flies.smap \
    -s ../test/flies.stree \
    -a ../test/0.nt.align \
    -o fly \
    -V 2 \
    --correct ../test/0.nt.tree $*
