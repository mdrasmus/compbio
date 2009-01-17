#!/bin/sh

ITER=1000
SAMPLES=1000

../../spidir \
    -p ../data/flies.nt.param \
    -S ../data/flies.smap \
    -s ../data/flies.stree \
    -a ../data/0.nt.align \
    -o ../output/fly_mcmc3 \
    --search mcmc \
    --lengths hky \
    --lkfunc none \
    -V 2 \
    -i $ITER \
    --correct ../data/0.nt.tree $*

# get posterior probabilities
echo "calculating posteriors..."
tail -n+1 ../output/fly_mcmc3.mcmc | thin $(( $ITER / $SAMPLES )) | \
    viewtree.py -S ../data/flies.smap -s ../data/flies.stree -i -  \
    > ../output/fly_mcmc3.counts
