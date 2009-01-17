#!/bin/sh

ITER=1000
SAMPLES=1000

../../spidir \
    -p ../data/flies.nt.param \
    -S ../data/flies.smap \
    -s ../data/flies.stree \
    -a ../data/0.nt.align \
    -o ../output/fly_mcmc \
    --search mcmc \
    --lengths spidir \
    --lkfunc hky \
    -V 2 \
    -i $ITER \
    --correct ../data/0.nt.tree $*

# get posterior probabilities
echo "calculating posteriors..."
cat ../output/fly_mcmc.mcmc | thin $(( $ITER / $SAMPLES )) |
    viewtree.py -i -  \
    > ../output/fly_mcmc.counts
