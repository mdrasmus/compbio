#!/bin/sh

ITER=100000
SAMPLES=1000
FILE=1

../../spidir \
    -p ~/work/spidir/train/flies.subset.spidir_ml.nt.param \
    -S ../data/flies.smap \
    -s ~/work/spidir/train/flies.subset.stree \
    -a ~/work/spidir/real/flies/ones-subset/data/$FILE/$FILE.nt.align \
    -o ../output/fly_mcmc4 \
    --search mcmc \
    --lengths spidir \
    --lkfunc hky \
    -V 2 \
    -i $ITER \
    --correct ~/work/spidir/real/flies/ones-subset/data/$FILE/$FILE.correct.tree $*

# get posterior probabilities
echo "calculating posteriors..."
# tail -n+1000
time cat ../output/fly_mcmc4.mcmc | thin $(( $ITER / $SAMPLES )) |
    viewtree.py -i -  \
    > ../output/fly_mcmc4.counts
