#!/bin/sh

OPTIONS=$*
QUIET=false


if [ "$1" = "-v" ]; then
    echo DATAPATH=
    echo -e ${DATAPATH//:/\\n}
    exit 0
fi


# reset datapath
export DATAPATH=.


# fungi
fungi() {
    export DATAPATH=$DATAPATH:$FILEPATH/data/fungi
    export DATAPATH=$DATAPATH:$FILEPATH/data/fungi/gene_coord
    export DATAPATH=$DATAPATH:$FILEPATH/data/fungi/pep
    export DATAPATH=$DATAPATH:$FILEPATH/data/fungi/nt
    export DATAPATH=$DATAPATH:$FILEPATH/data/fungi/blast_pep
}

fly_flybase() {
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/gene_coord
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/pep
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/nt
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/blast_pep
}

fly_flybase_single() {
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/gene_coord_single
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/pep
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/nt
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/flybase/blast_pep
}

fly_genscan() {
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/genscan
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/genscan/gene_coord
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/genscan/pep
    export DATAPATH=$DATAPATH:$FILEPATH/data/flies/genscan/blast_pep
}


# always include fungi
fungi


# process options
for OPTION in $OPTIONS; do
    if [ $OPTION = "-q" ]; then
        QUIET=true
    else
        $OPTION
    fi
done



# print out variable settings
if [ $QUIET = false ]; then
    echo DATAPATH=
    echo -e ${DATAPATH//:/\\n}
fi
