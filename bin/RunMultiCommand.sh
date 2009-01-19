#!/bin/sh
# adopted from pouya (Tue Sep  2 12:59:38 EDT 2008)

# Takes a file, breaks it into chunks, sends each chunk to a different
# queue with command.  If file is not supplied, simply passes on an
# integer to each of the runs (it is assumed that the command
# internally parallelizes)

CHUNKS=10
SEP=" "
QUEUE=broad
OUTCMD="cat"
COMBCMD="cat"   
FASTAINPUT=0
ILINPUT=0

if [ $# -eq 0 ]; then
		echo >&2 "USAGE: $0 [OPTIONS] 
 -c     Command to use [required]
 -i     Input file, if empty then no splitting and index is passed into command
 -o     Output file [required]
 -d     Temporary directory to use [required]
 -n     Number of jobs to spawn [default: $CHUNKS]
 -s     Separator to use [default: (space)]
 -q     Queue to use [default: $QUEUE]
 -f     FASTA input to be broken up by sequence [default: $FASTAINPUT]
 -I     Allow input to subprocesses to be interleaved (rather than kept in 
        contiguous blocks) [default: $ILINPUT]
 -t     Command to run on the final output [default: $OUTCMD]
 -T     Command used to combine the final output files [default: $COMBCMD]
"
	exit 1
fi

while getopts c:i:o:d:n:s:q:t:T:f:I: o
do      case "$o" in
		c)		CMD="$OPTARG";;
		i)		INFILE="$OPTARG";;
		o)		OUTFILE="$OPTARG";;
		d)		DIR="$OPTARG";;
		n)		CHUNKS="$OPTARG";;
		s)		SEP="$OPTARG";;
		q)		QUEUE="$OPTARG";;
		t)		OUTCMD="$OPTARG";;
		T)		COMBCMD="$OPTARG";;
		f)		FASTAINPUT="$OPTARG";;
		I)		ILINPUT="$OPTARG";;
        [?])    echo >&2 "ERROR: command line parameter not recognized."; 
	exit 1;;
        esac
done

if [ -z "$CMD" -o -z "$OUTFILE" -o  -z "$DIR"  ]; then
        echo >&2 "ERROR: -c -o -d are required."; exit 1
fi

mkdir -p $DIR

# if supplied, split the input file
if [ ! -z "$INFILE" ]; then
	# use these to build the awk commands on the basis of FASTAINPUT and 
        # ILINPUT
	# if FASTAINPUT then make X indicate the current FASTA sequence number 
        # and use when splitting rather than the line number
	if [ "$FASTAINPUT" = 0 ]; then
		AWKSPLITPRE=""
		AWKSPLITVAR="NR"
	else
		AWKSPLITPRE="/^>/{X++};"
		AWKSPLITVAR="X"
	fi

	# if ILINPUT = 0, then split into nearly equal sized chunks 
	# (1,1,1,2,2,2,3,3,3), otherwise interleave the output 
	# (1,2,3,1,2,3,...)
	if [ "$ILINPUT" = 0 ]; then
		AWKSPLITNUM="(int(($AWKSPLITVAR-1)*NC/NL)+1)"
	else
		AWKSPLITNUM="((($AWKSPLITVAR-1)%NC)+1)"
	fi

	# count the number of "lines"
	LINES=$(awk ''$AWKSPLITPRE' END{print '$AWKSPLITVAR'+0}' $INFILE)

	# the number of chunks cannot be greater than the number of "lines"
	if [ $LINES -lt $CHUNKS ]; then
		CHUNKS=$LINES
	fi

	# break file into chunks
	awk -vNC=$CHUNKS -vNL=$LINES \
	''$AWKSPLITPRE' {print > ("'$DIR'/"'$AWKSPLITNUM'".in")}' $INFILE
fi

for i in `seq 1 $CHUNKS`
do
	if [ ! -z "$INFILE" ]; then
		CMDINPUT=$DIR/$i.in
	else
		CMDINPUT=$i
	fi

	bsub -o $DIR/$i.bout -e $DIR/$i.berr -q $QUEUE -K 
		"`printf "$CMD$SEP%s" $CMDINPUT` > $DIR/$i.out" &
	OUTFILES="$OUTFILES $DIR/$i.out"
	PIDS="$PIDS $!"
done

wait $PIDS

eval $COMBCMD $OUTFILES | eval $OUTCMD > $OUTFILE

for i in `seq 1 $CHUNKS`
do
	rm -rf $DIR/$i.in $DIR/$i.out
done

