---------------------------------------------------------------------------------
Fri Jun 30 12:35:25 EDT 2006
Matt Rasmussen

To use my 'rasmus' scripts, you will need to set the following environment
variables (add them to your ~/.bashrc file if you use bash):

# for my binaries
export PATH=$PATH:<checked-out-dir>/rasmus/bin

# for python libraries
export PYTHONPATH=$PYTHONPATH:<checked-out-dir>/rasmus/python

# for summon visualization
export PYTHONPATH=$PYTHONPATH:/afs/csail.mit.edu/group/compbio/software/summon/lib
export PATH=$PATH:/afs/csail.mit.edu/group/compbio/bin


<checked-out-dir> refers to the location where you checked out the CVS files.


To checkout my CVS directory use the following CVS command:


export CVS_RSH=ssh
export CVSROOT=compbio.mit.edu:/afs/csail.mit.edu/group/compbio/cvsroot
cvs co rasmus



