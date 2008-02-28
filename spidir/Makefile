#
# SPIDIR (SPecies Informed DIstance-based Reconstruction) 
# Matt Rasmussen
# copyright 2007
#
# Makefile
#


# install prefix paths
prefix = /usr
python_prefix = /usr/lib/python2.4/site-packages/


CXX = g++
MEX = mex

CFLAGS = \
    -Wall -fPIC \
    -I/usr/include/python2.4 \
    -I/util/include/python2.4


# matlab options
MATLAB_DIR = /afs/csail/i386_linux24/matlab/2007a
#MATLAB_DIR = /afs/csail/i386_linux24/matlab/6.5r13sp1
MATLAB_CFLAGS = \
    -g \
    -I$(MATLAB_DIR)/extern/include/cpp 
MEX_EXT = mexglx

#=============================================================================
# target names

# program files
SPIDIR_PROG = spidir

SPIDIR_SRC = \
    spidir.cpp \
    mldist.cpp \
    common.cpp \
    likelihood.cpp \
    parsimony.cpp \
    phylogeny.cpp \
    search.cpp \
    Tree.cpp \

SPIDIR_OBJS = \
    spidir.o \
    mldist.o \
    common.o \
    likelihood.o \
    parsimony.o \
    phylogeny.o \
    search.o \
    Tree.o \
    Sequences.o 

PROG_SRC = spidir_main.cpp 
PROG_OBJS = spidir_main.o $(SPIDIR_OBJS)
PROG_LIBS = 
#`gsl-config --libs`


# C-library files
LIBSPIDIR = lib/libspidir.a
LIBSPIDIR_SHARED = lib/libspidir.so
LIBSPIDIR_OBJS = $(SPIDIR_OBJS)


# python files
PYTHON_MODULE = pyspidir.so
PYTHON_MODULE_OBJS = \
    pyspidir.o \
    $(SPIDIR_OBJS)  
PYTHON_MODULE_LIBS = 
#-lpython2.4 
#`gsl-config --libs`


# matlab files              
MATLAB_OBJS = matlab/spidir_treelk.$(MEX_EXT) \
              matlab/spidir_display_tree.$(MEX_EXT) \
              matlab/spidir_genbranches.$(MEX_EXT) \
              matlab/spidir_mlhkydist.$(MEX_EXT) \
              matlab/spidir_neighborjoin.$(MEX_EXT) \
              matlab/spidir_reconcile.$(MEX_EXT) \
              matlab/spidir_readtree.$(MEX_EXT)

MATLAB_COMPILE = spidir_matlab_compile.m
MATLAB_COMPILE_RULES = \
              matlab/spidir_treelk.rule \
              matlab/spidir_display_tree.rule \
              matlab/spidir_genbranches.rule \
              matlab/spidir_mlhkydist.rule \
              matlab/spidir_neighborjoin.rule \
              matlab/spidir_reconcile.rule \
              matlab/spidir_readtree.rule

MATLAB_SRC = $(SPIDIR_SRC) matlab_interface.cpp


#=============================================================================
# optional CFLAGS

# profiling
ifdef PROFILE
	CFLAGS := $(CFLAGS) -pg
endif

# debugging
ifdef DEBUG
	CFLAGS := $(CFLAGS) -g
else
	CFLAGS := $(CFLAGS) -O3
endif


#=============================================================================
# targets

all: $(SPIDIR_PROG) $(LIBSPIDIR) $(PYTHON_MODULE) test_spidir

# stand-alone program
$(SPIDIR_PROG): $(PROG_OBJS)
	$(CXX) $(CFLAGS) $(PROG_OBJS) $(PROG_LIBS) -o $(SPIDIR_PROG)

maxml: maxml.o $(SPIDIR_OBJS)
	$(CXX) $(CFLAGS) maxml.o $(SPIDIR_OBJS) $(PROG_LIBS) -o maxml

# C-library
spidirlib: $(LIBSPIDIR) $(LIBSPIDIR_SHARED)

$(LIBSPIDIR): $(LIBSPIDIR_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBSPIDIR) $(LIBSPIDIR_OBJS)

$(LIBSPIDIR_SHARED): $(LIBSPIDIR_OBJS)
	mkdir -p lib
	$(CXX) -o $(LIBSPIDIR_SHARED) -shared $(LIBSPIDIR_OBJS)

# testing program
test_spidir: $(SPIDIR_OBJS) test.o
	$(CXX) $(SPIDIR_OBJS) $(CFLAGS) test.o $(PROG_LIBS) -o test_spidir

# python module
#pyspidir: $(PYTHON_MODULE)

$(PYTHON_MODULE): $(PYTHON_MODULE_OBJS)
	$(CXX) -shared $(PYTHON_MODULE_OBJS) $(PYTHON_MODULE_LIBS) -o $(PYTHON_MODULE)


# matlab interface
matlab: $(MATLAB_OBJS) $(MATLAB_COMPILE)


$(MATLAB_OBJS): %.$(MEX_EXT): %.cpp
	$(MEX) $(MATLAB_CFLAGS) $(MATLAB_SRC) $< -o $@

# generate compile rules for windows
$(MATLAB_COMPILE): $(MATLAB_COMPILE_RULES)
$(MATLAB_COMPILE_RULES): %.rule: %.cpp
	echo "display('compiling $<...');" >> $(MATLAB_COMPILE)
	echo $(MEX) $(MATLAB_SRC) $< -o $(@:%.rule=%) >> $(MATLAB_COMPILE)
	touch $@

#=============================================================================
# basic rules

$(PYTHON_MODULE_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


install: $(SPIDIR_PROG)
	cp $(SPIDIR_PROG) $(prefix)/bin
	
installpy:
	cp $(PYTHON_MODULE) $(prefix_python)
	

myinstall: $(SPIDIR_PROG) $(PYTHON_MODULE) test_spidir maxml
	cp $(SPIDIR_PROG) test_spidir maxml ../bin
	cp $(PYTHON_MODULE) ../python

myinstall64: $(SPIDIR_PROG) $(PYTHON_MODULE) test_spidir maxml
	cp $(SPIDIR_PROG) test_spidir maxml ../bin64
	cp $(PYTHON_MODULE) ../python64

clean:
	rm -f $(PROG_OBJS) $(SPIDIR_PROG) $(LIBSPIDIR) \
              $(PYTHON_MODULE_OBJS) $(PYTHON_MODULE) \
              $(MATLAB_OBJS) maxml maxml.o \
              $(MATLAB_COMPILE) $(MATLAB_COMPILE_RULES) \
	      test.o test_spidir
