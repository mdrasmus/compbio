#
# SPIDIR (SPecies Informed DIstance-based Reconstruction) 
# Matt Rasmussen
# copyright 2007
#
# Makefile
#


CXX = g++
MEX = mex

CFLAGS = \
    -Wall \
    -I/usr/include/python2.4 \
    -I/util/include/python2.4


# matlab options
MATLAB_DIR = /afs/csail/i386_linux24/matlab/6.5r13sp1
MATLAB_CFLAGS = \
    -g \
    -I$(MATLAB_DIR)/extern/include/cpp \


#=============================================================================
# program files
SPIDIR_PROG = spidir

SPIDIR_SRC = \
    branchlen.cpp \
    common.cpp \
    likelihood.cpp \
    parsimony.cpp \
    search.cpp \
    Tree.cpp \

SPIDIR_OBJS = \
    branchlen.o \
    common.o \
    likelihood.o \
    parsimony.o \
    search.o \
    Tree.o \
    Sequences.o 

PROG_SRC = spidir.cpp 
PROG_OBJS = spidir.o $(SPIDIR_OBJS)
PROG_LIBS = 

LIBSPIDIR = lib/libspidir.a
LIBSPIDIR_OBJS = $(SPIDIR_OBJS)

# python files
PYTHON_MODULE = pyspidir.so
PYTHON_MODULE_OBJS = \
    pyspidir.o \
    $(SPIDIR_OBJS)

# matlab files
MATLAB_FUNC = matlab_spidir
MATLAB_SRC = $(SPIDIR_SRC) matlab_spidir.cpp
MATLAB_OBJ = matlab_spidir.mexglx


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

# C-library
$(LIBSPIDIR): $(LIBSPIDIR_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBSPIDIR) $(LIBSPIDIR_OBJS)

# testing program
test_spidir: $(SPIDIR_OBJS) test.o
	$(CXX) $(SPIDIR_OBJS) $(CFLAGS) test.o $(PROG_LIBS) -o test_spidir

# python module
$(PYTHON_MODULE): $(PYTHON_MODULE_OBJS)
	$(CXX) -shared $(PYTHON_MODULE_OBJS) -o $(PYTHON_MODULE)

# matlab function
$(MATLAB_FUNC): $(MATLAB_SRC)
	$(MEX) -v $(MATLAB_CFLAGS) $(MATLAB_SRC) -o $(MATLAB_FUNC)



#=============================================================================
# basic rules

$(PYTHON_MODULE_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


install: $(SPIDIR_PROG) $(PYTHON_MODULE) test_spidir
	cp $(SPIDIR_PROG) test_spidir ../bin
	cp $(PYTHON_MODULE) ../python


clean:
	rm -rf $(PROG_OBJS) $(SPIDIR_PROG) $(LIBSPIDIR) \
                $(PYTHON_MODULE_OBJS) $(PYTHON_MODULE) \
                $(MATLAB_OBJ) \
	        test.o test_spidir
