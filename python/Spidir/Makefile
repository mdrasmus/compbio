PYTHON_MODULE = pyspidir.so
PROG_OBJS = \
    branchlen.o \
    common.o \
    likelihood.o \
    parsimony.o \
    search.o \
    Tree.o \

PYTHON_MODULE_OBJS = \
    pyspidir.o \
    $(PROG_OBJS)


CC = g++

CFLAGS = \
    -O3 \
    -Wall \
    -I/usr/include/python2.4 \
    -I/util/include/python2.4

all: $(PYTHON_MODULE) spidir test_spidir

# stand along program
spidir: $(PROG_OBJS) spidir.o
	g++ $(PROG_OBJS) $(CFLAGS) spidir.o -o spidir

# testing program
test_spidir: $(PROG_OBJS) test.o
	g++ $(PROG_OBJS) $(CFLAGS) test.o -o test_spidir

# python module
$(PYTHON_MODULE): $(PYTHON_MODULE_OBJS)
	g++ -shared $(PYTHON_MODULE_OBJS) -o $(PYTHON_MODULE)




# basic compile rule
$(PYTHON_MODULE_OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<


install: spidir
	cp spidir ../../bin/

clean:
	rm -rf $(PYTHON_MODULE_OBJS) $(PYTHON_MODULE) test_spidir
