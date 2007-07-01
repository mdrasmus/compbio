CC = g++
MODULE = pyspidir.so
PROG_OBJS = \
    branchlen.o \
    common.o \
    likelihood.o \
    parsimony.o \
    search.o \
    Tree.o \

MODULE_OBJS = \
    pyspidir.o \
    $(PROG_OBJS)



CFLAGS = \
    -g \
    -Wall \
    -I/usr/include/python2.4 \
    -I/util/include/python2.4

all: $(MODULE) test_spidir



$(MODULE): $(MODULE_OBJS)
	g++ -shared $(MODULE_OBJS) -o $(MODULE)

test_spidir: $(PROG_OBJS) test.cpp
	g++ $(PROG_OBJS) -g test.cpp -o test_spidir


# basic compile rule
$(MODULE_OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<


install: test_spidir
	cp test_spidir ../../bin/

clean:
	rm -rf $(MODULE_OBJS) $(MODULE) test_spidir
