CC = g++
MODULE = pyspidir.so
OBJS = \
    pyspidir.o \
    common.o \
    likelihood.o \
    search.o \
    parsimony.o

CFLAGS = \
    -O3 \
    -Wall \
    -I/usr/include/python2.4 \
    -I/util/include/python2.4

$(MODULE): $(OBJS)
	g++ -shared $(OBJS) -o $(MODULE)


# basic compile rule
$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) -o $@ $<


clean:
	rm -rf $(OBJS) $(MODULE)
