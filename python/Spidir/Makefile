CC = g++
MODULE = pyspidir.so
OBJS = \
    pyspidir.o \
    common.o \
    likelihood.o \
    parsimony.o \
    search.o \
    Tree.o \


CFLAGS = \
    -g \
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
