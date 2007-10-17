LIB = compbio.a
C_OBJS = dnatrie.o
CXX_OBJS = 
PROGS = grep-overlap

#CC=g++
#CXX=g++

CFLAGS = -Wall 

ifdef DEBUG
	CFLAGS := $(CFLAGS) -g
else
	CFLAGS := $(CFLAGS) -O3
endif

all: $(LIB) grep-overlap test_compbio 


# basic compile rule
$(C_OBJS): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(CXX_OBJS): %.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

test_compbio: $(LIB) test.cpp
	$(CXX) $(CFLAGS) test.cpp $(LIB) -o test_compbio 

grep-overlap: grep-overlap.c
	$(CC) -Wall -O3 grep-overlap.c -o grep-overlap

# basic library building
$(LIB): $(C_OBJS) $(CXX_OBJS)
	$(AR) -r $(LIB) $(C_OBJS) $(CXX_OBJS)


clean:
	rm -f $(PROGS) $(LIB) $(C_OBJS) $(CXX_OBJS)
