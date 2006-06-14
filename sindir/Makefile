CC = g++
LIBS = 
PROG = sindir
OBJS = common.o \
       sinder.o \
       Tree.o
      

COMPILE = -g -Wall


all: $(PROG)




$(PROG): $(OBJS) 
	$(CC) -o $(PROG) $(OBJS) $(LIBS)

# basic compile rule
$(OBJS): %.o: %.cpp
	$(CC) $(COMPILE) -o $@ -c $<

# basic clean rules
clean:
	rm -rf $(OBJS) $(LIB) $(PROG)

cleanobj:
	rm -rf $(OBJS) $(LIB)
