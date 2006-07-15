#
# Scripting engine selection
#


CC = g++
AR = ar
COMPILE = -c -g

OBJS = formatfa.o verifyfa.o

all: formatfa verifyfa

formatfa: formatfa.o
	$(CC) formatfa.o -o formatfa

verifyfa: verifyfa.o
	$(CC) verifyfa.o -o verifyfa


# basic compile rule
$(OBJS): %.o: %.cpp
	$(CC) $(COMPILE) -o $@ $<

# basic library building
$(LIB): $(OBJS)
	$(AR) -r $(LIB) $(OBJS)

# basic clean rule
clean:
	rm -rf $(OBJS) $(LIB) $(PROG)
