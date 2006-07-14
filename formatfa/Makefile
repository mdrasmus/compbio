#
# Scripting engine selection
#


CC = g++
AR = ar
COMPILE = -c -g

OBJS = formatfa.o

formatfff: $(OBJS)
	$(CC) $(OBJS) -o formatfa

# basic compile rule
$(OBJS): %.o: %.cpp
	$(CC) $(COMPILE) -o $@ $<

# basic library building
$(LIB): $(OBJS)
	$(AR) -r $(LIB) $(OBJS)

# basic clean rule
clean:
	rm -rf $(OBJS) $(LIB) $(PROG)
