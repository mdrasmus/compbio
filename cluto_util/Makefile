CC = g++
CPPFLAGS = -I. -g
LIBS = -L. -lcluto
STATS_OBJS = cluto_stats.o cluto_io.o
REORDER_OBJS = cluto_reorder.o cluto_io.o cluto_tree.o
REORDER_TREE_OBJS = cluto_reorder_tree.o cluto_io.o cluto_tree.o
GRAPH_OBJS = cluto_graphify.o cluto_io.o 
APPS = cluto_stats cluto_reorder cluto_graphify 


all: $(APPS)

cluto_stats: $(STATS_OBJS)
	$(CC) $(STATS_OBJS) $(CPPFLAGS) $(LIBS) -o cluto_stats

cluto_reorder: $(REORDER_OBJS)
	$(CC) $(REORDER_OBJS) $(CPPFLAGS) $(LIBS) -o cluto_reorder

cluto_reorder_tree: $(REORDER_TREE_OBJS)
	$(CC) $(REORDER_TREE_OBJS) $(CPPFLAGS) $(LIBS) -o cluto_reorder_tree

cluto_graphify: $(GRAPH_OBJS)
	$(CC) $(GRAPH_OBJS) $(CPPFLAGS) $(LIBS) -o cluto_graphify

clean:
	rm -f $(APPS) *.o
