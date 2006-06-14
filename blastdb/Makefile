CC = g++
LIBS = -ldb_cxx


blastdb: blastdb.cpp
	$(CC) $(LIBS) blastdb.cpp -o blastdb

install: blastdb
	cp blastdb ../bin

clean:
	rm -rf blastdb
