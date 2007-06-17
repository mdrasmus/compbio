pyspidir.so: pyspidir.cpp common.cpp common.h likelihood.cpp 
	g++ -shared -O3 \
                -I/usr/include/python2.4 \
                -I/util/include/python2.4 \
                pyspidir.cpp common.cpp likelihood.cpp -o pyspidir.so
