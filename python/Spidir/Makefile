extension.so: extension.cpp common.cpp common.h likelihood.cpp 
	g++ -shared -g \
                -I/usr/include/python2.4 \
                -I/util/include/python2.4 \
                extension.cpp common.cpp likelihood.cpp -o extension.so
