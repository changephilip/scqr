CXX=g++
CC=gcc

CXXFLAGS=--std=c++11 -g -fopenmp -lzstd -lm -lz  -mpopcnt

RFLAGS=--std=c++11  -g -fopenmp -lz -mpopcnt -DDEBUG

TOP_DIR=.
INCLUDE_DIR=-I$(TOP_DIR)/include
GZFILES = include/gzclose.c include/gzlib.c include/gzread.c include/gzwrite.c

main:
	#$(CC) -c include/zstd_zlibwrapper.c include/zstd_zlibwrapper.h -DZWRAP_USE_ZSTD=1
	#$(CC) -c $(GZFILES)
	#$(CXX) src/main.cpp  gzclose.o gzlib.o gzread.o gzwrite.o zstd_zlibwrapper.o -o bin/scqr $(CXXFLAGS) $(INCLUDE_DIR) 
	$(CXX) src/main.cpp   -o bin/scqr $(RFLAGS) $(INCLUDE_DIR) 

.PHONY: clean
clean:
	rm *.o
	rm bin/scqr
