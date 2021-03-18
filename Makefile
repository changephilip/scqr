CC=g++

CXXFLAGS=--std=c++11 -g -fopenmp -lz -lm 

TOP_DIR=.
INCLUDE_DIR=-I$(TOP_DIR)/include

main:
	$(CC) src/main.cpp -o bin/main $(CXXFLAGS) $(INCLUDE_DIR)

.PHONY: clean
clean:
	rm bin/main
