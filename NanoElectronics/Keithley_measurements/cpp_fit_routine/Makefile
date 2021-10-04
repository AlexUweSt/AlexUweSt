#CC=x86_64-w64-mingw32-g++
CC=g++
CFLAGS= -lboost_iostreams -lgsl -lgslcblas -lm  -lboost_system -lboost_filesystem -lboost_thread


all: fit run

fit:  fit_boost.cpp
	$(CC) -O3 -o  fit_boost.bin fit_boost.cpp $(CFLAGS)

run: fit_boost.bin
	./fit_boost.bin #GSL_IEEE_MODE = "round-to-zero"

#.PHONY: all run
