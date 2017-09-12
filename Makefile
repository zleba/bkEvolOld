CFLAGS = $(shell root-config --cflags)
LIBS := $(shell root-config --libs)

HELPER=/home/radek/Dropbox/patrick/undolding/PlottingHelper
CC=g++

#iter: iterate.cpp
	#g++ -g -O3 $(CFLAGS)  $^ $(LIBS) -lgsl -lgslcblas -fopenmp -o $@ 


SRCS = src/iterate.cpp src/kernels.cpp src/main.cpp src/integration.cpp
OBJS = obj/iterate.o obj/kernels.o obj/main.o  obj/integration.o


obj/%.o: src/%.cpp inc/Solver.h inc/integration.h
	$(CC) -g  -fopenmp  -I./inc -c   -o $@ $< $(CFLAGS)

iter: $(OBJS) 
	$(CC) -g $^ $(LIBS)  -fopenmp -larmadillo  -o $@
