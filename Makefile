CFLAGS = $(shell root-config --cflags)
LIBS := $(shell root-config --libs)


armaLib=/home/zlebcr/libs/arma/armalib/lib64
armaInc=/home/zlebcr/libs/arma/armalib/include

HELPER=/home/radek/Dropbox/patrick/undolding/PlottingHelper
CC=g++
CC=mpic++

#iter: iterate.cpp
	#g++ -g -O3 $(CFLAGS)  $^ $(LIBS) -lgsl -lgslcblas -fopenmp -o $@ 


SRCS = src/iterate.cpp src/kernels.cpp src/main.cpp src/integration.cpp
OBJS = obj/iterate.o obj/kernels.o obj/main.o  obj/integration.o


#obj/%.o: src/%.cpp inc/Solver.h inc/integration.h
	#$(CC) -g  -fopenmp  -I$(armaInc) -DARMA_DONT_USE_WRAPPER   -I./inc -c   -o $@ $< $(CFLAGS)
#iter: $(OBJS) 
	#$(CC) -g  $^ $(LIBS)  -fopenmp  /opt/anaconda/3/pkgs/mkl-2017.0.3-0/lib/libmkl_rt.so -o $@

obj/%.o: src/%.cpp inc/Solver.h inc/integration.h
	$(CC) -g  -fopenmp  -I$(armaInc)   -I./inc -c   -o $@ $< $(CFLAGS)
iter: $(OBJS) 
	$(CC) -g  $^ $(LIBS)  -fopenmp   -Wl,-R$(armaLib) -L$(armaLib) -larmadillo  -o $@
