ARMAlib=/afs/desy.de/user/z/zlebcr/h1/TMD/Krakow/armadillo/armatemp
ARMAinc=/afs/desy.de/user/z/zlebcr/h1/TMD/Krakow/armadillo/armatemp/tmp/include/

convol: convol.cpp ../src/alphaSpline.cpp
	g++47 -fopenmp -std=c++11 -g -O3 -I$(ARMAinc) -I../inc $^  -Wl,-R$(ARMAlib)  -L$(ARMAlib)  -larmadillo -o $@

merge: merge.cpp
	g++47 -fopenmp -std=c++11 -g -O3 -I$(ARMAinc)  $<  -Wl,-R$(ARMAlib)  -L$(ARMAlib)  -larmadillo -o $@
