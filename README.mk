The program is using several ways to boost the the computation speed.

When fitting the HERA data there are 3 elementary steps which can be cumputionally intence:
1) Calculation of the evolution kernel (initialPDF -> PDF)
2) Calculation of the convolution kernel (PDF -> F2, FL, sigmaR)
3) Applying convolution and evolution kernel (called several hundered of times when fitting)

The paralelisation of the step 1) and 2) is trivial as the kernel cubes (3D matrixes) can be calculated independently for each z (evolution kernel) and each Q2 (convolution kernel).
These kernels can be consequently stored into the files (currently HDF5 is used).
Fitting then represents only loading the kernels form hdf5 files and applying them multiple times to the inputPDF to find the minimum.

Paralelisation applyed in each step:
1) Evolution kernel
Currently openMP is coded as default to use all CPU corres.
In addition the openMPI paralelisation is also impremented for systems with shared memory (~100 cores).  In case that only one CPU is availabe this paralelisation has no effect.
The batch caculation of the evolution kernel can be easily implremented but as there is not large time issue, not used yet.

2) Convolution kernel
The calculation requires multidimensional integration which is quite time consuming.
To boost it, each q2 node is calculated as one job at batch system (~50jobs).
Each such job employs openMP to make use of all cores of the machine.

3) Application of the kernels
This last step is purely about linear algebra but cannot be paralised trivialy.
The PDF always depends only on PDF at higher x-values.
Therefore, the next, lower x-value depens on the privous, higher, one.
To boost such calculation the special library for linear algebra (armadillo using BLAS precedures) is used.
The BLAS procedures internaly use openMP and SIMD paralelelisations. 
Futhermore, the evaluation of PDF for single x-value is paralised by openMPI.

As an alternative, also the GPU-based kernel application is implemented.
In this case, the BLAS procedures are evaluated at GPU within CUDA.
It can speed up the calculation substantionaly but nVidia GPU is required.

Instalation
The code is writen in such a way to be able to be executed on standard PC.
