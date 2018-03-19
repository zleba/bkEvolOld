The program is using several ways to boost the computation speed.

When fitting the HERA data there are 3 elementary steps which can be computationally intense:
1) Calculation of the evolution kernel (initialPDF -> PDF)
2) Calculation of the convolution kernel (PDF -> F2, FL, sigmaR)
3) Applying convolution and evolution kernel (called several hundred of times when fitting)

The parallelization of the step 1) and 2) is trivial as the kernel cubes (3D matrices) can be calculated independently for each z (evolution kernel) and each Q2 (convolution kernel).
These kernels can be consequently stored into the files (currently HDF5 is used).
Fitting then represents only loading the kernels form hdf5 files and applying them multiple times to the inputPDF to find the minimum.

Parallelization applied in each step:
1) Evolution kernel
Currently openMP is coded as default to use all CPU cores.
In addition the openMPI parallelization is also implemented for systems with shared memory (~100 cores).  In case that only one CPU is available this parallelization has no effect.
The batch calculation of the evolution kernel can be easily implemented but as there is not large time issue, not used yet.

2) Convolution kernel
The calculation requires multidimensional integration which is quite time consuming.
To boost it, each q2 node is calculated as one job at batch system (~50 jobs).
Each such job employs openMP to make use of all cores of the machine.

3) Application of the kernels
This last step is purely about linear algebra but cannot be parallelized trivially.
The PDF always depends only on PDF at higher x-values.
Therefore, the next, lower x-value depends on the previous, higher, one.
To boost such calculation the special library for linear algebra (armadillo using BLAS procedures) is used.
The BLAS procedures internally use openMP and SIMD parallelization. 
Furthermore, the evaluation of discretized "integrals" for single x-value is parallelized by openMPI.

As an alternative, also the GPU-based kernel application is implemented.
In that case, the BLAS procedures are evaluated at GPU within CUDA.
It can speed up the calculation substantially but nVidia GPU is required.

Installation
The code is written in such a way to be able to be executed on usual laptop.
It requires:
-- ROOT, as Minuit is used for minimization.
-- HDF5 library for writing kernels to the files
-- openMPI libraries
