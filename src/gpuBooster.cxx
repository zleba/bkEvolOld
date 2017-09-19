
/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>

//#include <helper_cuda.h>

#include <iostream>
#include <armadillo>
#include <vector>
#include <iomanip>
#include <cassert>

using namespace std;

#include "gpuBooster.h"

//typedef double  Double;


Double * CopyCube(const arma::cube &cube)
{
    Double *d_first;

    size_t size = size_t(cube.n_rows)*size_t(cube.n_cols)*size_t(cube.n_slices) * sizeof(*cube.begin());

    cout << "RADEK size " << size << endl;
    //exit(0);
    if (cudaMalloc((void **)&d_first, size) != cudaSuccess) {
        fprintf(stderr, "!!!! device memory allocation error (allocate A)\n");
        exit(1);
    }

    //cublasStatus_t status = cublasSetVector(cube.n_elem, sizeof(Double), cube.memptr(), 1, d_first, 1);
    cudaError_t status = cudaMemcpy(d_first, cube.memptr(), size, cudaMemcpyHostToDevice);


    //if (status != CUBLAS_STATUS_SUCCESS) {
    if (status != cudaSuccess) {
        fprintf(stderr, "!!!! device access error (write A)\n");
        exit(1);
    }
    return d_first;
}


Double * CopyAll(const arma::cube &cubeEvol, const arma::cube &cubeF2, const arma::cube &cubeFL)
{
    Double *d_first;

    auto GetSize = [](const arma::cube &cube) {
        return size_t(cube.n_rows)*size_t(cube.n_cols)*size_t(cube.n_slices) * sizeof(*cube.begin());
    };
    size_t size = GetSize(cubeEvol) + GetSize(cubeF2) + GetSize(cubeFL);

    size_t sliceSize = cubeEvol.slice(0).n_elem + cubeF2.slice(0).n_elem + cubeFL.slice(0).n_elem;
    size_t rowSize = cubeEvol.n_rows + cubeF2.n_rows + cubeFL.n_rows;


    cout << "RADEK size " << size << endl;
    //exit(0);
    if (cudaMalloc((void **)&d_first, size) != cudaSuccess) {
        fprintf(stderr, "!!!! device memory allocation error (allocate A)\n");
        exit(1);
    }

    //cublasStatus_t status = cublasSetVector(cube.n_elem, sizeof(Double), cube.memptr(), 1, d_first, 1);
    //cudaError_t status = cudaMemcpy(d_first, cube.memptr(), size, cudaMemcpyHostToDevice);

    for(int y = 0; y < cubeEvol.n_slices; ++y) 
    for(int i = 0; i < cubeEvol.n_cols; ++i) {
        double *base = d_first+y*sliceSize + i*rowSize;
        assert(!cudaMemcpy(base, cubeEvol.slice(y).colptr(i),
                   cubeEvol.n_rows*sizeof(*cubeEvol.begin()), cudaMemcpyHostToDevice));

        assert(!cudaMemcpy(base + cubeEvol.n_rows, cubeF2.slice(y).colptr(i),
                   cubeF2.n_rows*sizeof(*cubeF2.begin()), cudaMemcpyHostToDevice));
        assert(!cudaMemcpy(base + cubeEvol.n_rows+cubeF2.n_rows, cubeFL.slice(y).colptr(i),
                   cubeFL.n_rows*sizeof(*cubeFL.begin()), cudaMemcpyHostToDevice));
    }

    return d_first;
}







void gpuBooster::Init(const arma::cube &cube) {
    int nDev;
    assert(cudaGetDeviceCount(&nDev) == 0);
    //nDev = 1;
    assert(nDev > 0);
    devs.resize(nDev);

    //cout << "nDevices is " << nDev << endl;

    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(nDev);

    //Init memory in devices
    #pragma omp parallel
    //for(int i = 0; i < nDev; ++i)
    {
        int i =omp_get_thread_num();
        devs[i].devId = i;
        assert(cudaSetDevice(i+0) == 0);
        cublasStatus_t status = cublasCreate(&devs[i].handle);
        assert(status == CUBLAS_STATUS_SUCCESS);

        //cout << "My thread id is " << i << endl;


        //Move cube to GPU
        devs[i].d_Cube = CopyCube(cube);
        devs[i].n = cube.n_rows;

        //Init sol in GPU to zero
        size_t size = size_t(cube.n_rows) * size_t(cube.n_slices) * sizeof(*cube.begin());
        if (cudaMalloc((void **)&devs[i].d_phi, size) != cudaSuccess) {
            fprintf(stderr, "!!!! device memory allocation error (phi)\n");
            exit(1);
        }
        //cout << "Initialized GPU mem from " << devs[i].d_phi << " "<< cube.n_rows * cube.n_slices << endl;
        const Double zero = 0, one = 1;
        assert(cublasDscal(devs[i].handle, cube.n_rows * cube.n_slices, &zero, devs[i].d_phi, 1) == 0);

    }
}


void gpuBooster::Convolute(int y) {

    const Double one = 1;

    int nDev = devs.size();
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(nDev);

    #pragma omp parallel
    {
        int id =omp_get_thread_num();
        int start, end;
        tie(start,end)= GetStartEnd(devs.size(), id, 1, y-1);

        int n = devs[id].n;
        //cout << "Id " << id << " " << start << " "<< end << endl;

        for(int d = start; d <= end; ++d) {
            //yTemp += matN.slice(d) * PhiRapN[y-d];
            cublasStatus_t status = cublasDgemv(devs[id].handle, CUBLAS_OP_N, n, n, &one, devs[id].d_Cube + n*n*d, n, devs[id].d_phi + n*(y-d), 1, &one, devs[id].d_phi + n*y, 1);
            assert(status == CUBLAS_STATUS_SUCCESS);
        }
    }
}

void Unit::GetResult(int y, arma::vec &res) {
    //cout << "Device id " << devId << endl;
    assert(cudaSetDevice(devId+0) == 0);
    //cout << "Reading GPU mem from " << d_phi <<" "<< d_phi /*+ n*y */ <<", size " << res.n_elem <<" "<<n<< endl;
    cublasStatus_t status = cublasGetVector(n, sizeof(Double), d_phi + n*y, 1, res.memptr(), 1);
    //cout << "Status read" << status << endl;
    assert(status == CUBLAS_STATUS_SUCCESS);
}

void Unit::SetPhi(int y, arma::vec &phi) {
    
    //cout << "RADEKK " << devId << endl;
    assert(cudaSetDevice(devId+0) == 0);
    cublasStatus_t status = cublasSetVector(n, sizeof(Double), phi.memptr(), 1, d_phi +n*y, 1);
    //cout << "RADEKK stat write " <<status <<" "<< devId << endl;
    assert(status == CUBLAS_STATUS_SUCCESS);
}




template<typename T>
vector<Double*> TransferObjects(vector<T> &mats)
{
    Double *d_first;
    long long nSingle = mats[0].n_elem;
    long long size = mats[0].n_elem*sizeof(Double) *mats.size();
    if (cudaMalloc((void **)&d_first, size) != cudaSuccess) {
        fprintf(stderr, "!!!! device memory allocation error (allocate A)\n");
        exit(1);
    }

    vector<Double*> ps;
    for(int i = 0; i < mats.size(); ++i) {
        /* Initialize the device matrices with the host matrices */
        cublasStatus_t status = cublasSetVector(nSingle, sizeof(Double), mats[i].memptr(), 1, d_first+nSingle*i, 1);
        if (status != CUBLAS_STATUS_SUCCESS) {
            fprintf(stderr, "!!!! device access error (write A)\n");
            exit(1);
        }
        ps.push_back(d_first+nSingle*i);
    }
    return ps;
}




/* Main */
//int main(int argc, char **argv)
//{ return 0; }
