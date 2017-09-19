#ifndef _GPUBOOSTER_
#define _GPUBOOSTER_

#include <cublas_v2.h>
#include <armadillo>
#include <vector>
#include <utility>
#include <cassert>



using namespace std;



typedef double  Double;

struct Unit {
    size_t n;
    size_t nr, nc;//nrows, ncols
    int devId;
    cublasHandle_t handle;
    Double *d_Cube;
    Double *d_phi;

    void GetResult(int y, arma::vec &res);
    void SetPhi(int y, arma::vec &phi);

};

class gpuBooster {

    public:
    vector<Unit> devs;


        void Init(const arma::cube &cube);
        void Convolute(int y);

        void InitAll(const arma::cube &cubeEvol, const arma::cube &cubeF2, const arma::cube &cubeFL);
        void ConvoluteAll(int y);

        void GetResult(int y, arma::vec &res) {

            vector<arma::vec> temp(devs.size(), arma::vec(devs[0].nc));

            #pragma omp parallel
            {
                int i =omp_get_thread_num();
                assert(i < 4);
                //arma::vec temp(devs[i].nc);
                devs[i].GetResult(y, temp[i]);

            }


            for(arma::vec &t : temp)
                res += t;


            /*
            for(auto &dev : devs) {
                dev.GetResult(y, temp[0]);
                res += temp[0];
            }
            */

        }
        void SetPhi(int y, arma::vec &phi) {
            //cout <<"nDevices RADEK " << devs.size() << endl;
            /*
            for(auto &dev : devs) {
                //cout << "DeviceID : "<<dev.devId << endl;
                dev.SetPhi(y, phi);
            }
            */

            #pragma omp parallel
            {
                int i =omp_get_thread_num();
                devs[i].SetPhi(y, phi);
            }

        }
};


inline pair<long long,long long> GetStartEnd(int nrank, int rank, long long Min, long long Max)
{
    if(Max < Min) return make_pair(Min, Max);

    long long n = Max - Min + 1;

    //long long start = rank/(nrank+0.) * n;
    //long long end = (rank+1)/(nrank+0.) * n- 1;

    long long start, end;

    if(nrank <= n) {
        start =  (rank*n)   /nrank;
        end   = ((rank+1)*n)/nrank - 1;
    }
    else {
        if(rank < n)
            start = end = rank;
        else {
            start = 0;
            end = -1;
        }

    }
    start += Min;
    end += Min;

    return make_pair(start, end);

}


#endif
