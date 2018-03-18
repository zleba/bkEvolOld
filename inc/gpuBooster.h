#ifndef _GPUBOOSTER_
#define _GPUBOOSTER_

#include <cublas_v2.h>
#include <armadillo>
#include <vector>
#include <utility>
#include <cassert>

#include "utils.h"


using namespace std;



typedef double  Double;

struct Unit {
    size_t n;
    size_t nr, nc, ns;//nrows, ncols
    int devId;
    cublasHandle_t handle;
    Double *d_Cube;
    Double *d_phi;

    void GetResult(int y, arma::vec &res);
    void SetPhi(int y, arma::vec &phi);

};

class gpuBooster {

    public:
        bool isInited = false;
    vector<Unit> devs;


        void Init(const arma::cube &cube);
        void Convolute(int y);

        void InitAll(const arma::cube &cubeEvol, const arma::cube &cubeF2, const arma::cube &cubeFL);
        void ResetVector();
        void ConvoluteAll(int y);

        void GetResults(int y, arma::vec &pdf, arma::vec &f2, arma::vec &fl) {

            vector<arma::vec> temp(devs.size(), arma::vec(devs[0].nr));

            #pragma omp parallel
            {
                int i =omp_get_thread_num();
                assert(i < 4);
                //arma::vec temp(devs[i].nc);
                devs[i].GetResult(y, temp[i]);

            }

            for(int i = 1; i < devs.size(); ++i)
                temp[0] += temp[i];

            //cout << "size is " << pdf.n_rows << endl;
            pdf += temp[0](arma::span(0,pdf.n_rows-1));
            f2 += temp[0](arma::span(pdf.n_rows, pdf.n_rows + f2.n_rows-1 ));
            fl += temp[0](arma::span(pdf.n_rows+f2.n_rows, pdf.n_rows+f2.n_rows+fl.n_rows-1));


            //for(arma::vec &t : temp)
                //res += t;


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




#endif
