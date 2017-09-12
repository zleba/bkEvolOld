
#include <armadillo>

#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <algorithm>
#include <iomanip>

int N = 5000;

using namespace std;

void armaTest()
{
    arma::mat M1(N,N), M2(N,N);
    for(int i = 0; i < N; ++i)
    for(int j = 0; j < N; ++j){
        M1(i,j) = 1.1;
        M2(i,j) = 0.1;
    }

    //arma::mat res = M1 * M2;
    //cout << res(0,0) << " "<< res(5,8) << endl;
    //return;

    arma::vec v(N);
    for(int i = 0; i < N; ++i) v(i) = i;

    cout << "Start " << endl;
    cout << "End " << endl;
    double sum = 0;
    for(int i = 0; i < N; ++i) {
        arma::vec rr = M1 * v;
        sum += rr(0);
    }
    cout << sum << endl;

}

void vectorTest()
{
    vector<vector<double>> M1(N), M2(N), res(N);

    for(int i = 0; i < N; ++i) {
        M1[i].resize(N);
        M2[i].resize(N);
        for(int j = 0; j < N; ++j){
            M1[i][j] = 1.1;
            M2[i][j] = 0.1;
        }
    }

    /*
    cout << "Start" << endl;
    for(int i = 0; i < N; ++i) {
        res[i].resize(N,0.0);
        for(int j = 0; j < N; ++j){
            for(int k = 0; k < N; ++k)
                res[i][j] += M1[i][k]*M2[k][j];
        }
    }
    cout << "End" << endl;
    cout << res[0][0] << " "<< res[5][8] << endl;

    //return;
    */

    vector<double> v(N);
    for(int i = 0; i < N; ++i) v[i] = i;

    cout << "Start " << endl;
    //arma::mat res = M1 * M2;
    cout << "End " << endl;
    double sum = 0;
    for(int i = 0; i < N; ++i) {
        vector<double> rr(N);
        #pragma omp parallel for
        for(int k = 0; k < N; ++k)
        for(int j = 0; j < N; ++j)
            rr[k] += M1[k][j] * v[j];
        sum += rr[0];
    }
    cout << sum << endl;




}


int main()
{
    armaTest();
    //vectorTest();

    //cout << res << endl;

    return 0;
}

