#ifndef _Fitter_
#define _Fitter_

#include <vector>
#include <armadillo>


using namespace std;

struct dataPoint {
    double x, Q2, y;
    double sigma, err;

    int q2ID;
    double theor;
};

class Fitter {
    vector<dataPoint> data;

    public:
    vector<dataPoint> LoadData(string fname);

    void AddTheory(vector<arma::vec> &f2, vector<arma::vec> &fl);
    void Init();

    double getChi2();

    void DoFit();

};

#endif
