#ifndef _Fitter_
#define _Fitter_

#include <vector>
#include <armadillo>

#include "Solver.h"

using namespace std;

struct dataPoint {
    double x, Q2, y;
    double sigma, err;

    int q2ID;
    double theor;
};

class Fitter {
    vector<dataPoint> data;
    Solver solver;

    public:

    Fitter() : solver(512+1) {}
    vector<dataPoint> LoadData(string fname);

    void AddTheory(vector<arma::vec> &f2, vector<arma::vec> &fl);
    void Init();

    double getChi2(int &nDF);
    double operator()(const double *p, const double *q);
    double Eval(const double *p, const double *q);

    void DoFit();

};

#endif
