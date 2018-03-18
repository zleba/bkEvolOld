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
    double theor0, extra0;
    double softFr;
};

class Fitter {
    vector<dataPoint> data;
    Solver solver;

    public:

    Fitter() : solver(512+1) {}
    Fitter(istream &Stream): solver(Stream)  {}
    vector<dataPoint> LoadData(string fname);

    void AddTheory(vector<arma::vec> &f2, vector<arma::vec> &fl);
    void Init(string dirName);

    double getChi2(int &nDF);
    double getChi2Corr(int &nDF);
    double operator()(const double *p, const double *q);
    double Eval(const double *p);

    void CalculateBasis(int nElem, string name);
    arma::mat getPoints();

};

#endif
