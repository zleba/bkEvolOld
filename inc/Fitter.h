#ifndef _Fitter_
#define _Fitter_

struct dataPoint {
    double x, Q2, y;
    double sigma, err;

    int q2ID;
    double theor;
};

class Fitter {
    vector<dataPoint> data;

    public:
    vector<dataPoint> LoadData();

    void AddTheory(vector<arma::vec> &f2, vector<arma::vec> &fl);
    void Init();

    double getChi2();


};

#endif
