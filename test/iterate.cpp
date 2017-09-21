#include <iostream>
#include <functional>
#include <vector>
#include <iomanip>

//#include <gsl/gsl_chebyshev.h>

#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <armadillo>

using namespace std;

vector<double> GetWeights(int Size)
{
    //const int Size = 33;
    const int N = Size - 1;
    assert(N%2 == 0);

    vector<vector<double>> coef(Size);
    for(auto & c: coef) c.resize(Size);

    for(int k = 0; k <= N/2; ++k) {
        double s = 0;
        coef[2*k][N] = 1./N;
        coef[2*k][0] = 1./N ;

        coef[2*k][N/2] = 2./N * (2* ((k+1)%2) -1);

        for(int n = 1; n <= N/2-1; ++n)
            coef[2*k][n] = coef[2*k][N-n] = 2./N * cos(n*k*M_PI*2/N);
    }

    vector<double> wgt(Size, 0.);

    for(int i = 0; i < Size; ++i) {
        wgt[i] += coef[0][i];
        wgt[i] += coef[N][i] / (1. - N*N);
        for(int k = 1; k <= N/2 - 1; ++k) {
           double w = 2./(1 - 4*k*k);
           wgt[i] += w*coef[2*k][i];
        }
    }
    return wgt;
    //for(int i = 0; i < Size; ++i)
        //cout << i <<" "<<setprecision(17)<< wgt[i] << endl;
}

arma::mat GetCoefs(int oldSize, bool isInverse = false)
{
    const int N = oldSize - 1;
    assert(N%2 == 0);

    arma::mat coef(oldSize,oldSize);

    double mul = 1;
    double C = 1./N;
    if(isInverse == true) {C = 1./2; }

    //isInverse = false;
    for(int k = 0; k <= N; ++k) {
        double s = 0;
        if(!isInverse) {
            coef(k,N) = C;
            coef(k,0) = C * (k % 2 == 1 ? -1 : 1);
        }
        else {
            mul = k % 2 == 1 ? -1 : 1;
            coef(N-k, N) = C * mul;
            coef(N-k, 0) = C ;
        }

        for(int n = 1; n <= N-1; ++n) {
            double el = cos(n*k*M_PI / N) * 2.*C * mul;
            if(!isInverse) coef(k,N-n) = el;
            else           coef(N-k,N-n) = el;
        }
    }
    
    return coef;
}

arma::vec GetWeights(const arma::mat &coef)
{
    int Size = coef.n_rows;
    int N = Size - 1;
    arma::vec wgt(Size, arma::fill::zeros);

    for(int i = 0; i < Size; ++i) {
        wgt(i) += coef(0,i);
        wgt(i) += coef(N,i) / (1. - N*N);
        for(int k = 1; k <= N/2 - 1; ++k) {
           double w = 2./(1 - 4*k*k);
           wgt(i) += w*coef(2*k,i);
        }
    }
    return wgt;
}








arma::vec GetVals(int n, function<double(double)> fun) 
{
    arma::vec vals(n);
    //vector<double> vals;
    for(int i = n-1; i >= 0; --i) {
        double x = cos(i /(n-1.) * M_PI);
        //vals.push_back(fun(x));
        vals(n-1 - i) = fun(x);
    }
    return vals;
}

double Eval(arma::vec coefs, double theta)
{
    double sum = 0;

    //sum = coefs[0]/2;

    //cout << "RADEK 0 " << coefs[0]/2 << endl;
    for(int i = 0; i < coefs.size(); ++i) {
        double fact = (i==0 || i == coefs.size()-1) ? 1./2 : 1;
        //int fact = 1;
        sum += fact * coefs(i) * cos(theta*i);
        //cout << "RADEK "<<i<<" " << fact * coefs[i] * cos(theta*i) << endl;
    }

    return sum;
}

arma::mat Proj(int to, int from)
{
    arma::mat proj(to, from, arma::fill::eye); 
    if(from > to)
        proj(to-1, to-1) *= 2;
    else if(from < to)
        proj(from-1, from-1) *= 1./2;

    return proj;
}

int main()
{
    auto fun = [](double x) { /*double x = acos(CosX);*/ return  /*1+abs(x);*/cos(4*acos(x)); };


    auto vals1 = GetVals(5, fun);
    auto vals2 = GetVals(pow(2,3)+ 1, fun);

    auto coefMat1 = GetCoefs(vals1.size());
    auto coefMat2 = GetCoefs(vals2.size());

    auto coefMat1in = GetCoefs(vals1.size(), true);
    auto coefMat2in = GetCoefs(vals2.size(), true);
    //arma::mat both = coefMat1in * coefMat1;
    //cout << both << endl;
    //cout << both(1,1) << endl;
    //return 0;

    arma::vec res = GetWeights(coefMat2);
    for(auto a : res)
        cout << a << endl;

    return 0;

    arma::vec coef1 = Proj(vals2.size(),vals1.size())*coefMat1 * vals1;
    arma::vec coef2 = coefMat2 * vals2;

    //convert to higher
    auto matExt = coefMat2in * Proj(vals2.size(),vals1.size()) * coefMat1;
    auto matRed = coefMat1in * Proj(vals1.size(),vals2.size()) * coefMat2;

    for(auto a: vals1)
        cout << a << endl;
    cout << endl;
    for(auto a: vals2)
        cout << a << endl;

    cout << endl << endl;

    cout << "Coefs 1" << endl;
    for(auto a: coef1)
        cout << a << endl;
    cout << endl;
    cout << "Coefs 2" << endl;
    for(auto a: coef2)
        cout << a << endl;

    auto Print = [&](arma::vec coef, int n) {
        //int n = coef.size();
        for(int i = n-1; i >= 0; --i) {
            double x = cos(i /(n-1.) * M_PI);
            cout << x << " "<<Eval(coef, acos(x)) << " " << fun(x) << endl;
        }
    };
    cout << "In nodes 1" << endl;
    Print(coef1, coef1.size());
    cout << "In nodes 2" << endl;
    Print(coef2, coef2.size());
    cout << " Vals 2 my" << endl;
    auto vals2my = coefMat2in * coef2;
    cout << vals2my << endl;

    //cout << Eval(coef1, 0.5) <<" "<< cos(0.5) << endl;


    return 0;
}
