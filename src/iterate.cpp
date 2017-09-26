#include <iostream>
#include <functional>
#include <vector>
#include <iomanip>

//#include <gsl/gsl_chebyshev.h>

//#include "TLegend.h"

#include "TGraph.h"
//#include "TCanvas.h"
//#include "TAxis.h"
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
//#include "TGraph2D.h"
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

vector<vector<double>> GetTransformMatrix(int oldSize, int newSize)
{
    assert((newSize - 1) % (oldSize - 1) == 0);

    const int N = oldSize - 1;
    assert(N%2 == 0);

    vector<vector<double>> coef(oldSize);
    for(auto & c: coef) c.resize(oldSize);

    for(int k = 0; k <= N; ++k) {
        double s = 0;
        coef[k][N] = 1./N;
        coef[k][0] = 1./N * (k % 2 == 1 ? -1 : 1);

        for(int n = 1; n <= N-1; ++n)
            coef[k][n] = cos(n*k*M_PI / N) * 2./N;
    }
    
    //GetNew values


        
    const int Nnew = newSize - 1;
    assert(Nnew%2 == 0);

    for(int n = 0; n <= N; ++n) {
        for(int i = 0; i < Nnew; ++i) {
            double theta = i /(Nnew-1.) * M_PI;

            double sum = coef[0][n]/2;
            for(int k = 1; k <= N; ++k)
                sum += coef[k][n] * cos(k* theta);

          //cout << "R " << exp(0.5*xi[i]) << endl;
        }
    }

    //coef

}

map<double, TGraph*> ReadFile(const char *fName)
{
    ifstream file(fName);

    map<double, TGraph*> grMap;
    int i = 0;
    while(1) {
        string str;
        getline(file, str);
        if(!file.good()) break;
        if(str[0] == '#') continue;
        if(str.size() < 8) continue;
        stringstream  sStream(str);
        double kT2, b, y, Phi;
        sStream >> kT2 >> b >> y >> Phi;
        kT2 *= kT2;

        if(grMap.count(y) == 0) {
            i = 0;
            grMap[y] = new TGraph();
        }
        grMap[y]->SetPoint(i++, kT2, Phi);
        //cout << kT2 <<" "<< b<<" "<<y << " "<< Phi << endl;

    }
    return grMap;

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

arma::mat Proj(int to, int from)
{
    arma::mat proj(to, from, arma::fill::eye); 
    if(from > to)
        proj(to-1, to-1) *= 2;
    else if(from < to)
        proj(from-1, from-1) *= 1./2;

    return proj;
}

pair<arma::mat, arma::mat> GetTransMatrices(int base, int ext, bool toTrivial)
{

    arma::mat coefMat1 = GetCoefs(base);
    arma::mat coefMat2 = GetCoefs(ext);

    arma::mat coefMat1in = GetCoefs(base, true);
    arma::mat coefMat2in = GetCoefs(ext, true);

    arma::mat matExt = coefMat2in * Proj(ext, base) * coefMat1;

    arma::mat matRed;
    if(!toTrivial)
        matRed = coefMat1in * Proj(base,ext) * coefMat2;
    else {
        matRed = arma::mat(base, ext, arma::fill::zeros);
        int fac = (ext-1)/(base-1);
        for(int i = 0; i < base; ++i) 
            matRed(i, i*fac) = 1;
    }

    return make_pair(matRed, matExt);
}
