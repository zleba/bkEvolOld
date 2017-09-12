#ifndef Solver_H_
#define Solver_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <algorithm>
#include <iomanip>

#include <armadillo>
//using namespace arma;

using namespace std;

vector<double> GetWeights(int Size);
pair<arma::mat, arma::mat> GetTransMatrices(int base, int ext, bool toTrivial);

struct Nodes {
    Nodes(int _N, double _a, double _b) : N(_N), a(_a), b(_b) {}
    int N;
    double a, b;
    vector<double> xi, wi;
    void CalcNodes(bool bkSolverMode = true) {
        /*
        wi.resize(N);
        for(int i = 0; i < N; ++i) {
            double Cos = cos((2*i+1)/(2.*N) * M_PI);
            xi[i] = (a + b +  Cos * (a-b) )/2.;
            //cout << i<<" "<< xi[i] << endl;
            wi[i] = M_PI / N * sqrt(1-Cos*Cos) * (b-a)/2.;

            //double Cosp = cos((i+1.)/(N+1) * M_PI);
            //xi[i] = (a + b + Cosp * (a-b) )/2.;
            //wi[i] = M_PI / (N+1) * (1-Cosp*Cosp) / sqrt(1-Cosp*Cosp);
            //cout << (1-Cosp*Cosp) << endl;
            //assert(isfinite(wi[i]));
        }
        */

        //Classical nodes ala Clenshaw–Curtis
        if(!bkSolverMode) {
            xi.resize(N);
            for(int i = 0; i < N; ++i) {
              double Cos = cos(i /(N-1.) * M_PI);
              xi[i] = (a + b +  Cos * (a-b) )/2.;
              //cout << "R " << exp(0.5*xi[i]) << endl;
            }
            wi = GetWeights(N);
            for(auto &w : wi) w *= (b-a)/2.;

        }
        //Nodes used in BKsolver
        else {
            xi.resize(2*N-1);
            for(int i = 0; i < N*2-1; ++i) {
                double Cos = cos(i /(2*N-2.) * M_PI);
                double an = a - (b-a);
                double bn = b;
                xi[i] = (an + bn +  Cos * (an-bn) )/2.;
                //cout << "R " << exp(0.5*xi[i]) << endl;
            }

            xi.erase(xi.begin(), xi.begin() + N - 1);
            //cout << "Begin " << endl; 
            //for(auto x : xi)
                //cout << "ahoj " <<x <<" "<< exp(0.5*x) << endl;
            //cout << "End " << endl; 

            wi = GetWeights(2*N-1);
            wi.erase(wi.begin(), wi.begin() + N - 1);
            for(auto &w : wi) w *= 2*(b-a)/2.;

            wi.front() *= 0.5;
        }
        //exit(0);
    }
    double Integrate(function<double(double)> fun) {
        double sum = 0;
        for(int i = 0; i < N; ++i) {
           sum += fun(xi[i]) * wi[i];
        }
        return sum;
    }
    double Integrate(vector<double> vals) {
        assert(vals.size() == wi.size());
        reverse(vals.begin(), vals.end());
        double sum = 0;
        for(int i = 0; i < N; ++i) {
           sum += vals[i] * wi[i];
        }
        return sum;
    }



};



struct Solver {
    const double as = 0.2;
    const double eps = 1e-7;
    const int Nint; // kT nodes in Nintegral
    const int N;// = 32*16 + 1; //must be 2*n+1
    const int Nrap = 2*129;
    const bool toTrivial = true;

    const double Lmin= log(1e-2), Lmax = log(1e6);
    //const double Lmin= log(1e-4), Lmax = log(1e8);
    const double mu2 = 1e-2;
    const double rapMax = log(1e6), rapMin = log(1);
    bool putZero = false;

    Nodes nod, nodBase;
    
    Solver(int N_) : Nint(1*(N_-1)+1), N(N_),
        nod(Nint, Lmin, Lmax), nodBase(N, Lmin, Lmax) {
        nod.CalcNodes(false);
        nodBase.CalcNodes(false);
        tie(redMat,extMat) = GetTransMatrices(N, Nint, toTrivial);
    }

    vector<vector<vector<double>>> mat, matDiag, matTest;
    //vector<double> Phi0; //function to evolve
    vector<vector<double>> PhiRap, Phi0;

    vector<arma::mat> matN, matNDiag;
    vector<arma::vec> PhiRapN, Phi0N;

    arma::mat extMat, redMat;

    pair<double,double> GetKerPar(double l, double lp);
    double Delta(double z, double k2, double q2);
    //Equation 79
    double Kernel79(double l, double lp, double z);
    double Kernel79Diag(double l, double lp, double z);
    double Kernel80(double l, double lp, double z);
    double Kernel80Diag(double l, double lp, double z);
    double Kernel81(double l, double lp, double z);
    double Kernel81Diag(double l, double lp, double z);
    double Kernel83(double l, double lp, double z);
    double Kernel83Diag(double l, double lp, double z);
    double Kernel84(double l, double lp, double z);
    double Kernel84Diag(double l, double lp, double z);

    double DGLAPterm(double l, double lp, double z);

    double Kernel85(double l, double lp, double z);
    double Kernel85Diag(double l, double lp, double z);
    double Kernel85zDiag(double l, double lp, double x);

    double Kernel86(double l, double lp, double z);
    double Kernel86Diag(double l, double lp, double z);
    double Kernel86zDiag(double l, double lp, double x);


    double Kernel87(double l, double lp, double z);
    double Kernel87Diag(double l, double lp, double z);
    double Kernel88(double l, double lp, double z); //TODO not finish

    double KernelSub79(double l, double lp, double z);
    double KernelSub79Diag(double l, double lp, double z);
    double KernelSub80(double l, double lp, double z);
    double KernelSub80Diag(double l, double lp, double z);
    double KernelSub81(double l, double lp, double z);
    double KernelSub81Diag(double l, double lp, double z);
    double KernelSub83(double l, double lp, double z);
    double KernelSub83Diag(double l, double lp, double z);
    double KernelSub84(double l, double lp, double z);
    double KernelSub84Diag(double l, double lp, double z);
    double KernelSub87(double l, double lp, double z);
    double KernelSub87Diag(double l, double lp, double z);

    double KernelSub88(double l, double lp, double z);



    //Kernel org
    pair<double,double> KernelDiag(double L, double Lp) {
        double Exp = exp(L - Lp);
        double resSing, resReg;
        if(L != Lp)
            resSing = as * (- Exp/(abs(1. - Exp)) );
        else
            resSing = 0.;
        resReg  = as * (+ Exp/sqrt(4 + Exp*Exp));
        //resReg  = 0;
        return {resSing, resReg};
    }
    double Kernel(double L, double Lp) {
        double Exp = exp(L - Lp);
        double res = as / (abs(1 - Exp));
        return res;
    }


    void InitMat() {
        mat.resize(Nrap); //Nrap * N*N, evolution matrix
        for(auto & m : mat) {
            m.resize(N);
            for(auto & r: m)
                r.resize(N, 0.);
        }
        matTest = mat; //set to ZEROS
        matDiag = mat; //matDiag for the z-diagonal part (set to ZEROS)

        matN.resize(Nrap);
        matNDiag.resize(Nrap);

        int fac = (Nint-1)/(N-1);

        double stepY = (rapMax - rapMin) / (Nrap-1);
        for(int y = 0; y < Nrap; ++y) { 
            double rap = rapMin + stepY * y;

            arma::mat mTemp(Nint,Nint,arma::fill::zeros);
            arma::mat mDiagTemp(Nint,Nint,arma::fill::zeros);
            cout << "Matrix init y " << y << endl;

            #pragma omp parallel for
            for(int i = 0; i < Nint; ++i)   //loop over L
            for(int j = 0; j < Nint; ++j) { //loop over L' (integral)
                double L  = nod.xi[i];
                double Lp = nod.xi[j];
                double w = nod.wi[j];

                double l  = exp(0.5*L);
                double lp = exp(0.5*Lp);
                double z = exp(-rap);

                //mat[y][i][j] += Kernel85(l, lp, z) * w;
                //mat[y][i][i] += Kernel85Diag(l, lp, z) * w;

                //mat[y][i][j] += KernelSub81(l, lp, z) * w;
                //mat[y][i][i] += KernelSub81Diag(l, lp, z) * w;

                //mat[y][i][j] += KernelSub88(l, lp, z) * w;

                //matDiag[y][i][j] = Kernel85zDiag(l, lp, z) * w; //Diag-z DGLAP part

                if(!toTrivial || i % fac == 0) {
                    mTemp(i,j) += KernelSub79(l, lp, z) * w;
                    mTemp(i,i) += KernelSub79Diag(l, lp, z) * w;
                }


                //mTemp(i,j) += Kernel85(l, lp, z) * w;
                //mTemp(i,i) += Kernel85Diag(l, lp, z) * w;
                //mDiagTemp(i,j) = Kernel85zDiag(l, lp, z) * w;

                //Clasicall BFKL
                /*
                if(i != j)
                    mat[y][i][j] += Kernel(L, Lp) * w;

                double kerDs, kerDr;
                tie(kerDs, kerDr) = KernelDiag(L, Lp);
                if(i != j)
                    mat[y][i][i] += (kerDs+kerDr) * w;
                else
                    mat[y][i][i] += (kerDr) * w;
                */

            }
            //cout << redMat << endl;
            //cout << extMat << endl;

            matN[y]     = redMat * mTemp * extMat;
            matNDiag[y] = redMat * mDiagTemp * extMat;

            /*
            for(int i = 0; i < Nint; ++i)  //loop over L
            for(int j = 0; j < Nint; ++j) { //loop over L' (integral)
                double m  = matDiag[y][i][j];
                double mt = mDiagTemp(i,j);
                assert(isfinite(mat[y][i][j]));
                assert(isfinite(matDiag[y][i][j]));
                //cout << "Ihned "<<y<<" : " <<i<<" "<<j<<" "<< m <<" "<< mt << " "<< 2*(m-mt)/(m+mt) <<endl;
                //cout << "Mdiag "<<y<<" : " <<i<<" "<<j<<" "<< matDiag[y][i][j] <<endl;
            }
            */

            //if(y == 1) exit(0);

        }


        /*
        for(int k = 0; k < 4; ++k) {
            const double stepY = (rapMax - rapMin) / (Nrap-1);
            for(int i = 0; i < N; ++i)
            for(int j = 0; j < N; ++j) {
                cout <<"Kernel "<<k<<" : "<< i <<" "<< j <<" "<<setprecision(10)<< mat[k][i][j] <<" "<< matDiag[k][i][j]<<" "<< stepY<< endl;
            }
        }
        exit(0);
        */



        /*
        exit(0);
        */
    }

    vector<double> GetLinSolution(const vector<vector<double>> &MatEq, const vector<double> &y) {
        //cout 
        arma::mat M(N,N);
        for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j){
            M(i,j) = MatEq[i][j];
            //M(i,j) = -factor*Mat[i][j];
            //if(i == j) M(i,j) += 1;
        }
        arma::vec v(N);
        for(int i = 0; i < N; ++i)
            v(i) = y[i];
        arma::vec res = arma::solve(M,v);

        vector<double> Res(N);
        for(int i = 0; i < N; ++i)
            Res[i] = res[i];
        return Res;
    }

    arma::vec GetLinSolution(const arma::mat &Mat, const arma::vec &y) {
       return arma::solve(Mat, y); 
    }

    vector<double> GetRegSolution(const vector<vector<double>> &MatEq, const vector<double> &y, const vector<double> &yReg) {
        //cout 
        arma::mat M(N,N);
        for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j){
            M(i,j) = MatEq[i][j];
        }

        arma::vec v(N), vReg(N);
        for(int i = 0; i < N; ++i) {
           v(i) = y[i]; 
           vReg(i) = 0;//yReg[i]; 
        }
        arma::mat L(N-1,N);

        for(int i = 0; i < N-1; ++i) 
        for(int j = 0; j < N; ++j) {
            if(i == j)
                L(i,j) = 1.;///(pow(yReg[i],2) + 1e-10);
            if(i == j - 1)
                L(i,j) =-1.;///(pow(yReg[i],2) + 1e-10);
            else
                L(i,j) = 0.;///(pow(yReg[i],2) + 1e-10);
        }
        arma::mat L2 = trans(L)*L;

        //double tau = 3e-1;
        double tau = 5e-1;
        arma::mat E = arma::trans(M)*M + tau*L2;
        //for(int i = 0; i < N; ++i)
           //E(i,i) += tau * L(i,i);
    
        //cout << " " << Ep(3,3) << " "<< Ep(3,4) << endl;
        //arma::mat E = arma::trans(M)*M + tau;
        //cout << " " << E(3,3) << " "<< E(3,4) << endl;
        //assert(0);
        arma::vec r = arma::trans(M)*v + tau*L2*vReg;

        arma::vec res = arma::solve(E,r);

        vector<double> Res(N);
        for(int i = 0; i < N; ++i)
            Res[i] = res[i];
        return Res;

    }



    vector<double> IterSolution(const vector<vector<double>> &Mat, double factor, const vector<double> &y) {
        vector<double> yNow = y;
        vector<double> yTemp(N, 0.);
        vector<double> ySum = y;
        double diff = 0;
        for(int i = 0; i < 20; ++i) {
            std::fill(yTemp.begin(), yTemp.end(), 0.0);
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l)
                yTemp[k] += factor * Mat[k][l] * yNow[l];
            yNow = yTemp;

            diff = 0;
            for(int k = 0; k < N; ++k)
                diff = max(diff, abs(yNow[k]) / (1e-13+abs(yNow[k]) + abs(ySum[k])));
            //cout <<"Diff " <<  i << " "<< diff << endl;

            for(int k = 0; k < N; ++k)
                ySum[k] += yNow[k];

            //cout << "Diff is " << diff << endl;
            if(diff < 1e-9) break;
        }
        assert(diff < 1e-5);


        return ySum;
    }


    void EvolveNew() {
        const double stepY = (rapMax - rapMin) / (Nrap-1);

        int start = 0;

        //Classical approach
        if(start == 0) {
            arma::mat MatEq = arma::mat(N,N,arma::fill::eye) -  matNDiag[0];
            if(putZero) {
                PhiRapN[0] = arma::vec(N, arma::fill::zeros);
            }
            else
                PhiRapN[0] = GetLinSolution(MatEq, Phi0N[0]);
        }
        else { //Dummy start for DGLAP
            for(int y = 0; y <= start; ++y)
                PhiRapN[y] = Phi0N[y];
        }


        for(int y = start+1; y < Nrap; ++y) {
            //Starting point of evol with 0.5 (Trapezius)
            arma::vec yTemp = 0.5 * matN[y] * PhiRapN[0];

            //Remaining without mult
            for(int d = 1; d < y; ++d)
                yTemp += matN[d] * PhiRapN[y-d];

            //Whole right hand side
            yTemp = stepY * yTemp + Phi0N[y];
            
            cout <<"Rap point " << y << endl;
            /*
            auto matNow = mat[0]; //Adding diagonal DGLAP term
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l)
                matNow[k][l] += matDiag[y][k][l];
            */

            arma::mat matEq = arma::mat(N,N,arma::fill::eye) - 0.5*stepY*matN[0] -matNDiag[y];

            //vector<double> Pred = PhiRap[y-1];
            //if(y >= 4)
                //for(int k = 0; k < N; ++k)
                    //Pred[k] += PhiRap[y-1][k] - PhiRap[y-2][k];


            //PhiRap[y] = GetRegSolution(matEq, yTemp, Pred);
            PhiRapN[y] = GetLinSolution(matEq, yTemp);

            bool isGood = true;

            //cout << "Id of y = " << y << endl;
            if(!isGood && 0) {

                for(int j = 0; j < 10; ++j)
                    cout << y <<" " << j <<" "<< PhiRapN[y](j) << endl;

            }
        }

    }

    void Evolve() {
        double stepY = (rapMax - rapMin) / (Nrap-1);

        //PhiRap[0] = Phi0;

        int start = 0;

        //Classical approach
        if(start == 0) {
            auto MatEq = matDiag[0];
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l) {
                MatEq[k][l] *= -1;
                if(k == l) MatEq[k][l] += 1;
                assert(isfinite(MatEq[k][l]));
            }
            if(putZero)
                PhiRap[0].resize(N, 0.);
            else
                PhiRap[0] = GetLinSolution(MatEq, Phi0[0]);


            //if(

            //for(int i = 0; i < N; ++i)
                //cout << i << " "<< PhiRap[0][i] << endl;
            //exit(0);
        }
        else { //Dummy start for DGLAP
            for(int y = 0; y <= start; ++y)
                PhiRap[y] = Phi0[y];
        }


        for(int y = start+1; y < Nrap; ++y) {
            vector<double> yTemp(N,0.);
            //Starting point of evol with 0.5 (Trapezius)
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l)
                yTemp[k] += mat[y][k][l] * PhiRap[0][l];
            for(auto &v : yTemp) v /= 2;

            //Remaining without mult
            for(int d = 1; d < y; ++d) {
                #pragma omp parallel for
                for(int k = 0; k < N; ++k)
                for(int l = 0; l < N; ++l)
                    yTemp[k] += mat[d][k][l] * PhiRap[y-d][l];
            }

            //Whole right hand side
            for(int k = 0; k < N; ++k)
                yTemp[k] = yTemp[k] * stepY + Phi0[y][k];
            
            
            /*
            auto matNow = mat[0]; //Adding diagonal DGLAP term
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l)
                matNow[k][l] += matDiag[y][k][l];
            */

            auto matEq = mat[0];
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l) {
                matEq[k][l] *= -0.5*stepY;
                matEq[k][l] -= matDiag[y][k][l];
                if(k == l) matEq[k][l] += 1;
                assert(isfinite(matEq[k][l]));
            }
            
            for(int k = 0; k < N; ++k)
                assert(isfinite(yTemp[k]));


            //PhiRap[y] = IterSolution(matNow, 0.5 * stepY, yTemp);
            //auto myComp2 = GetLinSolution(matEq, yTemp);

            vector<double> Pred = PhiRap[y-1];

            //if(y >= 4)
                //for(int k = 0; k < N; ++k)
                    //Pred[k] += PhiRap[y-1][k] - PhiRap[y-2][k];


            //PhiRap[y] = GetRegSolution(matEq, yTemp, Pred);
            PhiRap[y] = GetLinSolution(matEq, yTemp);

            bool isGood = true;
            for(int k = 0; k < N; ++k)
                if(!isfinite(PhiRap[y][k])) isGood = false;


            //cout << "Id of y = " << y << endl;
            if(!isGood && 0) {
                cout << "Id of y = " << y << endl;
                
                arma::mat M(N,N);
                for(int i = 0; i < N; ++i)
                for(int j = 0; j < N; ++j)
                    M(i,j) = matEq[i][j];
                //cout <<"Determinant is " <<  arma::det(M) << endl;
                //cout <<"Rank is " <<  arma::rank(M) << endl;

                double s = 0;
                double sr = 0;
                double Min = 1e100;
                double Max =-1e100;
                for(int j = 0; j < N; ++j) {
                    s += 0;//myComp[j];
                    sr += PhiRap[y][j];
                    Min = min(Min, PhiRap[y][j]);
                    Max = max(Max, PhiRap[y][j]);
                }

                //cout <<"Condition is " <<  arma::cond(M) <<" "<< s<<" "<<sr <<":" << Min <<" "<< Max << endl;

                for(int j = 0; j < 10; ++j)
                    cout << y<<" " << j <<" "<< PhiRap[y][j] << endl;

                /*
                for(int j = 0; j < N; ++j)
                    cout << j <<" "<< PhiRap[y][j] << endl;
                cout << "Input vec " << endl;
                for(int j = 0; j < N; ++j)
                    cout << j <<" "<< yTemp[j] << endl;
                */

            }
            //for(int i = 0; i < N; ++i)
                //cout << i << " "<< PhiRap[y][i] <<" "<< test[i] << endl;

            //for(int j = 0; j < 10; ++j)
                //cout <<"Pusa "<< y<<" " << j <<" "<< PhiRap[y][j] << endl;



            //PhiRap[y] = yTemp;
        }

        //for(int j = 0; j < N; ++j)
            //cout << "RADEK "<<" " << j <<" "<< PhiRap[0][j] << endl;
        //exit(0);

    }

    void DoIteration() {
        const double stepY = (rapMax - rapMin) / (Nrap-1);

        vector<vector<double>> PhiRapNew(Nrap);

        PhiRapNew[0] = Phi0[0];
        for(int k = 0; k < N; ++k)
        for(int l = 0; l < N; ++l)
            PhiRapNew[0][k] += matDiag[0][k][l] * PhiRap[0][l];


        for(int y = 1; y < Nrap; ++y) {
            vector<double> yTemp(N,0.); //kT spectrum for particular bin y
            //Starting point of evol with 0.5 (Trapezius)
            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l) {
                yTemp[k] += mat[y][k][l] * PhiRap[0][l];
                yTemp[k] += mat[0][k][l] * PhiRap[y][l];
            }
            for(auto &v : yTemp) v /= 2;

            //Remaining without mult
            for(int d = 1; d < y; ++d) {
                #pragma omp parallel for
                for(int k = 0; k < N; ++k)
                for(int l = 0; l < N; ++l)
                    yTemp[k] += mat[d][k][l] * PhiRap[y-d][l];
            }


            //Whole right hand side
            for(int k = 0; k < N; ++k)
                yTemp[k] = yTemp[k] * stepY + Phi0[y][k];


            for(int k = 0; k < N; ++k)
            for(int l = 0; l < N; ++l)
                yTemp[k] += matDiag[y][k][l] * PhiRap[y][l];



            
            PhiRapNew[y] = yTemp;
        }
        PhiRap = PhiRapNew;
    }
    
    void RunIterations(int Niter, bool init = true) {
        //Init PhiRap 
        if(init) {
            assert(Phi0.size() == PhiRap.size());
            for(int y = 0; y < PhiRap.size(); ++y)
                PhiRap[y] = Phi0[y];
        }

        //Do itrerations
        for(int i = 0; i < Niter; ++i) {
            DoIteration();
            cout << "Iteration " << i <<" done." << endl;
            cout << "Phi[kT0] = " << PhiRap[Nrap-1][0] << endl;
        }
    }


    vector<double> GetRHS(const vector<double> &PHI) {
        vector<double> dPhi(N,0.);
        //#pragma omp parallel for
        for(int i = 0; i < N; ++i) { 
            for(int j = 0; j < N; ++j) {
               dPhi[i] += mat[0][i][j] * PHI[j]; 
            }
        }
        return dPhi;
    }

    void Step(double delta) {
        static int y = 0;
        if(y == 0) 
            PhiRap[0] = Phi0[0];

        vector<double> dPhi = GetRHS(PhiRap[0]);
        ++y;
        for(int i = 0; i < N; ++i) 
            PhiRap[0][i] += delta * dPhi[i];

        //for(int k = 0; k < 30; ++k)
            //cout << " y i " << y <<" "<< k <<" "<< PhiRap[0][k] << endl;
    }


    //Function of x and kT2
    void InitF(function<double(double, double)> fun) {
        const double stepY = (rapMax-rapMin) / (Nrap-1);
        
        //Old version
        Phi0.resize(Nrap);
        for(int y = 0; y < Nrap; ++y) {
            Phi0[y].resize(N);
            double x = exp(-y*stepY);
            for(int i = 0; i < N; ++i) {
                double kT2 = exp(nodBase.xi[i]);
                Phi0[y][i] = fun(x, kT2);
            }
        }
        PhiRap.resize(Nrap);
        for(auto &ph : PhiRap)
            ph.resize(N, 0.);

        //New version
        Phi0N.resize(Nrap);
        for(int y = 0; y < Nrap; ++y) {
            arma::vec temp(Nint);
            double x = exp(-y*stepY);
            for(int i = 0; i < Nint; ++i) {
                double kT2 = exp(nod.xi[i]);
                temp(i) = fun(x, kT2);
            }
            Phi0N[y] = redMat * temp;
        }

        PhiRapN.resize(Nrap, arma::vec(N, arma::fill::zeros));
    }

    double Interpolate(double y, double L, bool isNew=false) {

        y = max(y, rapMin);
        y = min(y, rapMax);

        double stepY = (rapMax-rapMin) / (Nrap-1);
        int yId = (y - rapMin) / stepY;
        int LId = 0;
        for(LId = N-1; LId >= 0; --LId)
            if(nodBase.xi[LId] <= L) break;
        //--LId;
        LId = max(0, LId);
        yId = min(Nrap-2, yId);

        double yLeft = rapMin + yId*stepY;
        double LLeft = nodBase.xi[LId];
        double LRight = nodBase.xi[LId+1];
        

        //if (LRight - L < 0 && L - LRight < 1e-7)
            //L = LRight;

        assert(y - yLeft >= 0);
        assert(yLeft+stepY - yLeft >= 0);
        


        if (LRight - L < 0 || L - LLeft < 0) {
            cout <<"Error " <<setprecision(24)<< L <<" : "<< LLeft <<" "<< LRight << endl;
            //L = LRight;
            //cout << "Ahojky " << nod.xi[0] << " "<< nod.xi[Nrap-1] << endl;
        }
        assert(L - LLeft >= 0);
        assert(LRight - L >= 0);

        if (yId > Nrap - 2) {
            cout << "Ahoj " << y << " :  " << yLeft << " "<< yLeft + stepY << endl;
        }

        assert(yId <= Nrap - 2);


        assert(LId <= N - 2);
        assert(yId >= 0);
        assert(LId >= 0);

        double fLL, fLR, fRL, fRR;

        if(isNew) {
            fLL = PhiRapN[yId](LId);
            fLR = PhiRapN[yId](LId+1);
            fRL = PhiRapN[yId+1](LId);
            fRR = PhiRapN[yId+1](LId+1);
        }
        else {
            fLL = PhiRap[yId][LId];
            fLR = PhiRap[yId][LId+1];
            fRL = PhiRap[yId+1][LId];
            fRR = PhiRap[yId+1][LId+1];
        }

        bool canBeLog = fLL > 0 && fLR > 0 && fRL > 0 && fRR > 0;
            
        if(canBeLog) {
            fLL = log(fLL);
            fLR = log(fLR);
            fRL = log(fRL);
            fRR = log(fRR);
        }

        //cout << exp(fLL) << " "<< exp(fLR) <<" : "<< exp(fRL) << " "<< exp(fRR)<<  endl;

        double sum = 
        (y - yLeft)*(L - LLeft)*fRR  + (yLeft + stepY - y)*(L - LLeft)*fLR +
        (y - yLeft)*(LRight - L)*fRL + (yLeft + stepY - y)*(LRight - L)*fLL;
        sum /= stepY * (LRight - LLeft);
        assert(stepY != 0);
        assert((LRight-LLeft) != 0);

        if(canBeLog)
            return exp(sum);
        else 
            return sum;

    }
    
    void PrintBaseGrid() {
        double stepY = (rapMax - rapMin) / (Nrap-1);
        double stepL = (Lmax - Lmin) / (N-1);


        for(int y = 0; y < Nrap; ++y)
        for(int l = 0; l < N; ++l) {
            double yNow = rapMin + y*stepY;
            double L  = nodBase.xi[l];
            
            double kt2 = exp(L);
            double x = exp(-yNow);

            cout << x << " "<< kt2 <<"  " << PhiRap[y][l] << " "<< PhiRapN[y](l) << endl;
        }
    }

    void PrintGrid() {
        double stepY = (rapMax - rapMin) / (Nrap-1);

        //cout <<"Interpolation "<< Interpolate(0, log(1) ) << endl;
        //return;

        double yStepMich = log(1e-2/1e-8) / 99;//(Nrap -1);
        double LStepMich = log(1e6/1e-2)  / 99;// (N -1);

        cout << "Radek " << endl;
        //for(int y = 0; y < 100; ++y) 
        for(int y = 99; y >= 0; --y) {
            double x = 1e-2 * pow(1e-8/1e-2, y/99.);
            double yNow = y*yStepMich;
            for(int i = 0; i < 100; ++i) {
                double kT2 = 0.01 * pow(1e6/0.01, i/99.);
                double LNow = i*LStepMich + log(0.01);
                if(i == 99) LNow = log(1e6)-1e-9;
                double res = Interpolate(yNow, LNow);
                double resNew = Interpolate(yNow, LNow, true);
                /*if(yNow == 0)*/ cout << x<<" "<< kT2<<" "<< res <<" "<< resNew <<  endl;
            }
        }
    }
};
#endif 
