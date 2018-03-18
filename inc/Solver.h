#ifndef Solver_H_
#define Solver_H_

#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <string>

#include <armadillo>
#include <mpi.h>

#include "gpuBooster.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "alphaSpline.h"

#include "Settings.h"

using namespace std;

vector<double> GetWeights(int Size);
pair<arma::mat, arma::mat> GetTransMatrices(int base, int ext, bool toTrivial);


inline pair<int,int> GetRankSize()
{
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return make_pair(world_rank, world_size);
}


inline pair<long long,long long> GetStartEnd(long long Min, long long Max)
{
    int rank, nrank;
    tie(rank,nrank) = GetRankSize();
    return GetStartEnd(nrank, rank, Min, Max);
}



struct Nodes {
    Nodes(int _N, double _a, double _b) : N(_N), a(_a), b(_b) {}
    Nodes() {}
    void Init(int _N, double _a, double _b) { N=_N; a=_a; b=_b; }
    int N;
    double a, b;
    vector<double> xi, SqrtExpXi,  wi;
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


        SqrtExpXi.resize(xi.size());
        for(int i = 0; i < xi.size(); ++i)
            SqrtExpXi[i] = exp(0.5*xi[i]);

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
    double asMZ = 0.2;
    double LnFreeze2 = 2*log(1);
    double eps = 1e-7;
    int Nint; // kT nodes in Nintegral
    int N;// = 32*16 + 1; //must be 2*n+1
    int Nrap = 1024;
    bool toTrivial = true;

    double Lmin= log(1e-2), Lmax = log(1e6);
    //const double Lmin= log(1e-4), Lmax = log(1e8);
    double mu2 = 1e-2;
    double rapMax = log(1e6), rapMin = log(1);
    bool putZero = true;

    string inputDir, outputDir;

    Nodes nod, nodBase;

    double alphaS(double l, double lp);
    
    Solver(int N_) : Nint(1*(N_-1)+1), N(N_),
        nod(Nint, Lmin, Lmax), nodBase(N, Lmin, Lmax) {
        nod.CalcNodes(false);
        nodBase.CalcNodes(false);
        tie(redMat,extMat) = GetTransMatrices(N, Nint, toTrivial);
    }


    void Init(std::istream &Stream) {

        /*
        boost::property_tree::ptree tree;
        boost::property_tree::ini_parser::read_ini(Stream, tree);

        bool bkSolverGrid;
        try {

            asMZ = tree.get<double>("Constants.alphaS");
            LnFreeze2 = 2*log(tree.get<double>("Constants.freezingScale"));

            eps = tree.get<double>("Constants.eps");
            mu2 = tree.get<double>("Constants.mu2");
            //Rapidity properties
            Nrap = tree.get<int>("RapiditySpace.Nrap");
            rapMax = -log( tree.get<double>("RapiditySpace.xMin") );
            rapMin = -log( tree.get<double>("RapiditySpace.xMax") );


            //Transverse properties
            N = tree.get<int>("TransverseSpace.NkT2");
            Nint = tree.get<int>("TransverseSpace.NkT2int");
            Lmin = log( tree.get<double>("TransverseSpace.kT2Min") );
            Lmax = log( tree.get<double>("TransverseSpace.kT2Max") );
            bkSolverGrid = tree.get<bool>("TransverseSpace.bkSolverGrid");
            toTrivial = tree.get<bool>("TransverseSpace.toTrivial");

            inputDir  = tree.get<string>("Files.inputDir");
            outputDir = tree.get<string>("Files.outputDir");
            
            //Fit Properties
            string funStr = tree.get<string>("Fit.function");
            bool isDone = false;
            int nPar = 0;
            for(int i = 0; i < 9; ++i)
                if(funStr.find("p["+to_string(i)+"]") != string::npos) {
                    ++nPar;
                    if(isDone) {
                        cout << "There is a gap between parameters" << endl;
                        assert(0);
                    }
                }
                else {
                    isDone = true;
                }
            for(int i = 0; i < nPar; ++i) {
                string par = tree.get<string>("Fit.p"+to_string(i));
                istringstream iss(par);
                vector<double> pars;
                double pNow;
                while(iss >> pNow) pars.push_back(pNow);
                assert(pars.size() == 1 || pars.size() == 3);

                double p, pmin, pmax;
                p = pars[0];
                pmin = pmax = 0;
                if(pars.size() == 3) {
                    pmin = pars[1];
                    pmax = pars[2];
                }

                cout << "Reading parameter "<< i <<" : " << p << " "<< pmin << " "<< pmax << endl;
                if(iss.good()) cout << "String is good " << p <<" "<< pmin <<" "<< pmax << endl;
            }
            exit(0);

        }
        catch(const std::exception& e) {
            cout << "Some of parameters in steering not defined:" << endl;
            cout << e.what() << endl;
            exit(1);
        }

        assert((N - 1) %2 == 0);
        assert((Nint - 1) %2 == 0);


        cout << "Used evolution parameters" << endl;
        cout << "Rapidity Space" << endl;
        cout << "xMin = " << exp(-rapMax) << endl;
        cout << "xMax = " << exp(-rapMin) << endl;

        cout << "Nrap = " << Nrap << endl;
        cout << "N = " << N << endl;
        cout << "bkSolverGrid = " << bkSolverGrid << endl;

        //exit(1);

        if(bkSolverGrid) assert(N == Nint);
        */
        Settings::I().Init(Stream);

        auto &S = Settings::I();
        asMZ = S.asMZ;
        LnFreeze2 = S.LnFreeze2;
        eps = S.eps;
        Nint = S.Nint; // kT nodes in Nintegral
        N = S.N;// = 32*16 + 1; //must be 2*n+1
        Nrap = S.Nrap;
        bool bkSolverGrid = S.bkSolverGrid;
        toTrivial = S.toTrivial;

        Lmin = S.Lmin;
        Lmax = S.Lmax;
        mu2 = S.mu2;
        rapMax = S.rapMax;
        rapMin = S.rapMin;

        putZero = S.putZero;

        inputDir = S.inputDir;
        outputDir= S.outputDir;



        nod.Init(Nint, Lmin, Lmax);
        nodBase.Init(N, Lmin, Lmax);
        nod.CalcNodes(bkSolverGrid);
        nodBase.CalcNodes(bkSolverGrid);
        tie(redMat,extMat) = GetTransMatrices(N, Nint, toTrivial);

        alphaSpline::FixMasses( 1e-8, 4.2,	1e21);
        alphaSpline::FixParameters(2, asMZ, 5, 91.2);
        //cout << "Radek " << alphaSpline::alphaS(2*log(91.2))<< endl;

    }
    Solver(std::istream &Stream)  { Init(Stream);  }
    //Solver()  {}



    //vector<arma::mat> matN, matNDiag;
    arma::cube matN, matNDiag, matNInv;
    vector<arma::vec> PhiRapN, Phi0N;

    arma::cube convF2, convFL;
    vector<arma::vec> F2rap, FLrap;
    
    //vector<arma::vec> &GetF2() {return F2rap;}
    //vector<arma::vec> &GetFL() {return FLrap;}


    arma::mat extMat, redMat;

    gpuBooster gpu;

    pair<double,double> GetKerPar(double l, double lp);
    double Delta(double z, double k2, double q2);

    double KernelBFKLDiag(double l, double lp, double z);
    double KernelBFKL(double l, double lp, double z);
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

    double DGLAPtermSimp(double l, double lp, double z);
    double Kernel9(double l, double lp, double z);
    double Kernel9Diag(double l, double lp, double z);
    double Kernel9zDiag(double l, double lp, double z);


    void InitMat();
    void SetSolution(function<double(double, double)> fun);

    void SaveEvolKernels(string file) {
        //string aStag = to_string(1000*asMZ);
        //file += "_as"+ to_string(lrint(1000*asMZ));
        cout << "Saving Evol Kernels to " << file << endl;
        matN.save(file+"/kernel_base.h5", arma::hdf5_binary);
        matNDiag.save(file+"/kernel_diag.h5", arma::hdf5_binary);
        matNInv.save(file+"/kernel_inv.h5", arma::hdf5_binary);
    }
    void LoadEvolKernels(string file) {
        matN.load(file+"/kernel_base.h5", arma::hdf5_binary);
        matNDiag.load(file+"/kernel_diag.h5", arma::hdf5_binary);
        cout << "RADEK size " << matNDiag.slice(0).n_rows << endl;
        matNInv.load(file+"/kernel_inv.h5", arma::hdf5_binary);
    }

    void LoadConvKernels(string file) {
        assert(convF2.load(file+"/conv_F2.h5", arma::hdf5_binary));
        assert(convFL.load(file+"/conv_FL.h5", arma::hdf5_binary));

        //MPI_Bcast(convF2.memptr(), convF2.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(convFL.memptr(), convF2.n_elem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }




    arma::vec GetLinSolution(const arma::mat &Mat, const arma::vec &y) {
       return arma::solve(Mat, y); 
    }

    vector<double> GetRegSolution(const vector<vector<double>> &MatEq, const vector<double> &y, const vector<double> &yReg);
    arma::vec IterSolution(const arma::mat &Mat, double factor, const arma::vec &y);


    void EvolveNew();
    void CalcF2L();
    void DoIteration();
    void RunIterations(int Niter, bool init = true);
    arma::vec GetRHS(const arma::vec &PHI);
    void Step(double delta);


    //Function of x and kT2
    void InitF(function<double(double, double)> fun);

    double Interpolate(double y, double L);
    void PrintBaseGrid();
    void PrintGrid();
    void PrintReduce();

    static arma::mat vector2matrix(vector<arma::vec> &Vec) {
        arma::mat matPhi(Vec.size(), Vec[0].n_rows);
        for(int i = 0; i < Vec.size(); ++i)
            matPhi.row(i) =  Vec[i].t();
        return matPhi;
    }

};
#endif 
