#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <tuple>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <set>
#include <armadillo>

using namespace std;

inline double pow2(double x) {return x*x;}

vector<double> LoadData(string fname)
{
    ifstream file(fname);
    char s[2000];
    file.getline(s,1990);
    file.getline(s,1990);

    int nPoints = 0, nSmallx=0;

    vector<double> Q2vals;
    while(1) {
        double Q2, x, y, Sigma, stat, uncor;
        double sys[162], delta[4];
        double tot_noproc, delta_rel;
        double delta_gp, delta_had;

        file >> Q2 >>  x >>  y >>  Sigma >>  stat >>  uncor;
        for(int i = 0; i < 162; ++i)
            file >> sys[i];
        file >> tot_noproc >> delta_rel;
        for(int i = 0; i < 4; ++i)
            file >> delta[i];
        file >> delta_gp >> delta_had;
        
        if(!file.good()) break;

        double err2 = pow2(tot_noproc) + pow2(delta_rel) + pow2(delta_gp) + pow2(delta_had);
        for(int i = 0; i < 4; ++i)
            err2 += pow2(delta[i]);

        ++nPoints;
        if(x <= 0.01) {++nSmallx;  }
        Q2vals.push_back(Q2);
        //cout << "Ahoj " << Q2 <<" "<< x <<" "<< Sigma <<" : " << stat <<" "<< sqrt(err2)<< endl;
    }
    cout << nPoints << " "<< nSmallx << endl;
    sort(Q2vals.begin(), Q2vals.end());
    auto last = unique(Q2vals.begin(), Q2vals.end());
    Q2vals.erase(last, Q2vals.end()); 

    for(auto q2 : Q2vals) {
        cout << q2 << endl;
    }
    //exit(0);

    return Q2vals;

}


vector<double> GetWeights(int Size);


struct Kernel {
    double mq2, eq2;
    double Q2, z, p2, p;
    
    void setMassCharge(double m2, double e2) {mq2 = m2; eq2 = e2;}

    void setQ2zP2(double _Q2, double _z, double _p2) {
        Q2 = _Q2;
        z = _z;
        p2 = _p2;
        p = sqrt(p2);
    }
    pair<double,double> Integrand(double k, double cosPhi);
    pair<double,double> SingleIntegrand(double beta, double k, double cosPhi);
    double DerivBeta(double beta, double k, double cosPhi);
    pair<double,double> GetBeta(double z, double k, double cosPhi);
    pair<double,double> MainTerm(double beta, double k, double cosPhi);

};














const double as = 0.2;


const double Lmin= log(1e-2), Lmax = log(1e6);



vector<double> GetNodes(int n, double a, double b)
{
    vector<double> nodes;
    for(int k = 0; k < n; ++k) {
        double Cos = cos(k*M_PI / (n-1.));

        double x = (a+b)/2 - Cos*(b-a)/2;
        nodes.push_back(x);
    }
    return nodes;
}

struct Integrator {
    vector<vector<double>> kNodes, cosPhiNodes;
    vector<vector<double>> kWeights, phiWeights;

    void Init() {
        double precMax = 12;

        kNodes.resize(precMax+1);
        cosPhiNodes.resize(precMax+1);
        kWeights.resize(precMax+1);
        phiWeights.resize(precMax+1);

        for(int i = 3; i <= precMax; ++i) {
            int n = pow(2,i) + 1;

            cosPhiNodes[i]  = GetNodes(n, 0, M_PI);
            for(auto &x : cosPhiNodes[i]) x = cos(x);

            kNodes[i] = GetNodes(n, Lmin/2., Lmax/2.);
            for(auto &x : kNodes[i]) x = exp(x);


            phiWeights[i]= GetWeights(n);
            kWeights[i] = phiWeights[i];

        }
    }

    pair<double,double> GetIntegral(int o1, int o2, double z, double Q2, double p2, double m2, double e2)
    {
        Kernel kern;
        kern.setMassCharge(m2, e2);
        kern.setQ2zP2(Q2, z, p2);


        double FT=0, FL=0;
        for(int i = 0; i < kNodes[o1].size(); ++i)
        for(int j = 0; j < cosPhiNodes[o2].size(); ++j) {
            double k = kNodes[o1][i];
            double cosPhi = cosPhiNodes[o2][j];
            double w = kWeights[o1][i] * phiWeights[o2][j];
            double ft, fl; 
            tie(ft, fl) = kern.Integrand(k, cosPhi);
            FT += ft*w; FL += fl*w;
        }

        double Kphi = 1./2; //average over angle
        double Kk   = (Lmax-Lmin)/2;
        double K = Kphi * Kk;

        return make_pair(K*FT, K*FL);
    }

};


pair<double,double> Kernel::Integrand(double k, double cosPhi)
{
    auto betas = GetBeta(z, k, cosPhi);
    //if( (betas.first < 1 && betas.first > 0) || (betas.second < 1 && betas.second > 0))
        //cout << "beta 12 is " << betas.first<<" "<<betas.second << endl;
    auto res1 = SingleIntegrand(betas.first,  k, cosPhi);
    auto res2 = SingleIntegrand(betas.second, k, cosPhi);

    return make_pair(res1.first+res2.first, res1.second+res2.second);
}

pair<double,double> Kernel::SingleIntegrand(double beta, double k, double cosPhi)
{
    if(beta <= 0 || beta >= 1) return make_pair(0.,0.);

    double der = 1./(z*abs(DerivBeta(beta, k, cosPhi)));
    return MainTerm(beta, k, cosPhi);
}

//Derivative dz^-1/dbeta
double Kernel::DerivBeta(double beta, double k, double cosPhi)
{
    double k2 = k*k;
    double mix = (k*p)*cosPhi;

    double c1 = k2 + mq2;
    double c2 = k2 + mq2 + p2 - 2*mix;


    double r = (c1 / pow2(1-beta) - c2/pow2(beta));

    return r;
}

pair<double,double> Kernel::GetBeta(double z, double k, double cosPhi)
{
    double k2 = k*k;
    double mix = k*p*cosPhi;
    
    double a = 1./z - 1;
    double b = (2*mix - p2)/Q2  - a;
    double c = (k2 - 2*mix + k2 + mq2) / Q2;

    if(z==1) {
        double res = -c/b;
        return make_pair(res, res);
    }
    else {
        double D = b*b - 4*a*c;
        if(D < 0) return make_pair(-1., -1.);

        double sqrtD = sqrt(D);
        return make_pair((-b - sqrtD)/(2*a), (-b + sqrtD)/(2*a));
    }
}

//function of k
pair<double,double> Kernel::MainTerm(double beta, double k, double cosPhi)
{
    double k2 = k*k;
    double mix = k*p*cosPhi;

    double c2  = beta*(1-beta)*Q2 + mq2;
    double kp2 = k2 + p2 - 2*mix;
    double D1  = k2 + c2;
    double D2  = kp2  + c2;
    double f   = beta*beta + (1-beta)*(1-beta);


    double r1 = pow2((D2 - D1) / (D1*D2));
    double r2 = k2/(D1*D1) - 2*(k2-mix)/(D1*D2) + kp2/(D2*D2);
    
    double FT = f * r2 + mq2 * r1;
    double FL = 4*Q2*pow2(beta*(1-beta)) * r1;

    //Overall constant
    double fact = as*eq2*Q2/(4*M_PI) * k2;




    return make_pair(fact*FT, fact*FL);
}



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

void CalculateGrid(string fname, int qid)
{
    const int Nrap = 1024;
    const int N = 512+1;

    int N2id = log2(N-1);
    assert(pow(2,N2id)+1 == N);

    const double Lmin= log(1e-2), Lmax = log(1e6);
    const double rapMax = log(1e6), rapMin = log(1);

    vector<double> Q2data =  LoadData(fname);
    int nQ2 = Q2data.size();
    cout <<"Calculating " << qid <<" from "<< nQ2 << endl;
    assert(qid < nQ2);

    const double rapStep = (rapMax-rapMin) / (Nrap-1.);


    //arma::cube conFT(nQ2,N, Nrap, arma::fill::zeros);
    //arma::cube conFL(nQ2,N, Nrap, arma::fill::zeros);
    arma::mat conF2(N, Nrap, arma::fill::zeros);
    arma::mat conFL(N, Nrap, arma::fill::zeros);

    Integrator quad;
    quad.Init();

    const double mL2 = 0;
    const double mC2 = 1.3*1.3;
    const double eL2 = 2* 1/3.*1/3. + 2/3.*2/3.;
    const double eC2 = 2/3.*2/3.;

    //for(int qid = 5; qid < 6; ++qid) {
        double Q2 = Q2data[qid];
        #pragma omp parallel for
        for(int y = 0; y < Nrap; ++y) {
            double z = exp(-rapMin - rapStep*y);
            for(int i = 0; i < N; ++i) {
                //double L = (Lmin + Lmax)/2 -  cos(M_PI*i/(N-1.))*(Lmax-Lmin)/2;
                //double p2 = exp(L);
                double p = quad.kNodes[N2id][i];
                //cout << i <<" "<< p << endl;

                auto resL =  quad.GetIntegral(9, 7, z, Q2, p*p, mL2, eL2);
                auto resC =  quad.GetIntegral(9, 7, z, Q2, p*p, mC2, eC2);

                conF2( i, y) = (resL.first + resC.first) + (resL.second + resC.second);
                conFL( i, y) = resL.second + resC.second;
                //cout << fixed<<setprecision(5) << res.first << " " ;

                //Applying weights for intregration over dz/z and dp2/p2
                conF2(i, y) *= rapStep * quad.kWeights[N2id][i] * (Lmax-Lmin)/2.;
                conFL(i, y) *= rapStep * quad.kWeights[N2id][i] * (Lmax-Lmin)/2.;
            }
            //cout <<"Z is " <<  z << endl;
        }
    //}

    //string address = "/nfs/dust/cms/user/zlebcr/Krakow/convMat/";
    string address = "";

    conF2.save(address+"convF2_"+to_string(qid)+".h5", arma::hdf5_binary);
    conFL.save(address+"convFL_"+to_string(qid)+".h5", arma::hdf5_binary);
}







int main(int argc, char **argv)
{
    //LoadData("heraTables/HERA1+2_NCep_920.dat");
    //LoadData("heraTables/HERA1+2_NCem.dat");
    //return 0;
    assert(argc == 2);
    int qid = stoi(argv[1]);
    CalculateGrid("/afs/desy.de/user/z/zlebcr/h1/TMD/Krakow/bkEvol/test/heraTables/HERA1+2_NCep_920.dat", qid);

    return 0;

    Integrator quad;
    quad.Init();

    double z =0.5, Q2 =5, p2 =1,  mq2 = 1;

    for(int i = 3; i <= 12; ++i) {
        for(int j = 3; j <= 12; ++j) {
            auto res =  quad.GetIntegral(i, j, z, Q2, p2, mq2, 1./9);
            cout << fixed<<setprecision(5) << res.first << " " ;
        }
        cout << endl;
    }

    return 0;
}
