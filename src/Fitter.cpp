#include "Fitter.h"
#include <cmath>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


vector<dataPoint> Fitter::LoadData(string fname)
{
    vector<dataPoint> points;
    ifstream file(fname);
    char s[2000];
    file.getline(s,1990);
    file.getline(s,1990);

    int nPoints = 0, nSmallx=0;

    map<double, int> q2vals;
    while(1) {
        dataPoint p;
        double Q2, x, y, Sigma, stat, uncor;
        double sys[162], delta[4];
        double tot_noproc, delta_rel;
        double delta_gp, delta_had;

        file >> p.Q2 >>  p.x >>  p.y >>  p.sigma >>  stat >>  uncor;
        for(int i = 0; i < 162; ++i)
            file >> sys[i];
        file >> tot_noproc >> delta_rel;
        for(int i = 0; i < 4; ++i)
            file >> delta[i];
        file >> delta_gp >> delta_had;
        
        if(!file.good()) break;

        double err2 = pow(tot_noproc,2) + pow(delta_rel,2) + pow(delta_gp,2) + pow(delta_had,2);
        for(int i = 0; i < 4; ++i)
            err2 += pow(delta[i],2);

        p.err = sqrt(err2); //in %!
        points.push_back(p);
        ++q2vals[p.Q2];
    }

    int i = 0;
    for(auto &x : q2vals) {
        x.second = i;
        ++i;
    }

    for(auto &p : points)
        p.q2ID = q2vals[p.Q2];

    return points;

}
void Fitter::Init()
{
    //Load data points
    data = LoadData("test/heraTables/HERA1+2_NCep_920.dat");

    //Load evoluton and convolution matrices
    //sol512.InitMat();
    solver.LoadEvolKernels("data/kernel");

    solver.LoadConvKernels("data/kernel");

    //MPI_Finalize();
    //return 0;


    /*
    ROOT::Math::Minimizer* minimiser =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    assert(minimiser);
    minimiser->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    minimiser->SetMaxIterations(10000);  // for GSL
    minimiser->SetTolerance(0.001);
    minimiser->SetPrintLevel(1);

    const int nPar = 2;
    ROOT::Math::Functor funct(*this, nPar);
    minimiser->SetFunction(funct);

    minimiser->SetVariable(0, "normalization",1,0.5);
    minimiser->SetVariable(1, "slope",1,0.5);
    minimiser->SetVariableLimits(0,0.0,5);
    minimiser->SetVariableLimits(1,0.0,10);
    */

    double p[] = {0.001,2};
    (*this)(p);
    //p[1] = 2;
    //(*this)(p);
    //p[1] = 1;
    //(*this)(p);
}


double Fitter::operator()(const double *p)
{
    cout << "Matrix initialised" << endl;

    solver.InitF([=](double x, double kT2) {
        //return pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        //return 1./pow(kT2,1);// pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        return p[0]*kT2 * exp(-p[1]*kT2);// * pow(max(0., 0.4-x), 2);
    });
    solver.EvolveNew();
    AddTheory(solver.F2rap, solver.FLrap);
    int nDF;
    double chi2 = getChi2(nDF);

    cout << "Chi2 is " <<chi2<< " / "<< nDF << endl;

    return chi2;
}


void Fitter::DoFit()
{
    //Do fitting



}

//Value in 
double getValue(vector<arma::vec> &f, int Q2id, double x)
{
    const double rapMin = log(1), rapMax = -log(1e-6);
    int nRap = f.size();
    
    double rap = -log(x);

    double pos = (rap-rapMin)/(rapMax-rapMin) * (nRap-1);
    int id = int(pos);
    double part = pos - id;

    //linear interpolation
    double val = f[id](Q2id)*(1-pos)  + f[id+1](Q2id)*pos;
    
    return val;
}

void Fitter::AddTheory(vector<arma::vec> &F2, vector<arma::vec> &FL)
{
    for(auto &p : data) {
       double f2 = getValue(F2, p.q2ID, p.x);
       double fl = getValue(FL, p.q2ID, p.x);

       double yTerm = p.y*p.y/(1+(1-p.y)*(1-p.y));
       double sRed = f2 - yTerm * fl;

       p.theor = sRed;
    }
}

double Fitter::getChi2(int &nDF)
{
    double chi2 = 0;
    nDF = 0;
    for(auto &p : data) {
        if(p.x > 0.01 || p.Q2 < 4) continue;

        chi2 += pow((p.sigma - p.theor) / (p.err*p.sigma*1e-2), 2);
        cout << p.Q2 <<" "<< p.x <<" : "<< p.sigma <<" "<< p.theor << endl;
        ++nDF;
    }
    return chi2;
}
