#include "Fitter.h"
#include <cmath>

//#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <gsl/gsl_multimin.h>

#include "TF2.h"

//Fitter *fitter;


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


    double
my_f (const gsl_vector *v, void *params)
{
    double x[10];
    Fitter *fit = (Fitter *)params;
    //cout <<"RADEKsize " << v->size << endl;
    for(int i = 0; i < v->size; ++i)
        x[i] = gsl_vector_get(v, i);

    //return pow(x-0.3, 2) + pow(y-0.5,2);
    double *q = nullptr;
    return fit->Eval(x);

    //return p[2] * (x - p[0]) * (x - p[0]) +
        //p[3] * (y - p[1]) * (y - p[1]) + p[4]; 
}






void Fitter::Init()
{
    //Load data points
    data = LoadData("test/heraTables/HERA1+2_NCep_920.dat");

    //Load evoluton and convolution matrices
    //sol512.InitMat();
    solver.LoadEvolKernels("data/kernel");

    solver.LoadConvKernels("data/kernel");


    CalculateBasis(40, "basis.dat");

    MPI_Finalize();
    return;

    //fitter = this;

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


                

    //TF2 *fun = new TF2("function", *this, 0.0, 1.0, 0.001, 10);
    //TF2 * fun = new TF2("fun",[&](double*x, double *q){int n; return getChi2(n); }, 0, 1, 0.0001, 10);
    /*
    TF2 * fun = new TF2("fun",this, &Fitter::Eval, 1e-10, 0.6001, 1e-10, 10, 0);

    //double p1, p2;
    //fun->GetMinimum(p1, p2);
    double p[2] = {0.1, 1};
    //x
    //
    p[0]=4.66136e-08; p[1]=0.00211218;
    //p[0]=0.000302243;
    //p[1]=1.79951;
    fun->GetMinimum(p);
    */

    //double p[] = {0.001,2};
    //(*this)(p);
    //p[1] = 2;
    //(*this)(p);
    //p[1] = 1;
    //(*this)(p);


    const gsl_multimin_fminimizer_type *T = 
        gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    /* Starting point */

    int nPar = 4;
    x = gsl_vector_alloc (nPar);

    vector<double> p;
    //p[0]=3.60325e-06; p[1]=0.0627581;
    p = {1e-4, 0, 0, 0};
    //p[0]=-9.58527, p[1]=7.06043; p[2]= -11.3249;
    //p[0]=-11.2986, p[1]=4.75103; p[2] =  -9.05604; p[3] =  -0.800283;
    p[0]=-11.2544; p[1]=4.05222; p[2] = -7.91101; p[3] =  -0.798476;

    Eval(p.data());
    return;

    for(int i = 0; i < nPar; ++i)
        gsl_vector_set (x, i, p[i]);

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (nPar);
    gsl_vector_set_all (ss, 10);


    /* Initialize method and iterate */
    minex_func.n = nPar;
    minex_func.f = my_f;
    minex_func.params = this;

    s = gsl_multimin_fminimizer_alloc (T, nPar);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) 
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-4);



        if (status == GSL_SUCCESS)
        {
            printf ("converged to minimum at\n");
        }

        printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
                iter,
                gsl_vector_get (s->x, 0), 
                gsl_vector_get (s->x, 1), 
                s->fval, size);
    }
    while (status == GSL_CONTINUE && iter < 1);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);


}

double Fitter::Eval(const double *p)
{
    cout << "Matrix initialised" << endl;
    solver.InitF([=](double x, double kT2) {
        //return pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        //return 1./pow(kT2,1);// pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        //return p[0]*kT2 * exp(-p[1]*kT2);// * pow(max(0., 0.4-x), 2);
        return exp(p[0] + p[1]*log(kT2) + p[2]*pow(log(kT2),2)  + p[3]*log(x)) ;
    });
    solver.EvolveNew();
    AddTheory(solver.F2rap, solver.FLrap);
    int nDF;
    double chi2 = getChi2(nDF);

    cout << "For parameters p[0]=" << p[0]<<", p[1]="<<p[1] <<" "<< p[2] <<" "<< p[3]<< endl;
    cout << "Chi2 is " <<chi2<< " / "<< nDF << endl;

    return chi2;

}

double Fitter::operator()(const double *p, const double *q)
{
}


void Fitter::CalculateBasis(int nElem, string name)
{
    arma::field<arma::mat> storage(nElem, 4);
    for(int iEl = 0; iEl < nElem; ++iEl) {
        //Do fitting
        cout << "Matrix initialised" << endl;
        solver.InitF([&](double x, double kT2) {
            //return pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
            //return 1./pow(kT2,1);// pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
            //return p[0]*kT2 * exp(-p[1]*kT2);// * pow(max(0., 0.4-x), 2);

            double t = 2.*(log(kT2)- solver.Lmin) / (solver.Lmax - solver.Lmin) - 1;
            return cos(iEl*acos(t));

        });
        solver.EvolveNew();
        AddTheory(solver.F2rap, solver.FLrap);

        storage(iEl, 0) = Solver::vector2matrix(solver.PhiRapN);
        storage(iEl, 1) = Solver::vector2matrix(solver.F2rap);
        storage(iEl, 2) = Solver::vector2matrix(solver.FLrap);
        storage(iEl, 3) = Fitter::getPoints();


        //int nDF;
        //double chi2 = getChi2(nDF);
        //cout << "Chi2 is " <<chi2<< " / "<< nDF << endl;

    }
    storage.save(name);

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
    double val = f[id](Q2id)*(1-part)  + f[id+1](Q2id)*part;
    assert(val >= min(f[id](Q2id), f[id+1](Q2id)));
    assert(val <= max(f[id](Q2id), f[id+1](Q2id)));
    
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

arma::mat Fitter::getPoints()
{
    arma::mat dataMat(data.size(), 5);
    for(int i = 0; i < data.size(); ++i) {
       dataMat(i, 0) = data[i].x; 
       dataMat(i, 1) = data[i].Q2; 
       dataMat(i, 2) = data[i].sigma; 
       dataMat(i, 3) = data[i].err; 
       dataMat(i, 4) = data[i].theor; 
    }
    return dataMat;
}
