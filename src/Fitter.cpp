#include "Fitter.h"
#include <cmath>

//#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <gsl/gsl_multimin.h>

#include "TF2.h"
#include "TMinuit.h"


Fitter *glFitter;


vector<dataPoint> Fitter::LoadData(string fname)
{
    vector<dataPoint> points;
    ifstream file(fname);
    if (!file.good())
    {
        cout << "File " << fname << " cant be open." << endl;
        assert(0);
    }


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

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *p, Int_t iflag)
{
        //npar = 3;
        npar = Settings::I().nPar;
        f = glFitter->Eval(p);
}







void Fitter::Init(string dirName)
{
    //Load data points
    data = LoadData("/home/zlebcr/prog/bkEvol/test/heraTables/HERA1+2_NCep_920.dat");

    //Load evoluton and convolution matrices
    //sol512.InitMat();
    solver.LoadEvolKernels(dirName);
    solver.LoadConvKernels(dirName);



    //New game

    glFitter = this;

    TMinuit *gMinuit = new TMinuit(Settings::I().nPar+1);  //initialize TMinuit with a maximum of 5 params
    gMinuit->SetFCN(fcn);
    int ierflg = 0;
    double arglist[10];
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    //static Double_t vstart[2] = {4.62153e+01, 3.69297e+00}; Old variant
    //static Double_t vstart[3] = {357.855, 6.23351, 1.86848};// New variant
    //static Double_t vstart[3] = {-128.786, -184.162, -71.7185};// Tchebyshev fit

    //static Double_t vstart[2] = {46, 3};// new hope
    //static Double_t vstart[2] = { 12.4613, 1.86526}; //  667.186, 14.3659,};

    //New play
    //static Double_t vstart[1] = { 12.4613}; //  667.186, 14.3659,};
//p = {667.269, 14.3668, 0, 0, 0, 0,}; 
//Chi2 is 1080.85 / 180
    //static Double_t vstart[2] = { 3, -3}; //  667.186, 14.3659,};
    static Double_t vstart[3] = {-21.9568, 2.488, 11.6075};





    /*
    //static Double_t vstart[2] = {3, 1 };
    static Double_t step[3] = {0.1 , 0.1, 0.1};
    gMinuit->mnparm(0, "a1", vstart[0], step[0], -40,1000,ierflg);
    gMinuit->mnparm(1, "a2", vstart[1], step[1], -40,1000,ierflg);
    gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);
    */

    auto &S = Settings::I();
    for(int i = 0; i < S.nPar; ++i) {
        gMinuit->mnparm(i, "p"+to_string(i), get<0>(S.pars[i]),  0.1,  get<1>(S.pars[i]),  get<2>(S.pars[i]), ierflg);
    }

    arglist[0] = 1500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    cout << "Success" << endl;
    cout << "Saving results" << endl;

    if(1) {
        //Save result
        arma::field<arma::mat> storage(1, 4);
        storage(0, 0) = Solver::vector2matrix(solver.PhiRapN);
        storage(0, 1) = Solver::vector2matrix(solver.F2rap);
        storage(0, 2) = Solver::vector2matrix(solver.FLrap);
        storage(0, 3) = Fitter::getPoints();

        string nTag = to_string(lrint(1000*solver.asMZ));
        storage.save("fit2Parm_" + nTag + ".dat");
        exit(0);
    }




    return;


    //CalculateBasis(40, "basis.dat");
    //MPI_Finalize();
    //return;

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


                
    cout << "Helenka " << __LINE__ << endl;

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


    vector<double> p;
    //p[0]=3.60325e-06; p[1]=0.0627581;
    //p = {1e-2, 0, 5, 0, 0.1, 0.5};

    //p = {13.3794, -6.7348, 5.05507, 1.19013, 4.70995, 1.39788,};
    //p = {11.8021, -5.97259, 1.8291, 3.60356, 5.13523, 1.24874,};
    //p = {11.7143, -6.24427, 0.0518992, 2.73242, 4.89552, 1.2672,};
    //p = {78.9713, -5.52102, 1e-7, 3.00587, 5.20282, 1.70011,};
    //p = {80.0684, -5.44104, 7.45582e-08, 3.01449, 5.23708, 1.70187,};
    //p = {80.1451, -5.44868, -2.29646e-06, 3.0168, 5.23553, 1.70169,};//for aS = 0.1118

    //p = {80.1505, -5.39695, -5.42062e-06, 2.99661, 5.27126, 1.70571,}; //for Krystof org
    //p = {46.1936, 3.69189}; //for Krystof - simple p0*kT2*exp(-p1*kT2)
    //p = {251.778, 5.77694, 1.71978}; //for Krystof - simple p0*kT2*exp(-p1*kT2)



    //p = {80.1505, +1.39695, -5.42062e-06, 0.39661, 5.27126, 1.70571,}; //for Krystof corr
    p = {80.12, -7.72392, -0.037154, 4.1288, 3.19872, 2.33265,};


    //p = {79.1046, -3.75694, 4.48014e-08, 3.06169, 8.21893, 1.64082,};//for aS = 0.118

    //p = {11.7357, -6.22287, 0.156637, 2.74329, 4.91692, 1.26597,};


    //p[0]=-9.58527, p[1]=7.06043; p[2]= -11.3249;
    //p[0]=-11.2986, p[1]=4.75103; p[2] =  -9.05604; p[3] =  -0.800283;
    //p[0]=-11.2544; p[1]=4.05222; p[2] = -7.91101; p[3] =  -0.798476;

    //p[0]=-2.26559; p[1]=8.30788; p[2] =  -7.52139; p[3] = 0;
    //p[0]=-1.66264; p[1]=8.35095; p[2] =  -9.09875; p[3] =  -0.483668;
    //p[0]=-0.5853; p[1]=5.02503; p[2] = -7.37439; p[3] =  -1.3508;

    //p[0]=0.268914; p[1]=8.01361; p[2]= -26.609; p[3] =  1.87533;

    //p[0]=0.264653; p[1]=-3.19877; p[2] = -8.82977; p[3] = 0.398698;

    //Eval(p.data());
    //return;

    int nPar = p.size();
    x = gsl_vector_alloc (nPar);

    for(int i = 0; i < nPar; ++i)
        gsl_vector_set (x, i, p[i]);

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (nPar);
    //gsl_vector_set_all (ss, 0.2);
    gsl_vector_set_all (ss, 0.7);
    //gsl_vector_set(ss, 2, 0.2);


    /* Initialize method and iterate */
    minex_func.n = nPar;
    minex_func.f = my_f;
    minex_func.params = this;

    s = gsl_multimin_fminimizer_alloc (T, nPar);
    cout << "Helenka start" << __LINE__ << endl;
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    cout << "Helenka end" << __LINE__ << endl;

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) 
            break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-3);



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
    while (status == GSL_CONTINUE && iter < -1000);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);


}

double Fitter::Eval(const double *p)
{
    auto fun = Settings::I().fitFun;
    //cout << "Matrix initialised" << endl;
    solver.InitF([=](double x, double kT2) {
        //return pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        //return 1./pow(kT2,1);// pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);

        //return p[0]*kT2 * exp(-abs(p[1])*kT2);

        //return kT2 * exp(-abs(p[0])*kT2);

        return fun(kT2, x, p);

        //return p[0]*pow(kT2,p[2]) * exp(-p[1]*kT2);// * pow(max(0., 0.4-x), 2);

        //return exp(p[0] + p[1]*log(kT2) + p[2]*pow(log(kT2),2) + p[3]*pow(log(kT2),3) );
        //return p[0] * pow(kT2,p[1]) * pow(1-x,abs(p[2])) * max(0.,1-p[3]*x) * exp(-p[4]*pow(log(kT2/p[5]),2));


        /*
        const double minkT2 = 1e-2;
        const double maxkT2 = 1e6;
        double y = -1 + 2.*(log(kT2) - log(minkT2)) / (log(maxkT2)-log(minkT2));
        double c0 = 1;
        double c1 = y;
        double c2 = 2*y*y - 1;
        */

        //double L = log(kT2);
        //return exp(p[0]*L + p[1]*L*L) * pow(1-x, abs(p[2]));


    });
    solver.EvolveNew();
    AddTheory(solver.F2rap, solver.FLrap);
    int nDF;
    double chi2 = getChi2Corr(nDF);

    //cout << "For parameters p[0]=" << p[0]<<", p[1]="<<p[1] <<" "<< p[2] <<" "<< p[3]<< endl;
    cout << "p = {" << p[0]<<", "<<p[1] <<", "<< p[2] <<", "<< p[3]<< ", "<<p[4]<< ", " << p[5]<<",}; "<< endl;
    cout << "Chi2 is " <<chi2<< " / "<< nDF << endl;


    if(0) {
        //Save result
        arma::field<arma::mat> storage(1, 4);
        storage(0, 0) = Solver::vector2matrix(solver.PhiRapN);
        storage(0, 1) = Solver::vector2matrix(solver.F2rap);
        storage(0, 2) = Solver::vector2matrix(solver.FLrap);
        storage(0, 3) = Fitter::getPoints();

        storage.save("fitTheb.dat");
        exit(0);
    }



    if(isfinite(chi2))
        return chi2;
    else
        return 1e40;

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

    assert(part >= 0);
    assert(part <= 1);
    //linear interpolation
    double val = f[id](Q2id)*(1-part)  + f[id+1](Q2id)*part;

    if(!isfinite(val)) return 0;

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

       p.theor0 = sRed;
       p.theor  = sRed;
    }
}



//Chi2 with correction for the extra term
double Fitter::getChi2Corr(int &nDF)
{
    auto Cond = [](double x, double Q2) {return x < 0.01 && Q2 > 4 && Q2 < 1200; };


    double A11, A12, A22; 
    double r1, r2;
    A11 = A12 = A22 = 0;
    r1 = r2 = 0;

    for(auto &p : data) {
        if(!Cond(p.x, p.Q2)) continue; 

        const double a = -0.08;
        const double b =  8;
        p.extra0 = pow(p.x, a) * pow(1 - p.x, b);

        double Vinv = pow(p.err*p.sigma*1e-2, -2);

        A11 += p.theor0*p.theor0 * Vinv;
        A12 += p.theor0*p.extra0 * Vinv;
        A22 += p.extra0*p.extra0 * Vinv;

        r1 += p.sigma * p.theor0 * Vinv;
        r2 += p.sigma * p.extra0 * Vinv;

    }

    double Disc = A11*A22 - A12*A12;
    double D1 = r1*A22 - r2*A12;
    double D2 = A11*r2 - r1*A12;
    double c1, c2;
    if(Disc != 0) {
        c1 = D1 / Disc;
        c2 = D2 / Disc;
    }
    else {
        c1 = r1 / A11;
        c2 = 0;
        //cout << "Is my case!!!" << endl;
    }

    double chi2 = 0;
    double hardSum = 0;
    double softSum = 0;
    nDF = 0;
    for(auto &p : data) {
        if(!Cond(p.x, p.Q2)) continue; 

        double Vinv = pow(p.err*p.sigma*1e-2, -2);
        p.theor = c1 * p.theor0 + c2 * p.extra0;

        hardSum += c1 * p.theor0;
        softSum += c2 * p.extra0;

        chi2 += pow(p.sigma - p.theor, 2) * Vinv;
        //cout << p.Q2 <<" "<< p.x <<" : "<< p.sigma <<" "<< p.theor << endl;
        ++nDF;
    }

    cout << "Soft fraction " << softSum / (softSum + hardSum) << endl;

    return chi2;

}


double Fitter::getChi2(int &nDF)
{
    double chi2 = 0;
    nDF = 0;
    for(auto &p : data) {
        if(p.x > 0.01 || p.Q2 < 4 || p.Q2 > 1200) continue;

        chi2 += pow((p.sigma - p.theor) / (p.err*p.sigma*1e-2), 2);
       // cout << p.Q2 <<" "<< p.x <<" : "<< p.sigma <<" "<< p.theor << endl;
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
