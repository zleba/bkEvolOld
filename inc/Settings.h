#ifndef _Settings_
#define _Settings_

#include <string>
#include <cassert>
#include <dlfcn.h>


using namespace std;

struct Settings {

    double asMZ = 0.2;
    double LnFreeze2 = 2*log(1);
    double eps = 1e-7;
    int Nint; // kT nodes in Nintegral
    int N;// = 32*16 + 1; //must be 2*n+1
    int Nrap = 1024;
    bool bkSolverGrid = false;
    bool toTrivial = true;

    double Lmin= log(1e-2), Lmax = log(1e6);
    //const double Lmin= log(1e-4), Lmax = log(1e8);
    double mu2 = 1e-2;
    double rapMax = log(1e6), rapMin = log(1);
    bool putZero = true;

    string inputDir, outputDir;

    string funStr;
    int maxIter;
    vector<tuple<double,double,double>> pars;
    int nPar;
    double (*fitFun)(double kT2, double x, const double *p);

    static Settings & I() {
        static Settings instance;
        return instance;
    }

    void Init(istream &Stream) {

        boost::property_tree::ptree tree;
        boost::property_tree::ini_parser::read_ini(Stream, tree);


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
            maxIter = tree.get<int>("Fit.maxIter");
            funStr = tree.get<string>("Fit.function");
            bool isDone = false;
            nPar = 0;
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
                vector<double> parsNow;
                double pNow;
                while(iss >> pNow) parsNow.push_back(pNow);
                assert(parsNow.size() == 1 || parsNow.size() == 3);

                double p, pmin, pmax;
                p = parsNow[0];
                pmin = pmax = 0;
                if(parsNow.size() == 3) {
                    pmin = parsNow[1];
                    pmax = parsNow[2];
                }
                pars.push_back(make_tuple(p, pmin, pmax));

                cout << "Reading parameter "<< i <<" : " << p << " "<< pmin << " "<< pmax << endl;
                //if(iss.good()) cout << "String is good " << p <<" "<< pmin <<" "<< pmax << endl;
            }
            //exit(0);

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

        cout << "Compiling function: " << endl;
        cout <<  funStr << endl;
        Compile(funStr);

        //cout << "f(1,0.01) = " << fitFun(1, 0.01) << endl;
        //cout << "Done " << endl;
        //exit(0);

    }

    void WriteFunction (string str)
    {
        ofstream f("fitFun.cpp");
        f << "#include <math.h>" << endl;
        f << "extern \"C\" {" << endl;
        f << "double fitFun(double kT2, double x, const double *p) {" << endl;
        f << "    return (" << str <<");"<< endl;
        f << "}" << endl;
        f << "}" << endl;
        //return 
    }
                                    

    void Compile (string str)
    {
        //double (*fun)(double x);

        string f = __FILE__;
        f.substr(0, f.size() -  10);//inc/Settings.h 

        #define STRINGIFY(x) #x
        #define TOSTRING(x) STRINGIFY(x)

        string n = TOSTRING(pwdDir);
        n.substr(0, f.size() - 14);
        n += "/obj/fitFun.so";
        cout <<"RADEK " <<  f <<" "<< n << endl;


        
        WriteFunction(str);
        system(("g++ fitFun.cpp -o " + n + " -shared -fPIC").c_str());
        //system("sleep 5");
        void *lib = dlopen(n.c_str(), RTLD_LAZY);
        assert(lib);

        fitFun = reinterpret_cast<double (*)(double,double,const double*)>(dlsym(lib, "fitFun"));
        //cout << dlerror() << endl;
        assert(dlerror() == NULL);
        assert(fitFun);

        //cout << "Fun " <<  fitFun(5, 0.5) << endl;

        return;
    }
                        

};

#endif
