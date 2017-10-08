#include <armadillo>
#include <string>
#include <cassert>
#include <iomanip>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TPad.h"

#include <functional>

using namespace std;

double Max(TGraph *gr) {
    double m = -1e20;
    for(int i = 0; i < gr->GetN(); ++i) {
        double x, y;
        gr->GetPoint(i, x, y);
        m = max(m, y);
    }
    return m;
}


struct LinearFitter {
    arma::mat A, Vinv;
    arma::vec xVec, yData;
    arma::field<arma::mat> fiel;

    function<bool(double,double)> Selector =[](double x, double Q2)
        {return x < 0.01 && Q2 > 4 && Q2 < 400;};

    void Load(string file, int nPolTrial) {
        fiel.load(file);

        int Npols = min<int>(fiel.n_rows, nPolTrial);

        A.zeros(500, Npols);
        yData.zeros(500);
        Vinv.zeros(500, 500);

        int nPoints = 0;



        for(int i = 0; i < Npols; ++i) {
            int iData = 0;
            for(int d = 0; d < fiel(i,3).n_rows; ++d) {
                double x  = fiel(i,3)(d,0);
                double Q2 = fiel(i,3)(d,1);
                double sig= fiel(i,3)(d,2);
                double err= fiel(i,3)(d,3);
                double theo= fiel(i,3)(d,4);

                if(!Selector(x,Q2)) continue;

                A(iData, i) = theo;

                yData(iData) = sig;
                double errAbs = err*sig * 1e-2;
                Vinv(iData, iData) = 1./errAbs/errAbs;

                ++iData;
            }
            nPoints = iData;
        }
        A.resize(nPoints, Npols);
        yData.resize(nPoints);
        Vinv.resize(nPoints, nPoints);

        //cout << "Migration matrix " << endl;
        //cout << yData << endl;
        //cout << Vinv << endl;

    }

    void getMinimum() {
        arma::mat E = A.t() * Vinv * A;
        cout << "Condition " << cond(E) << endl;

        assert(E.n_rows == E.n_cols);
        assert(E.n_rows == E.n_cols);
        xVec = inv(E) *  A.t() * Vinv * yData;

        //cout <<setprecision(15)<< xVec << endl;
        //for(int i = 0; i < xVec.n_rows; ++i)
            //cout << xVec(i) << endl;
    }

    double getChi2() {
        //xVec(xVec.n_rows -1) = 0;


           /*
            xVec(0) = 0.0763032345566899;
            xVec(1) =-0.212280978448689;
            xVec(2) = 0.155851881485432;
            xVec(3) =-0.217567288316786;
            xVec(4) = 0.150753792375326;
            xVec(5) =-0.205685093998909;
            xVec(6) = 0.112218377180398;
            xVec(7) =-0.216600357554853;
            xVec(8) = 0.184430223889649;
            xVec(9) =-0.261057717725635;
            xVec(10) = 0.0819770824164152;
            xVec(11) =-0.102573012933135;
            xVec(12) = 0.0981829538941383;
            xVec(13) =-0.219560103490949;
            xVec(14) = 0.145735025871545;
            xVec(15) =-0.0764536219649017;
            xVec(16) =0;
           */



        arma::vec yTheor = A*xVec;

        double sum = 0;
        for(int i = 0; i < yData.n_rows; ++i) {
            double d = yData(i);
            double t = yTheor(i);
            //cout <<"HU "<< i <<" "<< d <<" "<< t << " "<< (d-t)*sqrt(Vinv(i,i)) << endl;
            sum += pow(d-t,2)*Vinv(i,i);
        }

        arma::mat res = (yTheor - yData).t() * Vinv * (yTheor - yData);

        //cout << "RADEK size " << res.n_rows << " "<< res.n_cols << endl;
        cout << "Helenka " << sum <<" "<< res(0,0) << endl;
        return res(0,0);
    }
    int getNdf() const { return yData.n_rows; }



    void PrintReduce() {
        
        const double rapMax = log(1e6), rapMin = log(1);
        int Nrap = fiel(0,1).n_rows;

        const double stepY = (rapMax - rapMin) / (Nrap-1);

        vector<double> q2Arr={0.15, 0.2, 0.25, 0.35, 0.4, 0.5, 0.65, 0.85, 1.2, 1.5, 2, 2.7, 3.5, 4.5,
        6.5, 8.5, 10, 12, 15, 18, 22, 27, 35, 45, 60, 70, 90, 120, 150, 200, 250,
        300, 400, 500, 650, 800, 1000, 1200, 1500, 2000, 3000, 5000, 8000, 12000, 20000, 30000};

        double sBeam = pow(318.12,2);

        //assert(F2rap.size() == matN.n_slices);
        //assert(FLrap.size() == matN.n_slices);

        for(int rapID = 0; rapID < Nrap; ++rapID) {
            double rap = rapMin + rapID*stepY;
            double x = exp(-rap);

            for(int i = 0; i < 46; ++i) {
                double Q2 = q2Arr[i];
                double y = Q2/(x*sBeam);
                double yPlus = 1+(1-y)*(1-y);

                //cout << " rapID " << rapID << " " << F2rap[rapID].n_rows << endl;
                //assert(F2rap[rapID].n_rows == 46);
                //assert(FLrap[rapID].n_rows == 46);

                double f2=0, fl=0;
                //cout << xVec.n_rows <<" "<< A.n_rows<< endl;
                //cout << fiel(0,1).n_rows <<" "<<  fiel(0,1).n_cols  << endl;
                for(int k = 0; k < A.n_cols; ++k) {
                    f2 += fiel(k, 1)(rapID, i) * xVec(k);
                    fl += fiel(k, 2)(rapID, i) * xVec(k);
                }

                double sRed = f2 - y*y/yPlus * fl;

                cout << x << " "<< Q2 <<" "<< sRed << endl;
                //cout << x << " "<< Q2 <<" "<< F2rap[rapID](i) <<" "<< FLrap[rapID](i) << endl;

            }
        }
    }

    void PrintDistrib(TString fname) {
        const double rapMax = log(1e6), rapMin = log(1);
        int Nrap = fiel(0,1).n_rows;

        const double stepY = (rapMax - rapMin) / (Nrap-1);

        for(int k = 0; k < fiel(0,0).n_rows && k < 1; ++k) {

            double rap = rapMin + k*stepY;
            double x = exp(-rap);

            for(int l = 0; l < fiel(0,0).n_cols; ++l) {
                double pdf = 0;
                for(int i = 0; i < A.n_cols; ++i) {
                    pdf += fiel(i,0)(500, l) * xVec(i);
                }
                cout << l << " "<< pdf << endl;
            }
        }
    }

    void PrintPDF(int idQ2, TString fname)
    {
        const double rapMax = log(1e6), rapMin = log(1);
        int Nrap = fiel(0,1).n_rows;

        const double stepY = (rapMax - rapMin) / (Nrap-1);

        vector<double> q2Arr={0.15, 0.2, 0.25, 0.35, 0.4, 0.5, 0.65, 0.85, 1.2, 1.5, 2, 2.7, 3.5, 4.5,
        6.5, 8.5, 10, 12, 15, 18, 22, 27, 35, 45, 60, 70, 90, 120, 150, 200, 250,
        300, 400, 500, 650, 800, 1000, 1200, 1500, 2000, 3000, 5000, 8000, 12000, 20000, 30000};

        double sBeam = pow(318.12,2);

        double Q2 = q2Arr[idQ2];



        TGraph *grTheor = new TGraph();
        TGraph *grF2 = new TGraph();

        for(int k = 0; k < fiel(0,1).n_rows; ++k) {

            double rap = rapMin + k*stepY;
            double x = exp(-rap);

            double y = Q2/(x*sBeam);
            double yPlus = 1+(1-y)*(1-y);


            double f2 = 0, fl = 0;

            for(int i = 0; i < A.n_cols; ++i) {
                f2 += fiel(i,1)(k, idQ2) * xVec(i);
                fl += fiel(i,2)(k, idQ2) * xVec(i);
            }

            
            double sRed = (y<1) ?  f2 - y*y/yPlus * fl : 0;

            grTheor->SetPoint(k, x, sRed);
            grF2->SetPoint(k, x, f2);
            //cout <<Q2<<" "<< x <<" "<< sRed << endl;
        }

        TGraphErrors *grData      = new TGraphErrors();
        TGraphErrors *grDataFitted = new TGraphErrors();

        TGraphErrors *grTheorPoint = new TGraphErrors();

        arma::vec yTheor = A*xVec;

        int iData = 0;
        for(int d = 0; d < fiel(0,3).n_rows; ++d) {
            double x  = fiel(0,3)(d,0);
            double Q2now = fiel(0,3)(d,1);
            double sig= fiel(0,3)(d,2);
            double err= fiel(0,3)(d,3);
            double theo= 0;

            for(int k = 0; k < A.n_cols; ++k)
                theo += fiel(k,3)(d,4) * xVec(k);

            //if(!Selector(x,Q2)) continue;
            bool isFitted = Selector(x,Q2);

            if(Q2now != Q2) continue;

             //sig;
            double errAbs = err*sig * 1e-2;
            if(isFitted) {
                int last = grDataFitted->GetN();
                grDataFitted->SetPoint(last, x, sig);
                grDataFitted->SetPointError(last, 0, errAbs);
            }
            else {
                int last = grData->GetN();
                grData->SetPoint(last, x, sig);
                grData->SetPointError(last, 0, errAbs);
            }


            //cout << "Data point added " << x <<" "<<sig << endl;
            grTheorPoint->SetPoint(iData, x, theo);

            ++iData;
        }

        TCanvas *can = new TCanvas("can", "canvas");
        gPad->SetLogx();

        double M = max({Max(grTheorPoint), Max(grData), Max(grDataFitted), Max(grTheor)});

        cout << "Max is " << M << endl;
        TH1F *fr = gPad->DrawFrame(1e-6, 0, 1, 1.6*M);
        fr->SetName(TString::Format("%g", M));
        fr->SetMaximum(1.6*M);


        grTheorPoint->SetLineColor(kRed);
        //grTheorPoint->Draw("l same");

        grData->SetMarkerStyle(24);
        grDataFitted->SetMarkerStyle(20);
        grData->Draw("pe same");
        grDataFitted->Draw("pe same");

        grTheor->Draw("l same");
        grF2->SetLineColor(kBlue);
        //grF2->Draw("l same");

        TLatex *lat = new TLatex();
        lat->DrawLatexNDC(0.5,0.5, TString::Format("%g", Q2));

        can->SaveAs(fname);

        delete can;
    }

};


int main()
{
    
    LinearFitter lfitter;
    for(int k = 3; k < 11; ++k) {
        lfitter.Load("../basis.dat", k);
        //return 0;
        lfitter.getMinimum();
        cout <<k<<" "<< lfitter.getChi2() << " / " << lfitter.getNdf()<< endl;
    }
    //return 0;

    //lfitter.PrintReduce();
    //lfitter.PrintPDF();

    lfitter.PrintDistrib("ahoj");
    return  0;

    int first = 0;
    int last = 46;

    for(int i = first; i < last; ++i) {
        if(i==first)  lfitter.PrintPDF(i, "haha.pdf(");
        if(i==last-1) lfitter.PrintPDF(i, "haha.pdf)");
        else          lfitter.PrintPDF(i, "haha.pdf");
    }

    return EXIT_SUCCESS;
}
