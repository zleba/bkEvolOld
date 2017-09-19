#include "TGraph.h"
#include "Solver.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "Fitter.h"

map<double, TGraph*> ReadFile(const char *fName);

int main(int argc, char **argv)
{

    //GetWeights(33);
    //return 0;

    /*
    Nodes nod(32+1, -log(1e6), +log(1e6));
    nod.CalcNodes();

    //cout <<setprecision(20)<< nod.Integrate([](double x) {return 1*1;}) << endl;
    cout <<"Hope "<< nod.Integrate(
    {
0.0316228 ,
0.031887 ,
0.0326923 ,
0.0340777 ,
0.0361114 ,
0.0388961 ,
0.0425774 ,
0.0473552 ,
0.0535009 ,
0.0613807 ,
 0.0714885 ,
 0.0844917 ,
 0.101296 ,
 0.123134 ,
 0.151694 ,
 0.189298 ,
 0.239155 ,
 0.305719 ,
 0.395203 ,
 0.516299 ,
 0.68122 ,
 0.907166 ,
 1.21842 ,
 1.64933 ,
 2.24845 ,
 3.08438 ,
 4.25362 ,
 5.89075 ,
 8.18218 ,
 11.3951 ,
 16.01 ,
 23.417 ,
 0,
    }



    ) << endl;
    return 0;
    */

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    assert(MPI_THREAD_FUNNELED == provided);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Print off a hello world message
    cout << "Processors " << processor_name << endl;
    



    auto pubSol32 = ReadFile("../BKsolver/bkresult32.dat");
    auto pubSol512 = ReadFile("../BKsolver/bkresult512.dat");
    //return 0;

    Fitter fitter;
    fitter.Init();

    return 0;

    double yNew = log(1e8/1e2);
    int Ny = 280;
    cout << "Starting " << endl;
    Solver sol512(512+1);
    sol512.InitF([](double x, double kT2) {
        //return pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        //return 1./pow(kT2,1);// pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        return kT2 * exp(-kT2);// * pow(max(0., 0.4-x), 2);
    });
    cout << "Weights calculated " << endl;
    //sol512.InitMat();

    //sol512.SaveEvolKernels("data/kernel");
    sol512.LoadEvolKernels("data/kernel");
    sol512.LoadConvKernels("data/kernel");
    //MPI_Finalize();
    //return 0;

    cout << "Matrix initialised" << endl;
    sol512.EvolveNew();
    //sol512.CalcF2L();

    cout << "Done " << endl;
    if(GetRankSize().first == 0) {
        //sol512.PrintBaseGrid();
        sol512.PrintReduce();
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    //sol512.PrintGrid();
    return 0;
    for(int i = 0; i < sol512.N; ++i) {
        const int Nnow = sol512.Nrap - 1;
        cout << "Haho " <<  exp(0.5*sol512.nod.xi[i]) <<" "<< sol512.PhiRapN[Nnow](i) << endl;
    }


    return 0;

    sol512.RunIterations(3, false);
    //for(int y = 0; y < Ny; ++y)
        //sol512.Step(yNew/Ny);


    cout << "Evolution done" << endl;
    //sol512.PrintGrid();


    /*
    return 0;


    Solver sol32(1024+1);
    sol32.InitF([](double kT2) {
        //return pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
        return 1./pow(kT2,0.25);// pow(1.0/sqrt(kT2) * exp(-pow(log(kT2/(1.*1.)),2)), 4);
    });
    sol32.InitMat();
    */

    TGraph *gr512 = new TGraph();
    //for(int i = 0; i < sol512.N; ++i) {
        //grOrg->SetPoint(i, exp(sol512.nod.xi[i]), sol512.Phi[i]);
    //}


    for(int i = 0; i < sol512.N; ++i) {
        const int Nnow = sol512.Nrap - 1;
        gr512->SetPoint(i, exp(sol512.nod.xi[i]), sol512.PhiRapN[Nnow](i));
        cout << "Haho " <<  exp(0.5*sol512.nod.xi[i]) <<" "<< sol512.PhiRapN[Nnow](i) << endl;
    }

    return 0;

    //TGraph *gr32 = new TGraph();
    //for(int i = 0; i < sol32.N; ++i) {
        //gr32->SetPoint(i, exp(sol32.nod.xi[i]), sol32.PhiRap[0][i]);
    //}


    TCanvas *can = new TCanvas("can","");
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogx();
    gPad->SetLogy();

    TGraph *grOrg = (TGraph*) pubSol512[0]->Clone("tabad");

    grOrg->Draw();
    grOrg->GetXaxis()->SetTitle("k^{2}_{T} [GeV^{2}]");
    grOrg->GetYaxis()->SetTitle("#phi");
    grOrg->SetMaximum(1e4);

    gr512->SetLineColor(kRed);
    gr512->Draw("same");



    pubSol32[0]->SetLineColor(kGreen);
    pubSol32[0]->SetLineStyle(kDashed);
    pubSol32[0]->Draw("same");

    pubSol32[5]->SetLineColor(kMagenta);
    pubSol32[5]->SetLineStyle(kDashed);
    pubSol32[5]->Draw("same");

    can->cd(2)->SetLogx();
    gPad->SetLeftMargin(0.23);
    TGraph *grRat32 = new TGraph();
    TGraph *grRat512 = new TGraph();
    //TGraph *grRatMy32 = new TGraph();
    TGraph *grRatMy512 = new TGraph();
    for(int i = 0; i < pubSol512[5]->GetN(); ++i) {
        double x, y;
        pubSol512[5]->GetPoint(i, x, y);
        double yPubl512 = pubSol512[5]->Eval(x);
        double yPubl32 = pubSol32[5]->Eval(x);
        //double yPubl = pubSol512[5]->Eval(x);

        double yMy512 = gr512->Eval(x);
        //double yMy32 = gr32->Eval(x);
        //double yPubl = 1./pow(x, 0.25);
        grRat512->SetPoint(i, x, yPubl512/yMy512);
        grRat32->SetPoint(i, x, yPubl32/yMy512);
        //grRatMy32->SetPoint(i, x, yMy32/yMy512);
        grRatMy512->SetPoint(i, x, 1.);
    }
    grRat32->GetYaxis()->SetTitle("#phi/#phi^{my}_{512}");
    grRat32->GetXaxis()->SetTitle("k^{2}_{T} [GeV^{2}]");
    grRat32->SetLineStyle(kDotted);
    grRat32->SetLineColor(kRed);
    grRat32->Draw("a c");
    grRat32->GetYaxis()->SetRangeUser(0.95, 1.05);

    grRat512->SetLineStyle(kDotted);
    grRat512->Draw("c same");
    //grRatMy32->SetLineColor(kRed);
    //grRatMy32->Draw("c same");
    grRatMy512->Draw("c same");

    TLegend *leg = new TLegend(0.3, 0.2, 0.7, 0.4);
    leg->SetHeader("#Delta y = 5");
    leg->SetBorderSize(0);
    leg->AddEntry(grRat32, "BK solver (32 points)", "l");
    leg->AddEntry(grRat512, "BK solver (512 points)", "l");
    //leg->AddEntry(grRatMy32, "My solver (Phi version) (512 points)", "l");
    leg->AddEntry(grRatMy512, "My solver (org version) (512 points)", "l");
    leg->Draw();
    

    can->SaveAs("ahoj.pdf");




    return 0;

#if 0

    gsl_cheb_series *cs = gsl_cheb_alloc(10);

    gsl_function F;

    F.function = f;
    F.params = 0;
    cout << "Here " << __LINE__ << endl;
    gsl_cheb_init(cs, &F, log(1e-5), log(1.0));
    cout << cs->c[0] << " "<< cs->c[1] << endl;

    double x = 0;//log(1e-5)/2.;
    double y = 2*(x - log(1e-5)) / ( log(1) - log(1e-5)) - 1;
    cout << " y " << y << endl;
    cout << "my gess " << 0.5*cs->c[0] +cs->c[1]*y   << endl;
    cout << "my gess " << 0.5*cs->c[0] +cs->c[1]*y + (2*y*y-1.0)*cs->c[2]
    + (4*y*y*y - 3*y)*cs->c[3]  
    
    << endl;
    double res = gsl_cheb_eval_n(cs, 3, x);
    cout << "res is : "<<res << endl;
    cout << "Here " << __LINE__ << endl;

    TGraph *gr = new TGraph();
    TGraph *grApp = new TGraph();


    int i = 0;
    cout << "Here " << __LINE__ << endl;
    for(double l = -5; l <=0; l += 0.01) {
        double x = pow(10,l);
        double res = gsl_cheb_eval(cs, log(x));

        grApp->SetPoint(i, x, res);
        gr->SetPoint(i, x, f(log(x)));
        //cout << "Here " << __LINE__ << endl;
        ++i;
    }

    /*
    TCanvas *can = new TCanvas("can","");
    gPad->SetLogx();
    gPad->SetLogy();

    //gr->Draw();

    grApp->SetLineColor(kRed);
    grApp->Draw();

    can->SaveAs("ahoj.pdf");
    */

    return 0;
    #endif
}
