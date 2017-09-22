#include "/home/zlebcr/libs/plottingHelper/plottingHelper.h"
//#include "/home/zlebcr/prog/bkEvol/inc/Fitter.h"
R__LOAD_LIBRARY(/home/zlebcr/libs/plottingHelper/libPlottingHelper.so)
//R__LOAD_LIBRARY(/home/zlebcr/prog/bkEvol/obj/Fitter.o)
using namespace PlottingHelper;

pair<map<double, TGraphErrors*>, map<double, TGraphErrors*>> ReadSigmaTheor(const char *fName)
{
    ifstream file(fName);

    map<double, TGraphErrors*> grMapX, grMapQ2;

    while(1) {
        string str;
        getline(file, str);
        if(!file.good()) break;
        if(str[0] == '#') continue;
        if(str.size() < 8) continue;
        stringstream  sStream(str);
        double Q2, sigma, x;
        sStream >> x >> Q2 >> sigma;

        if(grMapX.count(x) == 0) {
            grMapX[x] = new TGraphErrors();
        }
        if(grMapQ2.count(Q2) == 0) {
            grMapQ2[Q2] = new TGraphErrors();
        }

        grMapX[x]->SetPoint(grMapX[x]->GetN(), Q2, sigma);
        grMapX[x]->SetTitle(SF("x=%g", x));
        grMapX[x]->SetName(SF("x=%g", x));
        grMapQ2[Q2]->SetPoint(grMapQ2[Q2]->GetN(), x, sigma);
        grMapQ2[Q2]->SetTitle(SF("k_{T}^{2}=%g", Q2));

        //cout << kT2 <<" "<< b<<" "<<  Phi << endl;

    }
    return make_pair(grMapX, grMapQ2);

}

int cols[] = {1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9};

void PlotOverview(vector<map<double, TGraphErrors*>> grMap, TString mode, double Min, double Max)
{
    gPad->SetLogx();
    gPad->SetLogy();

    vector<double> vals;
    if(mode == "Q2") {
    for(auto &g : grMap[0])
        vals.push_back(g.first);
    }
    else {
    vals = {
    1.58e-05, 2e-05, 3.2e-05, 3.98e-05, 5e-05, 8e-05, 9.86e-05, 0.0001, 0.00013, 0.0002, 0.000251,
    0.00032, 0.0005, 0.0008, 0.001, 0.0011, 0.0013, 0.0015, 0.002, 0.0032, 0.005, 0.008, 0.013,
    0.02, 0.032, 0.05, 0.08, 0.13, 0.18, 0.25, 0.4, 0.65};
    }







    gPad->SetRightMargin(0.2);

    int start = vals.size()*0.00;
    int end = vals.size()*1.00;

    int step = 1;
    //if(vals.size() > 100) step = 50;
    //else                  step = 1;
    int iCol =0;
    cout <<"Vals size is " <<  vals.size() << endl;
    for(int i = start; i < end; i+=step) {
        cout <<"huhu " << vals[i] <<" "<< grMap[0].count(vals[i])<< endl;
        assert(grMap[0].count(vals[i]) != 0);

        if(mode == "Q2") {
            if(iCol == 0) {
                    grMap[0][vals[i]]->Draw("al");
            }
            else
                grMap[0][vals[i]]->Draw("l same");
        }
        else {
            if(iCol == 0)
                gPad->DrawFrame(0.1, Min,  3.3e4, Max);
            grMap[0][vals[i]]->Draw("l same");
        }

        grMap[0][vals[i]]->SetLineColor(cols[iCol%9]);

        if(grMap.size() == 2) {
            grMap[1][vals[i]]->SetLineColor(cols[iCol%9]);
            grMap[1][vals[i]]->SetLineStyle(2);
            grMap[1][vals[i]]->SetMarkerColor(cols[iCol%9]);
            grMap[1][vals[i]]->SetMarkerStyle(20);
            grMap[1][vals[i]]->Draw("ple same");

            for(int k = 0; k < grMap[1][vals[i]]->GetN(); ++k) {
                double x,y;
                grMap[1][vals[i]]->GetPoint(k, x, y);
                cout <<vals[i]<<" "<< k <<" "<< x<<" " << y << endl;
            }
        }
        ++iCol;
    }


    //grMap[vals[0]]->SetMinimum(Min);
    //grMap[vals[0]]->SetMaximum(Max);
    GetYaxis()->SetRangeUser(Min, Max);
    GetFrame()->SetTitle("");

    GetYaxis()->SetTitle("#sigma_{r}");
    if(mode == "Q2")
        GetXaxis()->SetTitle("x");
    else
        GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");

    //gPad->Update();



    TLegend *leg = new TLegend(0.8, 0.3, 1.0, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry((TObject*)0, "#sigma_r (HERA e^{+}p)", "h");

    leg->AddEntry((TObject*)0, "Solid - Data", "h");
    leg->AddEntry((TObject*)0, "Dashed - Fit", "h");

    for(int i = start; i < end; i+=step) {
        double val =  vals[i];
        if(mode == "Q2")
            leg->AddEntry(grMap[0][vals[i]], TString::Format("Q^{2} = %3.3g", val), "l");
        else
            leg->AddEntry(grMap[0][vals[i]], TString::Format("x = %3.3g", val), "l");
    }
    leg->Draw();


}

pair<map<double, TGraphErrors*>, map<double, TGraphErrors*>> ReadSigmaData(string fname);

void plotRed()
{

    map<double, TGraphErrors*> theorX, theorQ2;
    map<double, TGraphErrors*> dataX, dataQ2;

    //tie(theorX, theorQ2) = ReadSigmaTheor("../ sigmaRed");
    tie(theorX, theorQ2) = ReadSigmaTheor("fitOut");
    tie(dataX, dataQ2) = ReadSigmaData("../test/heraTables/HERA1+2_NCep_920.dat");

    TCanvas *can = new TCanvas("can", "canvas");
    PlotOverview({dataX}, "x", 0.09, 2);

    TCanvas *dan = new TCanvas("dan", "canvas");
    PlotOverview({theorQ2, dataQ2}, "Q2", 1e-1, 2);

}


struct dataPoint {
    double x, Q2, y;
    double sigma, err;

    int q2ID;
    double theor;
};


vector<dataPoint> LoadData(string fname)
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

        p.err = sqrt(err2);
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




pair<map<double, TGraphErrors*>, map<double, TGraphErrors*>> ReadSigmaData(string fname)
{
    auto data = LoadData(fname);

    map<double, TGraphErrors*> grMapX, grMapQ2;

    for(auto &p : data) {
        double x = p.x;
        double Q2 = p.Q2;
        if(grMapX.count(x) == 0) {
            grMapX[x] = new TGraphErrors();
            cout << "New x val " << x << endl;
        }
        if(grMapQ2.count(Q2) == 0) {
            grMapQ2[Q2] = new TGraphErrors();
        }

        double err = p.sigma * p.err * 1e-2;
        int lastX = grMapX[x]->GetN();
        //cout << " p.err " << p.err << endl;
        grMapX[x]->SetPoint(lastX, Q2, p.sigma);
        grMapX[x]->SetPointError(lastX, 0., err);
        grMapX[x]->SetTitle(SF("x=%g", x));
        grMapX[x]->SetName(SF("x=%g", x));

        int lastQ2 = grMapQ2[Q2]->GetN();
        grMapQ2[Q2]->SetPoint(lastQ2, x, p.sigma);
        grMapQ2[Q2]->SetPointError(lastQ2, 0., err);
        grMapQ2[Q2]->SetTitle(SF("k_{T}^{2}=%g", Q2));

    }
    //cout << kT2 <<" "<< b<<" "<<  Phi << endl;
    cout << "Holka test " << grMapX.size() << " "<< grMapQ2.size() << endl;

    for(auto &gr : grMapX ) { 
        cout <<"Helenka " <<setprecision(10)<<  gr.first <<" "<< gr.second->GetN() << endl;
    }


    return make_pair(grMapX, grMapQ2);

}


