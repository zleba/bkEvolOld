vector<dataPoint> Fitter::LoadData(string fname)
{
    vector<dataPoint> points;
    ifstream file(fname);
    char s[2000];
    file.getline(s,1990);
    file.getline(s,1990);

    int nPoints = 0, nSmallx=0;

    map<double> q2vals;
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

        double err2 = pow2(tot_noproc) + pow2(delta_rel) + pow2(delta_gp) + pow2(delta_had);
        for(int i = 0; i < 4; ++i)
            err2 += pow2(delta[i]);

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
        p.q2ID = q2vals[Q2];

    return points;

}
void Fitter::Init()
{
    //Load data points
    data = LoadData();

    //Load evoluton and convolution matrices



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

void AddTheory(vector<arma::vec> &F2, vector<arma::vec> &FL)
{
    for(auto &p : points) {
       double f2 = getValue(F2, p.q2ID, p.x);
       double fl = getValue(FL, p.q2ID, p.x);

       double yTerm = p.y*p.y/(1+(1-p.y)*(1-p.y));
       double sRed = f2 - yTerm * fl;

       p.theor = sRed;
    }
}

double getChi2()
{
    double chi2 = 0;
    for(auto &p : points) {
        if(p.x > 0.01 || p.Q2 < 4) continue;

        chi2 += pow((p.sigma - p.theor) / p.err, 2);
    }
    return chi2;
}
