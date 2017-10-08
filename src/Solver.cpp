#include "Solver.h"


//--------------------------------------------------------------------
void Solver::InitMat()
{

    //matN.resize(Nrap, arma::mat(N,N,arma::fill::zeros));
    //matNDiag.resize(Nrap, arma::mat(N,N,arma::fill::zeros));



    matN.zeros( N, N, Nrap);
    matNDiag.zeros( N, N, Nrap);
    matNInv.zeros( N, N, Nrap);

    int fac = (Nint-1)/(N-1);

    int start, end;
    tie(start,end) = GetStartEnd(0, Nrap-1);




    cout << "Start+end|nrap " << start <<" "<< end <<"|" << Nrap<< endl;

    double stepY = (rapMax - rapMin) / (Nrap-1);
    //for(int y = 0; y < Nrap; ++y) { 
    for(int y = start; y <= end; ++y) { 
        double rap = rapMin + stepY * y;
        double z = exp(-rap);

        arma::mat mTemp(Nint,Nint,arma::fill::zeros);
        arma::mat mDiagTemp(Nint,Nint,arma::fill::zeros);
        cout << "Matrix init y " << y << endl;

#pragma omp parallel for
        for(int i = 0; i < Nint; ++i)   //loop over L
            for(int j = 0; j < Nint; ++j) { //loop over L' (integral)
                double L  = nod.xi[i];
                double Lp = nod.xi[j];
                double w = nod.wi[j];

                //double l  = exp(0.5*L);
                //double lp = exp(0.5*Lp);
                double l  = nod.SqrtExpXi[i];
                double lp = nod.SqrtExpXi[j];

                //mat[y][i][j] += Kernel85(l, lp, z) * w;
                //mat[y][i][i] += Kernel85Diag(l, lp, z) * w;

                //mat[y][i][j] += KernelSub81(l, lp, z) * w;
                //mat[y][i][i] += KernelSub81Diag(l, lp, z) * w;

                //mat[y][i][j] += KernelSub88(l, lp, z) * w;

                //matDiag[y][i][j] = Kernel85zDiag(l, lp, z) * w; //Diag-z DGLAP part

                if(!toTrivial || i % fac == 0) {

                    mTemp(i,j) += KernelBFKL(l, lp, z) * w;
                    mTemp(i,i) += KernelBFKLDiag(l, lp, z) * w;
                    //mTemp(i,j) += KernelSub79(l, lp, z) * w;
                    //mTemp(i,i) += KernelSub79Diag(l, lp, z) * w;
                }


                //mTemp(i,j) += Kernel85(l, lp, z) * w;
                //mTemp(i,i) += Kernel85Diag(l, lp, z) * w;
                //mDiagTemp(i,j) = Kernel85zDiag(l, lp, z) * w;

                //Clasicall BFKL
                /*
                   if(i != j)
                   mat[y][i][j] += Kernel(L, Lp) * w;

                   double kerDs, kerDr;
                   tie(kerDs, kerDr) = KernelDiag(L, Lp);
                   if(i != j)
                   mat[y][i][i] += (kerDs+kerDr) * w;
                   else
                   mat[y][i][i] += (kerDr) * w;
                   */

            }
        //cout << redMat << endl;
        //cout << extMat << endl;


        //cout << extMat << endl;
        //cout << redMat << endl;
        //exit(0);

        matN.slice(y)     = stepY * redMat * mTemp * extMat;
        matNDiag.slice(y) = redMat * mDiagTemp * extMat;




        /*
           for(int i = 0; i < Nint; ++i)  //loop over L
           for(int j = 0; j < Nint; ++j) { //loop over L' (integral)
           double m  = matDiag[y][i][j];
           double mt = mDiagTemp(i,j);
           assert(isfinite(mat[y][i][j]));
           assert(isfinite(matDiag[y][i][j]));
        //cout << "Ihned "<<y<<" : " <<i<<" "<<j<<" "<< m <<" "<< mt << " "<< 2*(m-mt)/(m+mt) <<endl;
        //cout << "Mdiag "<<y<<" : " <<i<<" "<<j<<" "<< matDiag[y][i][j] <<endl;
        }
        */

        //if(y == 1) exit(0);

    }

    cout << "Reduce start" << endl;
    //Merge things together
    MPI_Allreduce(MPI_IN_PLACE, matN.memptr(), matN.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, matNDiag.memptr(), matNDiag.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for(int y = start; y <= end; ++y) { 
        matNInv.slice(y) = inv(arma::mat(N,N,arma::fill::eye) - 0.5*matN.slice(0) -matNDiag.slice(y));
    }
    MPI_Allreduce(MPI_IN_PLACE, matNInv.memptr(), matNInv.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    cout << "Reduce done" << endl;




    //matN.save("ahoj.hdf5", arma::hdf5_binary);
    //MPI_Finalize();
    //exit(0);


    //for(int y = 0; y < Nrap; ++y) { 
    //MPI_Allreduce(MPI_IN_PLACE, matN.slice(y).memptr(),N*N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE, matNDiag.slice(y).memptr(),N*N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //}

    /*
       for(int y = start; y <=end; ++y) {
       int rank = GetRankSize().first;
       MPI_Bcast(matN[y].memptr(), N*N, MPI_DOUBLE, rank, MPI_COMM_WORLD);
       MPI_Bcast(matNDiag[y].memptr(), N*N, MPI_DOUBLE, rank, MPI_COMM_WORLD);
       }
       */

    //exit(0);

    /*
       if(start == 0) {
       for(int y = 0; y < Nrap; ++y) {
       const double stepY = (rapMax - rapMin) / (Nrap-1);
       for(int i = 0; i < N; ++i)
       for(int j = 0; j < N; ++j) {
       cout <<"Kernel "<<y<<" : "<< i <<" "<< j <<" "<<setprecision(10)<< matN[y](i,j) <<" "<< matNDiag[y](i,j)<<" "<< stepY<< endl;
    //cout <<"Kernel "<<y<<" : "<< i <<" "<< j <<" "<<setprecision(10)<< matMPI[y](i,j) <<" "<< matMPIDiag[y](i,j)<<" "<< stepY<< endl;
    }
    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
    */




    /*
       exit(0);
       */
}






//--------------------------------------------------------------------
void Solver::EvolveNew()
{
    const double stepY = (rapMax - rapMin) / (Nrap-1);

    F2rap.resize(Nrap);
    FLrap.resize(Nrap);

    //In case of convF2, convFL not filled
    if(convF2.n_rows == 0) {
        cout << "Dummy init of the convolution cube" << endl;
        convF2.zeros(46, N,  Nrap);
        convFL.zeros(46, N,  Nrap);
    }


    bool doGPU = true;
    if(doGPU & !gpu.isInited) gpu.InitAll(matN, convF2, convFL);
    if(doGPU) gpu.ResetVector();

    int start = 0;

    //Classical approach
    if(start == 0) {
        arma::mat MatEq = arma::mat(N,N,arma::fill::eye) -  matNDiag.slice(0);
        if(putZero) {
            PhiRapN[0] = arma::vec(N, arma::fill::zeros);
        }
        else
            PhiRapN[0] = GetLinSolution(MatEq, Phi0N[0]);

        F2rap[0].zeros(convF2.n_rows);
        FLrap[0].zeros(convFL.n_rows);
    }
    else { //Dummy start for DGLAP
        for(int y = 0; y <= start; ++y)
            PhiRapN[y] = Phi0N[y];
    }


    for(int y = start+1; y < Nrap; ++y) {
        //Starting point of evol with 0.5 (Trapezius)

        arma::vec yTemp(N, arma::fill::zeros);
        arma::vec myVec(N, arma::fill::zeros);

        F2rap[y] = 0.5*convF2.slice(y) * PhiRapN[0];// + convF2.slice(0) * PhiRapN[y]);
        FLrap[y] = 0.5*convFL.slice(y) * PhiRapN[0];// + convFL.slice(0) * PhiRapN[y]);

        //openMPI treatment
        int start=1, end;

        if(!doGPU) {
            tie(start,end) = GetStartEnd(1, y-1); //from 1 to y-1
            //cout << "Start+end|nrap " << start <<" "<< end <<"|" << y-1<< endl;
            //Remaining without mult
            for(int d = start; d <= end; ++d)
                yTemp += matN.slice(d) * PhiRapN[y-d];

            MPI_Allreduce(MPI_IN_PLACE, yTemp.memptr(), yTemp.n_elem,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }
        else {
            if(y > 1) {
                gpu.ConvoluteAll(y);
                //gpu.GetResult(y, yTemp);
                gpu.GetResults(y, yTemp, F2rap[y], FLrap[y]);

            }

            /*
               cout << "Radek start "<<y << endl;
            //cout << myVec << yTemp << endl;
            for(int i = 0; i < myVec.n_elem; ++i)
            cout << i <<" " << myVec(i) <<" "<< yTemp(i) << endl;

            cout << "Radek end " <<y<< endl;
            if(y > 2) exit(0);
            */
        }






        yTemp += 0.5 * matN.slice(y) * PhiRapN[0];


        //Whole right hand side
        //yTemp = stepY * yTemp + Phi0N[y];
        yTemp += Phi0N[y];

        if(y % 100 == 0)
            cout <<"Rap point " << y << endl;
        /*
           auto matNow = mat[0]; //Adding diagonal DGLAP term
           for(int k = 0; k < N; ++k)
           for(int l = 0; l < N; ++l)
           matNow[k][l] += matDiag[y][k][l];
           */

        //arma::mat matEq = arma::mat(N,N,arma::fill::eye) - 0.5*stepY*matN.slice(0) -matNDiag.slice(y);

        //vector<double> Pred = PhiRap[y-1];
        //if(y >= 4)
        //for(int k = 0; k < N; ++k)
        //Pred[k] += PhiRap[y-1][k] - PhiRap[y-2][k];


        //PhiRap[y] = GetRegSolution(matEq, yTemp, Pred);
        //PhiRapN[y] = GetLinSolution(matEq, yTemp);
        PhiRapN[y] = matNInv.slice(y) * yTemp;

        F2rap[y] += 0.5*convF2.slice(0) * PhiRapN[y];
        FLrap[y] += 0.5*convFL.slice(0) * PhiRapN[y];

        if(doGPU) gpu.SetPhi(y, PhiRapN[y]);

        //PhiRapN[y] =  IterSolution(matEq, double factor, const arma::vec &y);


        bool isGood = true;

        //cout << "Id of y = " << y << endl;
        if(!isGood && 0) {

            for(int j = 0; j < 10; ++j)
                cout << y <<" " << j <<" "<< PhiRapN[y](j) << endl;

        }
    }

}

//--------------------------------------------------------------------
void Solver::CalcF2L()
{
    F2rap.resize(Nrap);
    FLrap.resize(Nrap);
    F2rap[0] = arma::vec(convF2.n_rows, arma::fill::zeros);
    FLrap[0] = arma::vec(convFL.n_rows, arma::fill::zeros);
    for(int y = 1; y < Nrap; ++y) {
        //kT spectrum for particular bin y
        //Starting point of evol with 0.5 (Trapezius)
        F2rap[y] = 0.5*(convF2.slice(y) * PhiRapN[0] + convF2.slice(0) * PhiRapN[y]);
        FLrap[y] = 0.5*(convFL.slice(y) * PhiRapN[0] + convFL.slice(0) * PhiRapN[y]);
        //Remaining without mult by 0.5, convF2 should contain deltaRap factor
        for(int d = 1; d < y; ++d) {
            F2rap[y] += convF2.slice(d) * PhiRapN[y-d];
            FLrap[y] += convFL.slice(d) * PhiRapN[y-d];
        }
    }
}







//--------------------------------------------------------------------
void Solver::DoIteration()
{
    const double stepY = (rapMax - rapMin) / (Nrap-1);

    //vector<vector<double>> PhiRapNew(Nrap);
    vector<arma::vec> PhiRapNew(Nrap);


    PhiRapNew[0] = Phi0N[0] + matNDiag[0] * PhiRapN[0];

    for(int y = 1; y < Nrap; ++y) {

        //kT spectrum for particular bin y
        //Starting point of evol with 0.5 (Trapezius)
        arma::vec yTemp = 0.5*(matN.slice(y) * PhiRapN[0] + matN.slice(0) * PhiRapN[y]);

        //Remaining without mult
        for(int d = 1; d < y; ++d) {
            yTemp += matN.slice(d) * PhiRapN[y-d];
        }

        //Whole right hand side
        yTemp = yTemp + Phi0N[y];

        //Diag part (=virtual DGLAP term)
        yTemp += matNDiag.slice(y) * PhiRapN[y];

        PhiRapNew[y] = yTemp;
    }
    PhiRapN = PhiRapNew;
}


//--------------------------------------------------------------------
void Solver::RunIterations(int Niter, bool init)
{
    //Init PhiRap 
    if(init) {
        assert(Phi0N.size() == PhiRapN.size());
        for(int y = 0; y < PhiRapN.size(); ++y)
            PhiRapN[y] = Phi0N[y];
    }

    //Do itrerations
    for(int i = 0; i < Niter; ++i) {
        DoIteration();
        cout << "Iteration " << i <<" done." << endl;
        cout << "Phi[kT0] = " << PhiRapN[Nrap-1](0) << endl;
    }
}


//--------------------------------------------------------------------
arma::vec Solver::GetRHS(const arma::vec &PHI)
{
    arma::vec dPhi = matN.slice(0) * PHI;
    return dPhi;
}


//--------------------------------------------------------------------
void Solver::Step(double delta) {
    static int y = 0;
    if(y == 0) PhiRapN[0] = Phi0N[0];

    arma::vec dPhi = GetRHS(PhiRapN[0]);
    ++y;
    PhiRapN[0] += delta * dPhi;

}


//Function of x and kT2
//--------------------------------------------------------------------
void Solver::InitF(function<double(double, double)> fun) {
    const double stepY = (rapMax-rapMin) / (Nrap-1);

    //New version
    Phi0N.resize(Nrap);
    for(int y = 0; y < Nrap; ++y) {
        arma::vec temp(Nint);
        double x = exp(-y*stepY);
        for(int i = 0; i < Nint; ++i) {
            double kT2 = exp(nod.xi[i]);
            temp(i) = fun(x, kT2);
        }
        Phi0N[y] = redMat * temp;
    }

    PhiRapN.resize(Nrap, arma::vec(N, arma::fill::zeros));
}

//SetSolution of x and kT2
//--------------------------------------------------------------------
void Solver::SetSolution(function<double(double, double)> fun)
{
    const double stepY = (rapMax-rapMin) / (Nrap-1);

    //New version
    PhiRapN.resize(Nrap);
    for(int y = 0; y < Nrap; ++y) {
        arma::vec temp(Nint);
        double x = exp(-y*stepY);
        for(int i = 0; i < Nint; ++i) {
            double kT2 = exp(nod.xi[i]);
            temp(i) = fun(x, kT2);
        }
        PhiRapN[y] = redMat * temp;
        //cout << redMat << endl;
        //exit(0);
    }
}



//--------------------------------------------------------------------
double Solver::Interpolate(double y, double L)
{

    y = max(y, rapMin);
    y = min(y, rapMax);

    double stepY = (rapMax-rapMin) / (Nrap-1);
    int yId = (y - rapMin) / stepY;
    int LId = 0;
    for(LId = N-1; LId >= 0; --LId)
        if(nodBase.xi[LId] <= L) break;
    //--LId;
    LId = max(0, LId);
    yId = min(Nrap-2, yId);

    double yLeft = rapMin + yId*stepY;
    double LLeft = nodBase.xi[LId];
    double LRight = nodBase.xi[LId+1];


    //if (LRight - L < 0 && L - LRight < 1e-7)
    //L = LRight;

    assert(y - yLeft >= 0);
    assert(yLeft+stepY - yLeft >= 0);



    if (LRight - L < 0 || L - LLeft < 0) {
        cout <<"Error " <<setprecision(24)<< L <<" : "<< LLeft <<" "<< LRight << endl;
        //L = LRight;
        //cout << "Ahojky " << nod.xi[0] << " "<< nod.xi[Nrap-1] << endl;
    }
    assert(L - LLeft >= 0);
    assert(LRight - L >= 0);

    if (yId > Nrap - 2) {
        cout << "Ahoj " << y << " :  " << yLeft << " "<< yLeft + stepY << endl;
    }

    assert(yId <= Nrap - 2);


    assert(LId <= N - 2);
    assert(yId >= 0);
    assert(LId >= 0);

    double fLL, fLR, fRL, fRR;

    fLL = PhiRapN[yId](LId);
    fLR = PhiRapN[yId](LId+1);
    fRL = PhiRapN[yId+1](LId);
    fRR = PhiRapN[yId+1](LId+1);

    bool canBeLog = fLL > 0 && fLR > 0 && fRL > 0 && fRR > 0;

    if(canBeLog) {
        fLL = log(fLL);
        fLR = log(fLR);
        fRL = log(fRL);
        fRR = log(fRR);
    }

    //cout << exp(fLL) << " "<< exp(fLR) <<" : "<< exp(fRL) << " "<< exp(fRR)<<  endl;

    double sum = 
        (y - yLeft)*(L - LLeft)*fRR  + (yLeft + stepY - y)*(L - LLeft)*fLR +
        (y - yLeft)*(LRight - L)*fRL + (yLeft + stepY - y)*(LRight - L)*fLL;
    sum /= stepY * (LRight - LLeft);
    assert(stepY != 0);
    assert((LRight-LLeft) != 0);

    if(canBeLog)
        return exp(sum);
    else 
        return sum;

}



//--------------------------------------------------------------------
void Solver::PrintBaseGrid()
{
    double stepY = (rapMax - rapMin) / (Nrap-1);
    double stepL = (Lmax - Lmin) / (N-1);


    for(int y = 0; y < Nrap; ++y)
        for(int l = 0; l < N; ++l) {
            double yNow = rapMin + y*stepY;
            double L  = nodBase.xi[l];

            double kt2 = exp(L);
            double x = exp(-yNow);

            cout << x << " "<< sqrt(kt2) <<"  " <<  PhiRapN[y](l) << endl;
        }
}

//--------------------------------------------------------------------
void Solver::PrintGrid()
{
    double stepY = (rapMax - rapMin) / (Nrap-1);

    //cout <<"Interpolation "<< Interpolate(0, log(1) ) << endl;
    //return;

    double yStepMich = log(1e-2/1e-8) / 99;//(Nrap -1);
    double LStepMich = log(1e6/1e-2)  / 99;// (N -1);

    cout << "Radek " << endl;
    //for(int y = 0; y < 100; ++y) 
    for(int y = 99; y >= 0; --y) {
        double x = 1e-2 * pow(1e-8/1e-2, y/99.);
        double yNow = y*yStepMich;
        for(int i = 0; i < 100; ++i) {
            double kT2 = 0.01 * pow(1e6/0.01, i/99.);
            double LNow = i*LStepMich + log(0.01);
            if(i == 99) LNow = log(1e6)-1e-9;
            double res = Interpolate(yNow, LNow);
            /*if(yNow == 0)*/ cout << x<<" "<< kT2<<" "<< res <<  endl;
        }
    }
}


//--------------------------------------------------------------------
void Solver::PrintReduce()
{
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
            double sRed = F2rap[rapID](i) - y*y/yPlus * FLrap[rapID](i);

            //cout << x << " "<< Q2 <<" "<< sRed << endl;
            cout << x << " "<< Q2 <<" "<< F2rap[rapID](i) <<" "<< FLrap[rapID](i) << endl;

        }
    }
}


vector<double> Solver::GetRegSolution(const vector<vector<double>> &MatEq, const vector<double> &y, const vector<double> &yReg)
{
    //cout 
    arma::mat M(N,N);
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; ++j){
            M(i,j) = MatEq[i][j];
        }

    arma::vec v(N), vReg(N);
    for(int i = 0; i < N; ++i) {
        v(i) = y[i]; 
        vReg(i) = 0;//yReg[i]; 
    }
    arma::mat L(N-1,N);

    for(int i = 0; i < N-1; ++i) 
        for(int j = 0; j < N; ++j) {
            if(i == j)
                L(i,j) = 1.;///(pow(yReg[i],2) + 1e-10);
            if(i == j - 1)
                L(i,j) =-1.;///(pow(yReg[i],2) + 1e-10);
            else
                L(i,j) = 0.;///(pow(yReg[i],2) + 1e-10);
        }
    arma::mat L2 = trans(L)*L;

    //double tau = 3e-1;
    double tau = 5e-1;
    arma::mat E = arma::trans(M)*M + tau*L2;
    //for(int i = 0; i < N; ++i)
    //E(i,i) += tau * L(i,i);

    //cout << " " << Ep(3,3) << " "<< Ep(3,4) << endl;
    //arma::mat E = arma::trans(M)*M + tau;
    //cout << " " << E(3,3) << " "<< E(3,4) << endl;
    //assert(0);
    arma::vec r = arma::trans(M)*v + tau*L2*vReg;

    arma::vec res = arma::solve(E,r);

    vector<double> Res(N);
    for(int i = 0; i < N; ++i)
        Res[i] = res[i];
    return Res;

}



arma::vec Solver::IterSolution(const arma::mat &Mat, double factor, const arma::vec &y)
{
    arma::vec yNow = y;
    arma::vec ySum = y;

    double diff = 0;
    for(int i = 0; i < 20; ++i) {
        yNow = factor * Mat * yNow;

        diff = 0;
        for(int k = 0; k < N; ++k)
            diff = max(diff, abs(yNow[k]) / (1e-13+abs(yNow(k)) + abs(ySum(k))));
        //cout <<"Diff " <<  i << " "<< diff << endl;

        ySum += yNow;

        //cout << "Diff is " << diff << endl;
        if(diff < 1e-9) break;
    }
    assert(diff < 1e-5);


    return ySum;
}

