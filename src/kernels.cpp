#include "Solver.h"
#include "integration.h"
#include <iomanip>

double IntegralCos(double a, double x)
{
    double p = sqrt(1 - a*a); 

    if(x == M_PI)
        return 2.*M_PI/2. / p; //In case of whole half cyrcle

    double res = 2*atan( (a+1)*tan(x/2) / p) / p;
    return res;
}

pair<double,double> Solver::GetKerPar(double l, double lp)
{
    double r = l / lp;
    double ker = 1./(2*M_PI) * 1./(1 + r*r);
    double par =  2*r / (1 + eps + r*r);
    return {ker, par};
}

//Only for l != lp
//Integral of  lp^2/2PI * 1/(l^2+lp^2 - 2*l*lp*cos(phi)) dphi
double IntegralPhi(double l, double lp, double angle)
{
    assert(l != lp);

    double Int;
    if (angle < M_PI)
        Int = 2 * atan((l+lp)/(l-lp) * tan(angle/2)) / (l*l - lp*lp);
    else
        Int = M_PI / abs(l*l - lp*lp);

    return Int * lp*lp/ (2*M_PI);
}

//Only for l == lp
//Integral of  l^2/2PI * 1/(l^2+l^2 - 2*l*l*cos(phi)) dphi
//Measured from angle to M_PI
double IntegralPhiDiag(double l, double lp, double angle)
{
    assert(l == lp);
    assert(angle <= M_PI);
    if(angle == M_PI) return 0;

    double Int = 1./tan(angle/2.);
    return Int / (4*M_PI);
}


double GetAngle(double l, double lp, double cond)
{
    double Cos = (l*l + lp*lp - cond) / (2*l*lp);
    if(Cos <= -1) return M_PI;
    if(Cos >=  1) return 0;
    return acos(Cos);
}

//Arguments are always momenta l and lp in GeV, and z
//
//

////////////////////////////////////////////////
//Equation 79 (with eps reg)
////////////////////////////////////////////////

//Form of equation with phi
double Solver::Kernel79(double l, double lp, double z)
{
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double Int = 2 * IntegralCos(par, M_PI);
    return as * ker *  Int; 
}

//Form of bfkl with phi
double Solver::Kernel79Diag(double l, double lp, double z)
{
    if(lp > 2*l) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = acos(lp/(2.*l));
    
    double Int = 2 * IntegralCos(par, angleMax);
    return -as * ker * Int;

}

//Form of equation with phi
double Solver::KernelSub79(double l, double lp, double z)
{
    if (l == lp) return 0;
    double Int = 2 * IntegralPhi(l, lp, M_PI); //for uper and lower "hemisphere"
    return as * Int; 
}

//Form of bfkl with phi
double Solver::KernelSub79Diag(double l, double lp, double z)
{
    double angleMax = GetAngle(l,lp, l*l);
    if(angleMax == 0) return 0;
    
    double Int;
    if (l != lp)
        Int =  2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l,lp, angleMax);

    return -as * Int;
}





////////////////////////////////////////////////
//Equation 80 (with eps reg)
//Take care of mu2 dependence
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Solver::Kernel80(double l, double lp, double z)
{
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double Int = 2 * IntegralCos(par, M_PI);

    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    return as * ker *  Int; 
}

//Diagonal kernel with F(k)
double Solver::Kernel80Diag(double l, double lp, double z)
{
    if(lp > 2*l) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = acos(lp/(2.*l));

    double Int = 2 * IntegralCos(par, angleMax);

    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    return -as * ker * Int;

}


//Off-diagonal kernel with F(k+q)
double Solver::KernelSub80(double l, double lp, double z)
{
    if (l == lp) {
        double angleMin = GetAngle(l,lp, mu2);
        return as * 2 * IntegralPhiDiag(l, lp, angleMin);
    }

    double Int = 2 * IntegralPhi(l, lp, M_PI);

    // theta(q2 - mu2) term
    if (pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    return as * Int; 
}

//Diagonal kernel with F(k)
double Solver::KernelSub80Diag(double l, double lp, double z)
{
    double angleMax = GetAngle(l,lp, l*l);
    if(angleMax == 0) return 0;

    if(l == lp) {
        double angleMin = GetAngle(l,lp, mu2);
        if(angleMax <= angleMin) return 0;
        return -2*as * (IntegralPhiDiag(l, lp, angleMin) - IntegralPhiDiag(l, lp, angleMax));
    }

    double Int = 2 * IntegralPhi(l, lp, angleMax);
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    return -as * Int;
}




////////////////////////////////////////////////
//Equation 81 (with eps reg)
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Solver::Kernel81(double l, double lp, double z)
{
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax = M_PI;
    // theta(k2/q2 - z)
    double Cos = (lp*lp + l*l -l*l/z) / (2*l*lp);
    if(Cos >= 1) return 0;
    if(Cos > -1) angleMax = acos(Cos);

    double Int = 2 * IntegralCos(par, angleMax);

    //cout << ker << " "<< Int << " "<< angleMax <<" "<< Cos << " "<< l<<" "<<lp<<" "<< l*l/pow(l-lp,2) -z<<  endl;
    //assert(isfinite(as*ker*Int));

    return as * ker *  Int; 
}

//Diagonal kernel with F(k)
double Solver::Kernel81Diag(double l, double lp, double z) {
    return Kernel79Diag(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double Solver::KernelSub81(double l, double lp, double z)
{
    double angleMax = GetAngle(l,lp, l*l/z);
    if(angleMax == 0) return 0;

    double Int;
    if(l != lp)
        Int = 2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l, lp, angleMax);

    return as * Int; 
}

//Diagonal kernel with F(k)
double Solver::KernelSub81Diag(double l, double lp, double z) {
    return KernelSub79Diag(l, lp, z);
}


////////////////////////////////////////////////
//Equation 82 (with eps reg) -- Equivalent of 81
////////////////////////////////////////////////


////////////////////////////////////////////////
//Equation 83 (with eps reg) 
//Maybe not properly defined
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Solver::Kernel83(double l, double lp, double z)
{
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    if(z == 1) return 0;
    double angleMax = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax = acos(Cos);

    double Int = 2 * IntegralCos(par, angleMax);

    return as * ker *  Int; 
}


//Diagonal kernel with F(k)
double Solver::Kernel83Diag(double l, double lp, double z) {
    if(z == 1) return 0;
    return Kernel79Diag(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double Solver::KernelSub83(double l, double lp, double z)
{
    if(z == 1) return 0;

    // theta(k2/q2 - z/(1-z))
    double angleMax = GetAngle(l,lp, l*l*(1-z)/z);

    double Int;
    if(l != lp)
        Int =  2 * IntegralPhi(l, lp, angleMax);
    else
        Int = -2 * IntegralPhiDiag(l, lp, angleMax);

    return as *  Int; 
}


//Diagonal kernel with F(k)
double Solver::KernelSub83Diag(double l, double lp, double z) {
    if(z == 1) return 0;
    return KernelSub79Diag(l, lp, z);
}




////////////////////////////////////////////////
//Equation 84 (with eps reg) 
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Solver::Kernel84(double l, double lp, double z)
{
    Kernel83(l, lp, z);
}

//Diagonal kernel with F(k)
double Solver::Kernel84Diag(double l, double lp, double z)
{
    if(lp > 2*l) return 0;
    if(z == 1) return 0;

    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    double angleMax1 = acos(lp/(2.*l));
    
    double angleMax2 = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax2 = acos(Cos);

    double Int = 2 * IntegralCos(par, min(angleMax1, angleMax2));
    return -as * ker * Int;
}

//Off-diagonal kernel with F(k+q)
double Solver::KernelSub84(double l, double lp, double z)
{
    KernelSub83(l, lp, z);
}

//Diagonal kernel with F(k)
double Solver::KernelSub84Diag(double l, double lp, double z)
{
    if(z == 1) return 0;

    double angleMax1 = GetAngle(l,lp, l*l);
    if(angleMax1 == 0) return 0;

    double angleMax2 = GetAngle(l,lp, l*l*(1-z)/z);
    if(angleMax2 == 0) return 0;

    double Int;
    if (l != lp)
        Int = 2 * IntegralPhi(l, lp, min(angleMax1, angleMax2));
    else
        Int = -2 * IntegralPhiDiag(l, lp, min(angleMax1, angleMax2));

    return -as * Int;
}





//Return z/6*Pgg - 1
double PggMod(double z)
{
    double reg = z*z*(1-z) - (2*z+1);
    //return reg;
    //reg = 0;
    if(z == 1) return reg;
    else       return reg + 1./(1-z);
}

double Solver::DGLAPterm(double l, double lp, double z)
{
    double ker = lp*lp/(2*M_PI) *  1/(l*l); //constant before

    //q < l
    if(lp > 2*l) return 0.;
    double angleMax = acos(lp/(2.*l));
    double Int = 2 * angleMax;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * angleMin;
    }
    double res = PggMod(z) * as * ker * Int;

    return res;
}

////////////////////////////////////////////////
//Equation 85 (with eps reg) 
//BFKL with DGLAP
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Solver::Kernel85(double l, double lp, double z)
{
    //Adding normal BFKL
    double res =  Kernel79(l, lp, z);

    /*
    double ker = lp*lp/(2*M_PI) *  1/(l*l); //constant before

    //q < l
    if(lp > 2*l) return res;
    double angleMax = acos(lp/(2.*l));
    double Int = 2 * angleMax;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * angleMin;
    }
    res += PggMod(z) * as * ker * Int;
    */
    res += DGLAPterm(l, lp, z);

    return res;
}

double Solver::Kernel85Diag(double l, double lp, double z)
{
    //return 0;
    return Kernel79Diag(l, lp, z);
}

//Off-diagonal z-diaginal kernel with F(k+q)
double Solver::Kernel85zDiag(double l, double lp, double x)
{
    putZero = true;
    double ker = lp*lp/(2*M_PI) *  1/(l*l); //constant before

    //q < l
    if(lp > 2*l) return 0;
    double angleMax = acos(lp/(2.*l));
    double Int = 2 * angleMax;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * angleMin;
    }

    const int nf = 4;

    double rap = log(1/x);
    double nStepReal = (rap-rapMin)/(rapMax-rapMin) * (Nrap - 1) + 1;
    int nStep = round(nStepReal);
    assert(abs(nStepReal-nStep) < 1e-2);
    double stepSize = (rapMax - rapMin) / (Nrap - 1);
    double sum = 0;
    double rapNow = rap;
    for(int i = 0; i < nStep; ++i, rapNow -= stepSize) {
        if(i == 0) continue;//skip z == 1
        double step = (i == nStep - 1) ? stepSize/2 : stepSize;
        double z = exp(rapNow - rap);
        sum += -z/(1-z) * step;
        //sum += -1/(1-z) * step;
        //cout << "RADEK " << x <<" "<< z << endl;
    }

    //assert(x != 1);
    //cout << "x val = " << x << endl;
    sum += (33 - 2*nf) / 36.0;

    if(nStep != 1)
        sum += log(1 - x/exp(-rapMin));
    else //For the begining?
        sum += log(1 - exp(-stepSize));

    double res = as * sum * ker * Int;

    return res;
}


////////////////////////////////////////////////
//Equation 86 (with eps reg) 
//BFKL with DGLAP
////////////////////////////////////////////////

double Solver::Kernel86(double l, double lp, double z)
{
    //Adding normal BFKL
    double res =  Kernel83(l, lp, z);

    res += DGLAPterm(l, lp, z);

    return res;
}

double Solver::Kernel86Diag(double l, double lp, double z)
{
    //return 0;
    return Kernel79Diag(l, lp, z);
}

//Off-diagonal z-diaginal kernel with F(k+q)
double Solver::Kernel86zDiag(double l, double lp, double x)
{
    return Kernel86zDiag(l, lp, x);
}

////////////////////////////////////////////////
//Equation 87 (with eps reg) 
//BFKL with DGLAP
////////////////////////////////////////////////


//Off-diagonal kernel with F(k+q)
double Solver::Kernel87(double l, double lp, double z)
{
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    if(z == 1) return 0;
    double angleMax1 = M_PI;
    // theta(k2/q2 - z/(1-z))
    double Cos = (lp*lp + l*l - l*l/(z/(1-z)) ) / (2*l*lp);
    if(Cos >= 1) return 0;
    else if(Cos > -1) angleMax1 = acos(Cos);

    //q < l
    if(lp > 2*l) return 0;
    double angleMax2 = acos(lp/(2.*l));

    double Int = 2 * IntegralCos(par, min(angleMax1, angleMax2));

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = acos((l*l + lp*lp - mu2) / (2*l*lp));
        Int -= 2 * IntegralCos(par, angleMin);
    }

    double res = M_PI * PggMod(z) * as * ker * Int; //Is Pi correct?

    //Adding normal BFKL
    res += Kernel83(l, lp, z);

    return res;
}

//Diagonal kernel with F(k)
double Solver::Kernel87Diag(double l, double lp, double z)
{
    if(z == 1) return 0;
    return Kernel79Diag(l, lp, z);
}


//Off-diagonal kernel with F(k+q)
double Solver::KernelSub87(double l, double lp, double z)
{
    double ker, par;
    tie(ker, par) = GetKerPar(l, lp);

    if(z == 1) return 0;
    // theta(k2/q2 - z/(1-z))
    double angleMax1 = GetAngle(l,lp, l*l*(1-z)/z);
    if(angleMax1 == 0) return 0;



    //q < l
    double angleMax2 = GetAngle(l,lp, l*l);
    if(angleMax2 == 0) return 0;

    double Int = 0;

    //q > mu
    // theta(q2 - mu2) term
    if(pow(lp-l,2) < mu2) { 
        double angleMin = GetAngle(l,lp, mu2);
        
        if(l == lp)
            Int += 2 * IntegralPhiDiag(l, lp, angleMin);
        else
            Int -= 2 * IntegralPhi(l, lp, angleMin);
    }

    if(l == lp)
        Int -= 2 * IntegralPhiDiag(l, lp, min(angleMax1, angleMax2));
    else
        Int += 2 * IntegralPhi(l, lp, min(angleMax1, angleMax2));

    double res = M_PI * PggMod(z) * as * Int; //Is Pi correct?

    //Adding normal BFKL
    res += KernelSub83(l, lp, z);

    return res;
}

//Diagonal kernel with F(k)
double Solver::KernelSub87Diag(double l, double lp, double z)
{
    return KernelSub79Diag(l, lp, z);
}





double Solver::Delta(double z, double k2, double q2)
{
    if(k2 <= (k2+q2)*z) return 1;
    if(k2 <= mu2) return 1;

    return exp( -as * log(k2/((k2+q2)*z)) * log(k2/mu2) );
}



////////////////////////////////////////////////
//Equation 88 (with eps reg) 
//resumed BFKL
////////////////////////////////////////////////

//Off-diagonal kernel with F(k+q)
double Solver::KernelSub88(double l, double lp, double z)
{
    if (z == 1) //Delta has no effect
        return KernelSub80(l, lp, z);


    //q > mu
    // theta(q2 - mu2) term
    double angleMin = 0;
    if(pow(lp-l,2) < mu2) { 
        angleMin = GetAngle(l,lp, mu2);
    }
    double angleMax = GetAngle(l,lp, l*l*(1-z)/z);


    //TODO integral
    
    auto fun = [&](double phi) {
        double q2=l*l + lp*lp - 2*l*lp*cos(phi);
        return lp*lp/(2*M_PI) * 1./q2 * Delta(z, l*l, q2);
    };

    auto funAlt = [&](double phiInv) {
        double phi = 1/phiInv;
        double q2=l*l + lp*lp - 2*l*lp*cos(phi);
        double res = lp*lp/(2*M_PI) * 1./q2 * Delta(z, l*l, q2);
        res *= phi*phi;
        return res;
    };
    auto funAlt2 = [&](double t) {
        double a = sqrt(pow(l-lp,2)/(l*lp));
        double phi = a * tan(a*t);
        double q2=l*l + lp*lp - 2*l*lp*cos(phi);
        double res = lp*lp/(2*M_PI) * 1./q2 * Delta(z, l*l, q2);
        res *= (phi*phi + a*a);
        return res;
    };

    auto funAlt3 = [&](double s) {
        double K = (l+lp)/abs(l-lp);
        double phi = 2*atan(tan(s)/K);

        double q2=l*l + lp*lp - 2*l*lp*cos(phi);
        double res = lp*lp/(2*M_PI) * 1./q2 * Delta(z, l*l, q2);
        res *= q2 *2./ abs(l*l - lp*lp);
        assert(res >= 0);
        return res;
    };

    double rat = (angleMax > angleMin) ? Delta(z, l*l, pow(l+lp,2) ) / Delta(z, l*l, pow(l-lp,2) ) : 1;
    double alpha = as * log(l*l/mu2);
    double kapa = 2*l*lp/ (2*l*l + lp*lp);

    auto funNest = [&](double s) {
        double K = (l+lp)/abs(l-lp);
        //double
        double phi = 2*atan(tan(s)/K);

        double q2=l*l + lp*lp - 2*l*lp*cos(phi);
        double res = lp*lp/(2*M_PI) * 1./q2 * Delta(z, l*l, q2);
        res *= q2 *2./ abs(l*l - lp*lp);
        assert(res >= 0);
        return res;
    };




    double Int = 0, err = 0;
    double Int2 = 0, err2 = 0;
    if(angleMax > angleMin) {
        //Int += 2*Integral61(fun, angleMin, angleMax, err);
        /*
        double err1=0, err2=0;
        double Int1 = Integral61(fun, angleMin, angleMax, err1);
        double Int2 = Integral61(funAlt, 1/angleMax, 1/angleMin, err2);
        cout << "Error  " <<  l<<" "<<lp<<" : "<< angleMin <<" "<< angleMax  <<" : "<<   Int1<<" "<<err1 << endl;
        cout << "Error2 " <<  l<<" "<<lp<<" : "<< angleMin <<" "<< angleMax  <<" : "<<   Int2<<" "<<err2 << endl;
        */
        if(l != lp) {
            double a = sqrt(pow(l-lp,2)/(l*lp));
            Int = 2*Integral61(funAlt2, 1/a * atan(angleMin/a), 1/a * atan(angleMax/a), err);
    
            Int2 = 2*Integral61(fun, angleMin, angleMax, err2);

            double K = (l+lp)/abs(l-lp);
            double Int3, err3=0;
            if(angleMax != M_PI) {
                double Min = atan(K*tan(angleMin/2));
                double Max = atan(K*tan(angleMax/2));
                assert(Max > Min);
                Int3 = 2*Integral61(funAlt3, atan(K*tan(angleMin/2)), atan(K*tan(angleMax/2)), err3);
                if( err3 / Int3 > 0.01) {
                cout << "Helenka K Kp " << l <<" "<< lp << endl;
                cout <<"Helenka " << Int <<" "<< Int2 << " "<<Int3 <<" : "
                                  << err <<" "<< err2 <<" "<< err3<< endl;
                cout << "Helenka Rat " << Delta(z, l*l, pow(l+lp,2)) / Delta(z, l*l, pow(l-lp,2))  << endl;
                }
            }
            //cout << "Error3 " <<  l<<" "<<lp<<" : "<< angleMin <<" "<< angleMax  <<" : "<<   Int3<<" "<<err3 << endl;

        }
        else {
            Int = 2*Integral61(funAlt, 1/angleMax, 1/angleMin, err);
            Int2 = 2*Integral61(fun, angleMin, angleMax, err2);

        }
        //cout << "Error " <<  l<<" "<<lp<<" : "<< angleMin <<" "<< angleMax  <<" : "<<   Int<<" "<<err << endl;
        if(Int > 1e-10 && err/Int > 1e-2) {
            cout << "Be Careful " <<l<<" "<<lp<<" : "<<angleMin<<" "<<angleMax<<" <> "<< Int <<" "<< err << " "<< err2<< endl;
            err = 0;
            //cout << "Check " << 2*Integral61(fun, angleMin, angleMax, err) <<" "<< err<< endl;
            //cout << "Check2 " << 2*Integral61(fun, 1/angleMax, 1/angleMin, err) <<" "<< err<< endl;
            cout << "Int2 " << Int2 << endl;
            double a = sqrt(pow(l-lp,2)/(l*lp));
            for(double x = angleMin; x < angleMax; x+=0.01)
                cout <<"list "<< x <<" "<<  fun(x) <<" "<< funAlt2(1/a *atan(x/a)) << endl;
        }
    }

    double maxTot = max(angleMin, angleMax);


    if(maxTot < M_PI) {
        if(l != lp)
            Int += 2*(IntegralPhi(l,lp, M_PI) - IntegralPhi(l,lp, maxTot));
        else {
            Int += 2*IntegralPhiDiag(l,lp, maxTot);
        }

    }



    double res = as * Int; //Is Pi correct?

    return res;
}

//The diagonal peace is missing here
