#include <iostream>

using namespace std;

//Integrals over phi/2pi
// 1/(a-b*cos(phi))
double intOverPow1(double a, double b)
{
    return 1/sqrt(a*a-b*b);
}

//Integrals over phi/2pi
// 1/(a-b*cos(phi))^2
double intOverPow2(double a, double b)
{
    return a/pow(a*a-b*b,3./2);
}

//Integrals over phi/2pi
// cos(phi)/(a-b*cos(phi))
double intCosOverPow1(double a, double b)
{
    //return a/b*1./sqrt(a*a-b*b) - 1/b;
    return 1/b* (a/sqrt(a*a-b*b) - 1);
}

//Integrals over phi/2pi
// cos(phi)/(a-b*cos(phi))^2
double intCosOverPow2(double a, double b)
{
    return b/pow(a*a-b*b,3./2);
}


//function of k
pair<double,double> integrand(double z, double mq2, double k2, double p2)
{
    double c2  = z*(1-z)*Q2 + mq2;
    double D1  = k2 + c2;
    double d2  = k2 + c2 + p2;
    double f   = z*z + (1-z)*(1-z);

    double D2a = sqrt(pow(d2,2) - 4*k2*p2);

    //Term (1/D1 - 1/D2)^2
    double r1 = 1/(D1*D1) - 2/(D1*D2a) + d2/pow(D2a,3);

    //Term (2*k*p/D1/D2 -2*k*p/D2^2 + k^2/D2^2) 
    double r2 = 1/D1*(d2/D2a - 1) - 4*k2*p2/pow(D2a,3) + k2*d2/pow(D2a,3);

    //Transverse part
    double FT = (f*k2 + mq2)*r1 + f*r2;

    //Longitudinal part
    double FL = 4*Q2*pow(z*(1-z),2)*r1;

    //Overall constant
    double fact = as*eq2*Q2/(4*M_PI) * k2;

    return make_pair(fact*FT, fact*FL);
}

int main()
{


    return 0;
}
