#include <vector>
#include <iostream>
#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>
#include <typeinfo>


double Z(double T, double J)
{
    double beta = 1 / T;
    return 12 + 4 * std::cosh(8*beta*J);
}

double expectation_value_e(double T, double J, double L, double N)
{
    double beta = 1 / T;
    double e = - (8*J)/Z(T, J) * std::sinh(8*beta*J);
    return e;
}

double expectation_value_ee(double T, double J, double L, double N)
{
    double beta = 1 / T;
    double ee = 16*J*J/Z(T, J) * std::cosh(8*beta*J);
    return ee;
}

double expectation_value_m(double T, double J, double L, double N)
{
    double beta = 1 / T;
    double m = (4 + 2*std::exp(8*beta*J))/Z(T, J);
    return m;
}

double expectation_value_mm(double T, double J, double L, double N)
{
    double beta = 1 / T;
    double mm = 2/Z(T, J) * (1 + std::exp(8*beta*J));
    return mm;
}

double heat_capacity(double T, double J, double L, double N)
{
    double e = expectation_value_e(T, J, L, N);
    double ee = expectation_value_ee(T, J, L, N);
    double cv = N * (1/(T*T)) * (ee - e*e); 
    return cv;
}

double susceptibility(double T, double J, double L, double N)
{
    double m = expectation_value_m(T, J, L, N);
    double mm = expectation_value_mm(T, J, L, N);
    double chi = N  * (mm - m*m)/T;
    return chi;
}

int main()
{
    double J = 1.;
    double T = 1.;
    double L = 2.;
    double N = L*L;

    double e = expectation_value_e(T, J, L, N);
    double ee = expectation_value_ee(T, J, L, N);
    double m = expectation_value_m(T, J, L, N);
    double mm = expectation_value_mm(T, J, L, N);
    double cv = heat_capacity(T, J, L, N);
    double chi = susceptibility(T, J, L, N);

    std::cout << "<e>: " << e << std::endl;
    std::cout << "<|m|>: " << m << std::endl;
    std::cout << "C_v/N: " << cv << std::endl;
    std::cout << "chi/N: " << chi << std::endl;

    return 0;
}