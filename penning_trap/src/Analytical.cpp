#include "../include/Analytical.hpp"

#include <vector>
#include <iostream>
#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>

arma::vec analytical(arma::vec r_input, arma::vec v_input, double q, double B_0, double V_0, double m, double d, double t)
{
    arma::vec r = arma::vec(3);

    double x_0 = r_input.at(0);
    double z_0 = r_input.at(2);
    double v_0 = v_input.at(1);
    

    double omega_0 = (q*B_0) / m;
    double omega_z = std::sqrt( (2*q*V_0) / (m*pow(d, 2)));
    double omega_plus = (omega_0 + std::sqrt( pow(omega_0, 2) - 2*pow(omega_z, 2) )) / (2.);
    double omega_minus = (omega_0 - std::sqrt( pow(omega_0, 2) - 2*pow(omega_z, 2) )) / (2.);

    double A_plus = (v_0 + omega_minus * x_0) / (omega_minus - omega_plus);
    double A_minus = -(v_0 + omega_plus * x_0) / (omega_minus - omega_plus);

    r.at(0) = A_plus*std::cos(omega_plus*t) + A_minus*std::cos(omega_minus*t);
    r.at(1) = A_plus*std::sin(omega_plus*t) + A_minus*std::sin(omega_minus*t);
    r.at(2) = z_0 * std::cos(omega_z * t);

    return r;
}