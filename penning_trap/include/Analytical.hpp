#ifndef __Analytical_hpp__  
#define __Analytical_hpp__

#include <vector>
#include <iostream>
#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>

arma::vec analytical(arma::vec r_input, arma::vec v_input, double q, double B_0, double V_0, double m, double d, double t);

#endif 