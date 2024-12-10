#ifndef __lattice_hpp__  
#define __lattice_hpp__
#include <vector>
#include <iostream>
#include <armadillo>
#include <numeric>
#include <sstream>
#include <string>
#include <iomanip>
#include <typeinfo>

// initialize L x L lattice for Ising-Lenz model
arma::mat lattice(int L, bool unordered);

// calculate the initial total energy of the lattice
double initial_energy(arma::mat S, int L);

// calculate the initial total magnetization of the lattice
double initial_magnetization(arma::mat S, int L);

// the five possible values for dE
std::vector<double> Boltzmann(double T);

#endif 