#include <vector>
#include <iostream>
#include <armadillo>
#include <numeric>
#include <sstream>
#include <string>
#include <iomanip>
#include <typeinfo>

#include "../include/lattice.hpp"

// initialize L x L lattice for Ising-Lenz model
arma::mat lattice(int L, bool unordered) 
{
    // the default is ordered initial state (all spin up -> 1)
    arma::mat M = arma::mat(L, L).fill(1); 

    // if unordered is true, then we randomize the initial state
    arma::mat dist = arma::mat(L, L).randn(); // L x L matrix filled with random numbers from U(0, 1)
    if (unordered)
    {
        for (int i = 0; i < L; i++)
        {   
            for(int j = 0; j < L; j++)
            {
                double y = dist(i, j); // spin state from dist
                if (y <= 0) // y is less than or equal to 0, then flip the spin. Else, keep it the same.
                {
                    M(i, j) = M(i, j) * (-1);
                }
            }        
        }    
    }
    
    return M;
}


// calculate the initial total energy of the lattice
double initial_energy(arma::mat S, int L)
{
    double E = 0;
    for (int ix = 0; ix < L; ix++)
    {
        for (int iy = 0; iy < L; iy++)
        {
            // NO DOUBLE COUNTING
            E -= S(ix, iy) * 
                (S(ix, ( (iy - 1) < 0 ? (L-1) : (iy - 1) )) +  // up
                S(( (ix + 1) == (L) ? 0 : (ix + 1) ), iy) ); // right
        }
    }
    return E;
}

// calculate the initial total magnetization of the lattice
double initial_magnetization(arma::mat S, int L)
{
    double M = 0;
    for (int ix = 0; ix < L; ix++)
    {
        for (int iy = 0; iy < L; iy++)
        {
            M += S(ix, iy);
        }
    }
    return M;
}


// the five possible values for dE
std::vector<double> Boltzmann(double T)
{
    double beta = 1 / T;
    std::vector<double> dE(5);
    std::vector<double> boltzmann(5);
    for(int i=0; i<5; i++){
        dE.at(i) = (4 * (i - 2));
        boltzmann.at(i) = exp(-dE.at(i)*beta);
    }
    return boltzmann;
}
