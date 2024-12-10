#include <vector>
#include <iostream>
#include <armadillo>
using namespace arma;
#include <numeric>
#include <sstream>
#include <string>
#include <iomanip>
#include <typeinfo>
#include <cassert>

#include "../include/lattice.hpp"

// Markov Chain Monte Carlo
void MCMC(arma::mat S, int L, int N_cycles_burn, int N_cycles, double& E, double& M, double T, std::string filename)
{
    std::mt19937 generator;
    generator.seed(std::clock());

    // five possible values for dE are stored in dE_vec
    std::vector<double> boltzmann_vec = Boltzmann(T);
    std::vector<double> energy_vec_cycles;
    std::vector<double> energy_vec_cycles_samples;
    std::vector<double> energy_vec_squared;
    std::vector<double> magnetic_vec;
    std::vector<double> magnetic_vec_squared;
    std::uniform_int_distribution<int> random(0, L - 1);
    double boltzmann_factor;
    double A;
    int dE;
    int samples = L*L;
    int ix = random(generator);
    int iy = random(generator);


    /* --------------------------------------------------------------------------------------------------------*/
    /*first we run the number of cycles we have defined as our burn in, without saving the results*/
    for (int i = 0; i < N_cycles_burn; i++)
    {
        for (int N = 0; N < samples; N++)
        {
            // sample a random spin (i, j) in lattice 
            ix = random(generator);
            iy = random(generator);
            
            std::uniform_real_distribution<double> r(0.0, 1.0); // acceptence probability

            //std::cout << S << std::endl;

            // dE if spin is flipped
            dE = 2*S(ix, iy) *                                             // sample
                     ((S(( (ix + 1) == (L) ? 0 : (ix + 1) ), iy)) +            // right
                     (S(ix, ( (iy + 1) == (L) ? 0 : (iy + 1) ))) +             // down
                     (S(( (ix - 1) < 0 ? (L-1) : (ix - 1) ), iy)) +            // left
                     (S(ix, ( (iy - 1) < 0 ? (L-1) : (iy - 1) ))));            // up
            
            //std::cout << "int dE is: " << dE << std::endl;
            // accept or reject flip
            boltzmann_factor = boltzmann_vec.at((dE / 4) + 2);
            A = std::min(1.0, boltzmann_factor);

            if (r(generator) < A)  //også godkjenne lavere energi, må fikses
            {  
                S(ix, iy) *= (-1); // flip spin
                E += dE;           // update energy
                M += 2*S(ix, iy);  // update magnetization
            }
        }  
          
    }

    /*----------------------------------------------------------------------------------------------------------*/
    /*here we start saving results, as we are past the burn-in time*/
    for (int i = 0; i < N_cycles; i++)
    {
        for (int N = 0; N < samples; N++)
        {
            // sample a random spin (i, j) in lattice 
            ix = random(generator);
            iy = random(generator);
            
            std::uniform_real_distribution<double> r(0.0, 1.0); // acceptence probability

            //std::cout << S << std::endl;

            // dE if spin is flipped
            dE = 2*S(ix, iy) *                                                 // sample
                     ((S(( (ix + 1) == (L) ? 0 : (ix + 1) ), iy)) +            // right
                     (S(ix, ( (iy + 1) == (L) ? 0 : (iy + 1) ))) +             // down
                     (S(( (ix - 1) < 0 ? (L-1) : (ix - 1) ), iy)) +            // left
                     (S(ix, ( (iy - 1) < 0 ? (L-1) : (iy - 1) ))));            // up
            
            //std::cout << "int dE is: " << dE << std::endl;
            // accept or reject flip
            boltzmann_factor = boltzmann_vec.at((dE / 4) + 2);
            A = std::min(1.0, boltzmann_factor);

            if (r(generator) < A)  //også godkjenne lavere energi, må fikses
            {  
                S(ix, iy) *= (-1); // flip spin
                E += dE;           // update energy
                M += 2*S(ix, iy);  // update magnetization
            }
            energy_vec_cycles_samples.push_back(E);
        }   

    energy_vec_cycles.push_back(E);
    energy_vec_squared.push_back(E*E);
    magnetic_vec.push_back(std::abs(M));
    magnetic_vec_squared.push_back(std::abs(M)*std::abs(M));
          
    }

    std::ofstream ofile;

    /*Use this when printing the results for the end states \epsilon*/

    // // print to file the values from each sample from each cycle
    // ofile.open(filename, std::ios_base::app); // append instead of overwrite
    // for (int i = 0; i < N_cycles; i++)
    // {    
    //     ofile << energy_vec_cycles.at(i)/(L*L) << std::endl;
    // }
    // ofile.close();


    /*-----------------------------------------------------------------------------------------------------------------------*/
    /*use this when saving the average of the individual samples to find the results for the end states <\epsilon>*/

    // std::vector<int> sums;
    // for (int i = 0; i < energy_vec_cycles_samples.size(); i += 400) {
    //     // Sum the next 100 elements
    //     int sum = std::accumulate(energy_vec_cycles_samples.begin() + i, energy_vec_cycles_samples.begin() + i + 400, 0);
    //     sums.push_back(sum);
    // }

    
    // //print to file the values from each cycle
    // ofile.open(filename, std::ios_base::app); // append instead of overwrite
    // for (int i = 0; i < N_cycles; i++)
    // {   
    //     ofile << sums.at(i)/(pow(L, 4)) << std::endl;
    // }
    // ofile.close();
    /*-----------------------------------------------------------------------------------------------------------------------*/




    double expected_energy_per_spin = arma::mean(arma::vec(energy_vec_cycles))/(pow(L,2));
    double expected_energy_per_spin_squared = arma::mean(arma::vec(energy_vec_squared))/(pow(L,4));
    double heat_capacity = (L*L)*(expected_energy_per_spin_squared - pow(expected_energy_per_spin, 2));

    
    double expected_magnetisation_per_spin = arma::mean(arma::vec(magnetic_vec))/(pow(L,2)); 
    double expected_magnetisation_per_spin_squared = arma::mean(arma::vec(magnetic_vec_squared))/(pow(L,4));
    double susceptibility = (L*L)*(expected_magnetisation_per_spin_squared - pow(expected_magnetisation_per_spin, 2)); 

    /*Use this to print a set of the quantities we can calculate*/

    // ofile.open(filename, std::ios_base::app); // append instead of overwrite
    // ofile << T << "," << expected_energy_per_spin << "," << heat_capacity << "," << expected_magnetisation_per_spin << "," << susceptibility << std::endl;
    // ofile.close();

}

int main(int argc, const char* argv[])
{   
    assert(argc == 6);                 // Check number of command line arguments
    double T = atof(argv[1]);          // temperature [J/k_B]
    int L = atoi(argv[2]);             // lattice size
    int N_cycles_burn = atoi(argv[3]); // number of burn-in MC cycles 
    int N_cycles = atoi(argv[4]);      // number of MC cycles
    std::string filename = argv[5];    // output filename

    std::cout << "Your input arguments are: " << std::endl;
    std::cout << "Temperature: " << T << "J/k_B" << std::endl;
    std::cout << "Lattice size: " << L << std::endl;
    std::cout << "Number of burn-in MC cycles: " << N_cycles_burn << std::endl;
    std::cout << "Number of MC cycles: " << N_cycles << std::endl;

    // initialize matrix
    int N = L*L;
    arma::mat S = lattice(L, true);

    double E = initial_energy(S, L);
    double M = initial_magnetization(S, L);

    MCMC(S, L, N_cycles_burn, N_cycles, E, M, T, filename);

    return 0;
} 