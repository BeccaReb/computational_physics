#include <vector>
#include <iostream>
#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>
#include <typeinfo>

#include "../include/PenningTrap.hpp"
#include "../include/Particle.hpp"
#include "../include/Analytical.hpp"

int main()
{
    // constants
    double k_e = 1.38935333 * pow(10, 5);   // Coulomb constant 
    double T = 9.64852558 * pow(10, 1);     // Magnetic field strenght
    double V = 9.64852558 * pow(10, 7);     // Electric potential
    
    bool particle_interaction = false;

    // Penning trap configuration
    double B0 = 9.65 * pow(10, 1);
    double V0 = 2.41 * pow(10, 6);
    double d = 500;
    double V0_div_d_squared = 9.65;

    // particle properties
    double q = 1.; // nothing specified in project description?
    double m = 40.078; // mass of single-ionised Calcium ion [u] 
    arma::vec r1 = arma::vec{20., 0., 20.}; // [mu*m]
    arma::vec r2 = arma::vec{25., 25., 0.}; // [mu*m]
    arma::vec v1 = arma::vec{0., 25., 0.}; // [(mu*m) / (mu*s)]
    arma::vec v2 = arma::vec{0., 40., 5.}; // [(mu*m) / (mu*s)]

    double t_end = 50; // total simualtion time [mu*s]
    double n_steps = 4000;
    double dt = t_end / n_steps ;  
    arma::vec t = arma::linspace(0, t_end, n_steps);
    int width = 12;
    int prec  = 5;

    
    // setting up Penning trap
    PenningTrap penningtrap = PenningTrap(B0, V0, d);
    
    //Particle 1:
    Particle Ca1 = Particle(q, m, r1, v1);
    penningtrap.add_particle(Ca1);
    
    
    //Particle 2:
    //the next two lines of code must be commented in if you want to add a second particle to the Penning trap.
    Particle Ca2 = Particle(q, m, r2, v2);
    penningtrap.add_particle(Ca2);
    
    
    
    /*SETTING FILENAMES FOR SINGLE PARTICLE SIMULATIONS:
    ENABLE WHEN WORKING WITH PARTICLE 1:*/
    //---------------------------------------------------------------------------------------------
    //std::vector<std::string> filenames = {"textfiles/single_particle_1_4000_rk.txt"};
    //std::vector<std::string> filenames = {"textfiles/single_particle_1_8000_rk.txt"};
    //std::vector<std::string> filenames = {"textfiles/single_particle_1_32000_rk.txt"};
    //std::vector<std::string> filenames = {"textfiles/single_particle_1_32000_rk.txt"};
    //---------------------------------------------------------------------------------------------

    
    
    /*SETTING FILENAMES FOR TWO PARTICLE SIMULATIONS:
    ENABLE WHEN WORKING WITH NO INTERACTIONS:*/
    //---------------------------------------------------------------------------------------------
    //std::vector<std::string> filenames = {"textfiles/particle_1_4000_no_interaction_rk.txt", "textfiles/particle_2_4000_no_interaction_rk.txt"};
    //---------------------------------------------------------------------------------------------

    /*ENABLE WHEN WORKING WITH INTERACTIONS:*/
    /*------------------------------------------------------------------------------------------------------------------------------------------------*/
    //std::vector<std::string> filenames = {"textfiles/particle_1_4000_with_interaction_rk.txt", "textfiles/particle_2_4000_with_interaction_rk.txt"};
    /*------------------------------------------------------------------------------------------------------------------------------------------------*/


    arma::vec initial_position = arma::vec{20., 0., 20., 25., 25., 0.};
    arma::vec initial_velocity = arma::vec{0., 25., 0., 0., 40., 5.};
    std::vector<arma::Mat<double>> info;

    for (int k = 0; k < 2; k++)
    {
        arma::Mat<double> zero_matrix(n_steps, 7, arma::fill::zeros);
        info.push_back(zero_matrix);
    }

   for (int i = 1; i < n_steps; i++)
    {   
        penningtrap.evolve_RK4(dt, particle_interaction);
        
        for (int j = 0; j < penningtrap.particle_.size(); j++)
        {   
            info.at(j).row(i) = arma::vec({t.at(i), penningtrap.particle_.at(j).r_.at(0), penningtrap.particle_.at(j).r_.at(1), penningtrap.particle_.at(j).r_.at(2), 
            penningtrap.particle_.at(j).v_.at(0), penningtrap.particle_.at(j).v_.at(1), penningtrap.particle_.at(j).v_.at(2)}).t();
        }


    }

    //printing to file
    for(int l = 0; l < penningtrap.particle_.size(); l++)
    {
        std::ofstream ofile;
        ofile.open(filenames.at(l));
        
        ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << std::setw(width) << std::setprecision(prec) 
        << std::scientific << t.at(0) << "," << initial_position.at((3*l)) << "," << initial_position.at((3*l)+1) << "," 
        << initial_position.at((3*l)+2) << "," << initial_velocity.at((3*l)) << "," << initial_velocity.at((3*l)+1) << "," << initial_velocity.at((3*l)+2) << std::endl;

        for(int m = 1; m < n_steps; m++)
        {
                ofile<< std::setw(width) << std::setprecision(prec) << std::scientific << std::setw(width) << std::setprecision(prec) << std::scientific
                << info[l].row(m)[0] << "," << info[l].row(m)[1] << "," << info[l].row(m)[2] << "," << info[l].row(m)[3] <<  ","
                << info[l].row(m)[4] << "," << info[l].row(m)[5] << "," << info[l].row(m)[6] << std::endl;
        }
        ofile.close();
    }

    return 0;
}