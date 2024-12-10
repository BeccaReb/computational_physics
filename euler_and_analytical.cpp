
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


    // setting up Penning trap
    Particle Ca1 = Particle(q, m, r1, v1);
    //Particle Ca2 = Particle(q, m, r2, v2);
    PenningTrap penningtrap = PenningTrap(B0, V0, d);
    penningtrap.add_particle(Ca1);
    //penningtrap.add_particle(Ca2);


    double t_end = 50; // total simualtion time [mu*s]
    double n_steps = 32000;
    double dt = t_end / n_steps ;  
    arma::vec t = arma::linspace(0, t_end, n_steps);
    int width = 12;
    int prec  = 5;


 
    std::ofstream ofile; 

    
    /*FORWARD EULER TEXT FILES  
    ---------------------------------------------------------------------------------------------*/
    //std::string filenames_euler = "textfiles/single_particle_1_4000_euler.txt";
    //std::string filenames_euler = "textfiles/single_particle_1_8000_euler.txt";
    //std::string filenames_euler = "textfiles/single_particle_1_16000_euler.txt";
    //std::string filenames_euler = "textfiles/single_particle_1_32000_euler.txt";
    /*---------------------------------------------------------------------------------------------*/


    ofile.open(filenames_euler); //here you have to make sure the number in the file name corresponds to the number of time steps you have chosen!
        

    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << t.at(0) << "," 
            << std::setw(width) << std::setprecision(prec) << std::scientific << r1.at(0) << "," 
            << std::setw(width) << std::setprecision(prec) << std::scientific << r1.at(1) << "," 
            << std::setw(width) << std::setprecision(prec) << std::scientific << r1.at(2) <<  "\n";

        
    for (int i = 1; i < n_steps; i++)
    {
        penningtrap.evolve_forward_Euler(dt, particle_interaction);
        

        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << t[i] << "," << std::setw(width) << std::setprecision(prec) 
            << std::scientific << penningtrap.particle_[0].r_[0] << "," << std::setw(width) << std::setprecision(prec) << std::scientific
            << penningtrap.particle_[0].r_[1]  << "," << std::setw(width) << std::setprecision(prec) << std::scientific  
            << penningtrap.particle_[0].r_[2] <<  "\n";    
        }
    ofile.close();

    

    
    /*ANALYTICAL TEXTFILES
    ---------------------------------------------------------------------------------------------*/
    //std::string filenames_analytical = "textfiles/single_particle_1_4000_analytical.txt";
    //std::string filenames_analytical = "textfiles/single_particle_1_8000_analytical.txt";
    //std::string filenames_analytical = "textfiles/single_particle_1_16000_analytical.txt";
    //std::string filenames_analytical = "textfiles/single_particle_1_32000_analytical.txt";
    //---------------------------------------------------------------------------------------------


    ofile.open(filenames_analytical);

    for (int i = 0; i < n_steps; i++)
    {
        // analytical
        arma::vec r_analytical = analytical(r1, v1, q, B0, V0, m, d, t.at(i));
        
        
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << t[i] << "," << std::setw(width) << std::setprecision(prec) 
            << std::scientific << r_analytical[0] << "," << std::setw(width) << std::setprecision(prec) << std::scientific
            << r_analytical[1]  << "," << std::setw(width) << std::setprecision(prec) << std::scientific  
            << r_analytical[2] <<  "\n";    
        }
    
    ofile.close();

    return 0;
}
