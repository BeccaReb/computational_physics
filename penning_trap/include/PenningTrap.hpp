#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <vector>
#include "../include/Particle.hpp"

class PenningTrap
{
    
public:
    double B0_; // Magnetic field strength
    double V0_; // Applied potential
    double d_;  // Characteristic dimension
    //arma::vec<Particle>particle_;
    
    std::vector<Particle>particle_; // Contains all particle objects (not sure if this is right way to declare)
    //this should also be an arma vec for consistency, but not sure how to do this

    // Constructor
    //PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_in);
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // External electric field at point r=(x,y,z)
    arma::mat external_E_field(arma::mat R);

    // External magnetic field at point r=(x,y,z)
    arma::mat external_B_field(arma::mat V);

    // Force on particle_i from particle_j
    arma::vec force_particle(arma::vec r_i, arma::vec r_j, double q_i, double q_j);

    // The total force on particle_i from the external fields
    arma::mat total_force_external(arma::mat R, arma::mat V);

    // The total force on particle_i from the other particles
    arma::mat total_force_particles(arma::mat R);

    // The total force on particle_i from both external fields and other particles
    arma::mat total_force(arma::mat R, arma::mat V, bool particle_interaction);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt, bool particle_interaction);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt, bool particle_interaction);

};



#endif