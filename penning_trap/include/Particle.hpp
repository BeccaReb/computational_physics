#ifndef __Particle_hpp__  
#define __Particle_hpp__

#include <armadillo>
#include <vector>

class Particle
{
protected:

public:
    double q_, m_;  // Charge q and mass m
    arma::vec r_;   // Position r 
    arma::vec v_;   // Velocity v

    // Constructor
    Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in);
    
    // Method that returns charge
    double q();

    // Method that returns mass
    double m();

    // Method that returns position
    arma::vec r();

    // Method that returns velocity
    arma::vec v();

};



#endif