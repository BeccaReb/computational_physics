#include "../include/Particle.hpp"

#include <armadillo>
#include <vector>

// Constructor
Particle::Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
    q_ = q_in;
    m_ = m_in;
    r_ = r_in;
    v_ = v_in;
}

// Method that returns charge
double Particle::q()
{
    return q_;
}

// Method that returns mass
double Particle::m()
{
    return m_;
}

// Method that returns position
arma::vec Particle::r()
{
    return r_;
}

// Method that returns velocity
arma::vec Particle::v()
{
    return v_;
}