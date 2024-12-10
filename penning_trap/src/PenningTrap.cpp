#include "../include/PenningTrap.hpp"
#include "../include/Particle.hpp"

#include <armadillo>
#include <vector>
#include <iostream>
#include <cmath>

double k_e = 1.38935333 * pow(10, 5);   // Coulomb constant 
double T = 9.64852558 * pow(10, 1);     // Magnetic field strenght
double V = 9.64852558 * pow(10, 7);     // Electric potential

double B0 = 9.65 * pow(10, 1);
double V0 = 2.41 * pow(10, 6);
double d = 500;
double V0_div_d_squared = 9.65;

//PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_in)
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
}

void PenningTrap::add_particle(Particle p_in)
{
    particle_.push_back(p_in);
}

arma::mat PenningTrap::external_E_field(arma::mat R)
{
    int N = particle_.size();
    arma::mat E_tot = arma::zeros(3, N);
    double C = V0_ / pow(d_, 2);

    for (int i = 0; i < N; i ++)
    {   
        arma::vec r = R.col(i);
        E_tot.col(i) = C * arma::vec{r.at(0), r.at(1), -2. * r.at(2)};
    }
    return E_tot;
}

arma::mat PenningTrap::external_B_field(arma::mat V)
{
    int N = particle_.size();
    arma::vec B = arma::vec{0., 0., B0_};
    arma::mat B_tot = arma::zeros(3, N);

    for (int i = 0; i < N; i ++)
    {
    double q = particle_[i].q();
    arma::vec v = V.col(i);
    B_tot.col(i) = arma::cross(q * v, B);
    }

    return B_tot;
}


arma::mat PenningTrap::total_force_external(arma::mat R, arma::mat V) //this works!!
{
    arma::mat E = PenningTrap::external_E_field(R);
    arma::mat B = PenningTrap::external_B_field(V);

    arma::mat F = E + B; 

    return F;
}

arma::vec PenningTrap::force_particle(arma::vec r_i, arma::vec r_j, double q_i, double q_j)
// Force on particle_i from particle_j
{
    arma::vec FP = k_e * q_j * (r_i - r_j) / ( pow(arma::norm(r_i - r_j), 3) );
    return FP;
}

arma::mat PenningTrap::total_force_particles(arma::mat R) // this function yields only zeros !!!
{
    int N = particle_.size();
    arma::mat F_interaction = arma::zeros(3, N);

    for (int i = 0; i < N; i ++)
    {
        double q_i = particle_[i].q();
        arma::vec r_i = R.col(i);

        for (int j = 0; j < N; j ++)
        {
            if (j == i)
            {
                //std::cout << F_tot << std::endl;
                // If particle i and particle j, there is no particle interaction
            }

            else
            {
                double q_j = particle_[j].q();
                arma::vec r_j = R.col(j);
                F_interaction.col(i) += PenningTrap::force_particle(r_i, r_j, q_i, q_j);
                
            }
        }

    }
    return F_interaction;
}

arma::mat PenningTrap::total_force(arma::mat R, arma::mat V, bool particle_interaction) // this function yields only zeros !!!
{
    int N = particle_.size();
    arma::mat F_tot = arma::mat(3,N);
    if (particle_interaction)
    {
        F_tot = PenningTrap::total_force_external(R, V) + total_force_particles(R);
    }
    else
    {
        F_tot = PenningTrap::total_force_external(R, V);
    }
    return F_tot;
}


// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, bool particle_interaction)
{
    int N = particle_.size();
    arma::mat R = arma::mat(3, N);
    arma::mat V = arma::mat(3, N);
    double m = particle_[0].m(); // assuming we are using same type of particles throughout, thus same mass

    for (int i = 0; i < N; i ++)
    {
        R.col(i) = particle_[i].r_;
        V.col(i) = particle_[i].v_;
    
        arma::mat a = total_force(R, V, particle_interaction) / m; // acceleration
        R += dt * V;
        V += dt * a;

        for (int i = 0; i < N; i ++)
        {
            particle_[i].r_ = R.col(i);
            particle_[i].v_ = V.col(i);

        }
    
    }
    
}

 // Evolve the system one time step (dt) using Runge-Kutta 4th order method
void PenningTrap::evolve_RK4(double dt, bool particle_interaction)
{
 
    int N = particle_.size();
    double m = particle_[0].m(); //means that we are assuming that the particle has the same mass throughout the simulation
    // possible to declare multiple variable in one line?? 

    arma::mat a_1; arma::mat a_2; arma::mat a_3; arma::mat a_4;
    
    //we start with empty matrices so we can fill them as we update them
    arma::mat R = arma::zeros(3, N);
    arma::mat V = arma::zeros(3, N);

    arma::mat kr_1 = arma::zeros(3, N);
    arma::mat kv_1 = arma::zeros(3, N);

    arma::mat kr_2 = arma::zeros(3, N);
    arma::mat kv_2 = arma::zeros(3, N);

    arma::mat kr_3 = arma::zeros(3, N);
    arma::mat kv_3 = arma::zeros(3, N);

    arma::mat kr_4 = arma::zeros(3, N);
    arma::mat kv_4 = arma::zeros(3, N);

    for (int i = 0; i < N; i++)
    {   
        R.col(i) = particle_[i].r_;
        V.col(i) = particle_[i].v_;  
    }

    a_1 = total_force(R, V, particle_interaction) / m;
    kr_1 = dt * V;
    kv_1 = dt * a_1; 


    a_2 = total_force(R + 0.5*kr_1, V + 0.5*kv_1, particle_interaction) / m;
    kr_2 = dt * (V + ((1/2)*kv_1));
    kv_2 = dt * a_2;


    a_3 = total_force(R + 0.5*kr_2, V + 0.5*kv_2, particle_interaction) / m;
    kr_3 = dt * (V + ((1/2)*kv_2)); 
    kv_3 = dt * a_3; 


    a_4 = total_force(R + kr_3, V + kv_3, particle_interaction) / m;
    kr_4 = dt * (V + kv_3);
    kv_4 = dt * a_4;

    //std::cout << (kr_1 + 2*kr_2 + 2*kr_3 + kr_4) << std::endl;
        
    R += 1/6. * (kr_1 + 2*kr_2 + 2*kr_3 + kr_4);
    V += 1/6. * (kv_1 + 2*kv_2 + 2*kv_3 + kv_4);

    //std::cout << R << std::endl;
    //std::cout << V << std::endl;

    for (int i = 0; i < N; i++)
    {
        particle_[i].r_ = R.col(i);
        particle_[i].v_ = V.col(i);
    }
}