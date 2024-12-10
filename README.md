Authors: 
- [Synne Sigstad Jørgensen](https://github.uio.no/synnesjo)
- [Rebecca Nguyen](https://github.uio.no/rebeccng)

A Penning trap utlize a static configuration of quadrupole electric and homogeneous magnetic fields to confine charge particles. Such instrument is useful in experiments, including CERN. 

We have two classes: One called ```Particle.cpp``` which store a single particle with given properties and ```PenningTrap.cpp``` which sets up the Penning Trap itself

## Simulations

The file ```main.cpp``` allows you to simualate the Penning trap with one or two Ca+ ions with/without interactions. In this file you can simulate particle 1 by itself or particle 1 and 2 together. 
- To add/remove particle 2 comment in/out relevant code block, see code.
- Turn on/off particle interaction by setting ```particle_interaction``` to true/false.
- Under the headings 'SETTING FILENAMES FOR SINGLE PARTICLE SIMULATIONS' and 'SETTING FILENAMES FOR TWO PARTICLE SIMULATIONS': there are multiple vectors filled with strings corresponding text filenames used for plotting. Comment out the relevant lines for the simulations you are running.

To compare RK4 with analytical solution and Euler:
1. Run ```euler_and_analytical.cpp```, choose number time steps
2. Run ```main.cpp``` with particle 1 and same number of time steps as point 1.
NOTE: Under the headings 'FORWARD EULER TEXT FILES' and 'ANALYTICAL TEXT FILES': there are multiple vectors filled with strings corresponding text filenames used for plotting. Comment out the relevant lines for the simulations you are running.


Build as follows for Windows
```sh
g++ filename.cpp src/*.cpp -I include -o filename_executable -larmadillo
```
and run the executable
```sh
./filename_executable
```

If there is an error message requiring you to enable a newer version of the compiler, build as follows
```sh
g++ -std=c++17 filename.cpp src/*.cpp -I include -o filename_executable -larmadillo
```

