Authors: 
- [Synne Sigstad Jørgensen](https://github.uio.no/synnesjo)
- [Rebecca Nguyen](https://github.uio.no/rebeccng)

Build with
```sh
g++ main.cpp src/*.cpp -I include -o main.exe -larmadillo
```
Run as
```sh
./main <T> <L> <N_cycles_burn> <N_cycles> <output_filename>
```
The terminal line arguments are as follows:
- T: System temperature [ $J/k_B$ ].
- L: Lattice size (this code generate square lattices).
- N_cycles_burn: Number of Monte Carlo cycles included in burn-in time
- N_cycles: Number of Monte Carlo cycles.
- output_filename: Self explanatory. Make sure output filename is in format: textfiles/filename.txt

See the comments over the parts in the code where we print to file, as it is described there when to use which print-out statement.