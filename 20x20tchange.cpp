// program that makes use of the IsingModel class defined in "IsingModel.hpp" and runs simulations for a lattice of size 20x20
// after compiling the files and running this simulation, the data can be visualized by running the python program 20x20.py

// we simulate the behaviour of a single instance of a 20x20 lattice, with 10^6 monte carlo cycles, for a changing temperature from 0 to 6, with dT=0.5

#include "IsingModel.hpp"

int main() {
    int L = 20; // size set to 20x20
    int cycles = 1000000; // 10^6 MC cycles
    IsingModel IsingModel(L, "20x20tchange", cycles); // defines environment, data file will have tha basename "20x20tchange"
    IsingModel.initialize();  // initializes spins in random directions
    for (double T = 0; T <= 6.0; T += 0.1){ // runs simulations and logs data for envitonments at different temperatures
        IsingModel.simulate_log(T, cycles);  // simualtes and logs our files for expectation values, energy and amount of spin-flips accepted
        IsingModel.reset();   // resets expectation values$
    } 
    return 0;
}