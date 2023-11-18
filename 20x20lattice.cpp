// program that makes use of the IsingModel class defined in "IsingModel.hpp" and runs simulations for a lattice of size 20x20
// after compiling the files and running this simulation, the data can be visualized by running the python program 20x20.py

// we simulate the behaviour of four different environments of the same size, for the temperature T=1 and T=2.4 with spins initialized in the same, and random directions
#include "IsingModel.hpp"

int main() {
    int L = 20; // size set to 20x20
    int cycles = 1000000; // 10^6 MC cycles

    //for T=1
    double T= 1;
    IsingModel isingt1(L, "20x20ordered", cycles); // isingmodel class environment is establised with a file name base "20x20ordered"
    isingt1.initialize_ordered(); // spins initialized in the same direction
    isingt1.simulate_log(T, cycles); //runs simulations and logs data
    isingt1.reset(); // resets expectation values

    IsingModel isingt1up(L, "20x20random", cycles); // isingmodel class environment is establised with a file name base 20x20random"
    isingt1up.initialize(); // spins initialized in raodnom directions, up or down
    isingt1up.simulate_log(T, cycles); //runs simulations and logs data
    isingt1up.reset();  // resets expectation values

    //process is repeated for T=2.4
    T= 2.4;
    IsingModel isingt2(L, "20x20ordered", cycles);
    isingt2.initialize_ordered();
    isingt2.simulate_log(T, cycles);
    isingt2.reset(); 

    IsingModel isingt2up(L, "20x20random", cycles);
    isingt2up.initialize();
    isingt2up.simulate_log(T, cycles);
    isingt2up.reset(); 
    return 0;
}