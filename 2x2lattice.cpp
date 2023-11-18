// program that makes use of the IsingModel class defined in "IsingModel.hpp" and runs simulations for a lattice of size 2x2
// after compiling the files and running this simulation, the data can be visualized by running the python program 2x2.py

#include "IsingModel.hpp"

int main() {
    int L = 2; // lattice size 2x2
    for (int i = 1; i <= 5; ++i) { // monte carlo cycles 10, 100, 1000, 10000, 100000
        IsingModel isingModel(L, "2x2lattice" + to_string(i), pow(10, i)); // creates a class that will write data into a .txt file with stem "2x2lattice" 
        isingModel.initialize(); // initializes spins in random directions
        // we do not need to log any other information than the expected values for this simulation
        for (double T = 0.5; T <= 4.0; T += 0.1) { // runs simulation for temperature in intervall 0.5 to 4.0 in intervalls of dT=0.1
            isingModel.simulate_log(T, pow(10, i));// logs expectation values, flipped spins and energy
            isingModel.reset(); // resets expectation values
        }
    }
    return 0;
}