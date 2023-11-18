// program that makes use of the IsingModel class defined in "IsingModel.hpp" and runs simulations for lattices of varying sizes
// after compiling the files and running this simulation, the data can be visualized by running the python program sizevariation.py

// we simulate the behaviour for environments of lattices of sizes in the interval 40-120 with dL 20, and varying temperatures from 2.0 to 2.5 with dT=0.01
#include "IsingModel.hpp"

int main() {
    int cycles = 100000; // 10^5 MC cycles
    #pragma omp parallel for //parallelizes the loop
    for (int L = 40; L <= 120; L += 20) { // size-loop
        IsingModel IsingModel(L, "mcc5zoom" + to_string(L), cycles); // defines environment, data file will have tha basename "mcc5"
        IsingModel.initialize(); // spins initialized in raodnom directions, up or down
        for (double T = 2.2; T <= 2.4; T += 0.005) { // temperature loop for each lattice size
            IsingModel.simulate_log(T, cycles); // runs simulations and logs data
            IsingModel.reset(); // resets expectation values  
        }  
    }
    return 0;
}
