// program that makes use of the IsingModel class defined in "IsingModel.hpp" and runs simulations for lattices of varying sizes, over varying temperatures without logging data 
// the purpose of this program is to log how long the simulation takes on average, for each size lattice when the code has been parallized by using #pragma omp parallel for 
// each simuilation is ran for a temperature intervall, the spesifics of this is not as important as we wish to show the difference in efficiency when parallelizing code and when not, so as long as this is done in the same intervall it does not matter
// after compiling the files and running this simulation, the data can be visualized by running the python program timer.py

#include "IsingModel.hpp"

int main() {
	int MCCs = 100000; // want to run the simulation for 10^5 monte carlo cycles
	string filename = "tmr5par"; // base name of the file containing time-data
    //defines the intervall of lattice sizes we wish to perform time tests for
    int Lstart = 10; 
	int Lend = 100;
	int dL = 5;
    // defines variables eneded needed simulation
    int Nthreads = 4;
	int Ntests = 3*Nthreads;
    int NL = (Lend-Lstart)/dL;

	ofstream outfile("data/"+filename+".txt", ios_base::app); // opens the file we log data in 
	for(int i = 0; i <= NL; i++) { // loops through a temperature intervall that we define
        int L = Lstart + i*dL; // defines lattice size L in this point
        outfile << L << " "; // writes L into the file
        // we start a timer for the further parallelized code
        omp_set_num_threads(Nthreads);
        double start = omp_get_wtime();
        #pragma omp parallel for //parallelizes the loop
        for(int N = 0; N < Ntests; N++) {
            // this loop runs our simulation 
            double T = 1 + N*0.1; 
            IsingModel problem = IsingModel(L, filename, MCCs); //creates environment 
            problem.initialize(); //initializes spins in randomly generated direction up or down
            problem.simulate(T, MCCs); // simulates the environment's development for the 10^5 monte carlo cycles, without logging data
            problem.reset();	 // resets expectation values
        }
        double end = omp_get_wtime(); // stops timer
        double totalTime = (end-start); //calculates the time intervall used by machine to finish this simulation 
        totalTime = double (totalTime/Ntests); //calculates the time intervall used by machine to finish this simulation 
        outfile << totalTime << " " <<endl; // writes the time it took for the simulation to be ran for the L size lattice
	}
	outfile.close(); // closes file
}