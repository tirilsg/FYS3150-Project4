#include "IsingModel.hpp" // calls the file containing the structure of the IsingModel class


IsingModel::IsingModel(int size, string save_name, int amount_of_data) : L(size), filename(save_name), data_lines(amount_of_data), spins(size, size), expect(5, fill::zeros), boltzmann(17, fill::zeros), generator(random_device{}()), dist(0.0, 1.0), L_dist(0, L - 1) {}
// defines the size of spin array, expectation value and boltzmann vectors, as well as the generator, [0,1] and [0,L]

void IsingModel::initialize(){ // initializes the spins randomly
    for (int i = 0; i < L; ++i){ // row loop
        for (int j = 0; j < L; ++j){ // column loop
            spins(i, j) = (dist(generator) > 0.5) ? 1 : -1; // if the number generated is larger tha 0.5 the spin is initialized as up (1), and if not, the spin is initialized as down (-1)
        }   
    }
}

void IsingModel::initialize_ordered() { // initializes every spin in the direction of up
    for (int i = 0; i < L; ++i){ // row loop
        for (int j = 0; j < L; ++j){ // column loop
            spins(i, j) = 1; // every spin is intitialized in positive direction 
        }   
    }
}

void IsingModel::metropolis() { // implements the metropolis algorithm
    int dE, x_i, y_i, sums; // defines variables needed for the algorithm
    for (int i = 0; i < L * L; ++i) { // loops through each spin in the LxL lattice
        x_i = L_dist(generator); //randomly generated positions between 0 and L, x direction 
        y_i = L_dist(generator); //randomly generated positions between 0 and L, y direction 
        sums = spins(y_i, (x_i + 1) % L) + spins(y_i, (x_i - 1 + L) % L) + spins((y_i + 1) % L, x_i) + spins((y_i - 1 + L) % L, x_i); // neighboring spin sum
        dE = 2 * sums * spins(y_i, x_i); // delta E calculation 
        if (dE < 0 || dist(generator) <= boltzmann(dE + 8)) { // checks if the energy brings the simulation closer to the critical temperature
            // if it does, then we update the energy, magnetization and flipped spin
            E += dE; // updates the energy 
            M -= 2 * spins(y_i, x_i); // updates the magnetization 
            spins(y_i, x_i) *= -1; // flips the spin
            ++flipcount; // adds a count to the flipcount, since the flip was 
        }
    }
}

void IsingModel::calculate_boltzmann(double T) { // defines the function calculating the boltmann constants
    double beta = 1.0 / T; // beta 1/kT
    // fills the vector containing only zeros to contain the correct boltzmann values
    boltzmann(0) = exp(8.0 * beta); 
    boltzmann(4) = exp(4.0 * beta);
    boltzmann(8) = 1.0;
    boltzmann(12) = exp(-4.0 * beta);
    boltzmann(16) = exp(-8.0 * beta);
    // every other instance of boltzmann constant will be set to 0, by definition
    // the reason we do this, is the energy is only not 0 in these instances, so it only makes sense to make calculations where it has an actual value
} 

void IsingModel::expval(){ // defines the function calculating the expectation values
    expect(0) += E; // updates the energy
    expect(1) += E * E ; // updates the energy^2
    expect(2) += M ; // updates the magnetization
    expect(3) += M * M; // updates the magnetization^2
    expect(4) += abs(M) ; // absoulte value of the magnetization
}

void IsingModel::energycalc(){ // function that calculates energy 
    E = 0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            E -= spins(i, j) * (spins(i, (j + 1) % L) + spins((i + 1) % L, j)); // function (1) in the report
        }
    }
}

void IsingModel::magcalc(){ // function that calculates magnetization 
    M= 0;
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            M += spins(i, j); // function (2) in the report
        }
    }
}

void IsingModel::spinfile(ofstream& file_spins){ // a simple function that writes a spin-array into a file 
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
           file_spins << spins(i, j) << " "; // writes a spin and a space into the file that is already opened
        }
    }
    file_spins << "\n"; // writes a return statement into the file
}

void IsingModel::expfile(int cycle , ofstream& filename){ // writes the expectation values into a file, devided by amounts of cycles performed to obtain the values, as well as the amount of cycles performed
    filename << expect(0) / (cycle + 1) << " " << expect(1)  / (cycle + 1)<< " " << expect(2)  / (cycle + 1)<< " " << expect(3)  / (cycle + 1)<< " " << expect(4) / (cycle + 1) << " " << cycle + 1 <<  "\n";
}

void IsingModel::simulate_log(double T, int cycles) { // is called to run a simulation that makes calculations and logs neccecary data in .txt files
    stringstream stream; // opens files to contain all data we want :
    stream << fixed << setprecision(3) << T;
    ofstream outfile_exp("data/" + filename + "_" + stream.str() + "_exp.txt"),  // expectation values collected in file ending with "_exp.txt"
             outfile_en("data/" + filename + "_" + stream.str() + "_en.txt"),  // energies collected in file ending with "_en.txt"
             outfile_flips("data/" + filename + "_" + stream.str() + "_flips.txt"); // accepted flips of spins collected in file ending with "_flips.txt"
   ofstream outfile_spins;
   if (spin_file == true){ // if we want to obtain a spin_file, the spin_file boolean has to be true
        ofstream outfile_spins("data/" + filename + "_" + stream.str() + "_spins.txt"); // the spin_file is opened
    }
    int N = L * L; // defines amount of points in lattice
    flipcount = 0; // sets the integer  flipcount to 0
    // runs simulations calculationg boltzmann values, energy and magnetization in initial condition
    calculate_boltzmann(T); //boltzmann constant
    energycalc(); //energy calculation for intial condition
    magcalc(); // magnetixation calc for initial condition
    spinfile(outfile_spins); // writes initial spin data into spinfile if wanted
    for (int cycle = 0; cycle <= cycles; ++cycle) {
        metropolis(); // performs metropolis algorithm for wanted amounts of cycles
        expval(); // evaluates expectation values for each instance of cycles performed
        if ((cycle) % (cycles / data_lines) == 0) { // writes data into files, one for expectation values, one for energy and one for amounts of accepted flips of spins
            cout << L << " " << cycle << endl; 
            outfile_en << E << "\n"; // writes energy into a file
            expfile(cycle,  outfile_exp); // writes the energy in this instance os cycle into the file
            outfile_flips << flipcount << "  " << (cycle + 1) * N << "\n"; // writes the amounts of accepted spin-flips into the document
            if (spin_file == true){
                spinfile(outfile_spins); // writes the spin-arrau in this instance os cycle into the file
            }
        }
    }
    // closes the files
    if (spin_file == true){
        outfile_spins.close();
    }
    outfile_exp.close();
    outfile_en.close();
    outfile_flips.close();
}

void IsingModel::simulate(double T, int cycles) { // runs the simulation without documenting data
    int N = L * L;
    calculate_boltzmann(T);
    energycalc();
    magcalc();
    for (int cycle = 0; cycle <= cycles; ++cycle) {
        metropolis();
        expval();
    }
}

void IsingModel::reset() { // resets the expectation values by setting them to 0
    expect.zeros();
}
