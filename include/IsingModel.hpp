// Header file containint the definition of the structure of our IsingModel class, that we edit and define in the file IsingModel.cpp
// Linked by include-statement in program files that wishes to make use of this class to run simulations

#ifndef __IsingModel_hpp__  
#define __IsingModel_hpp__

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>
#include <omp.h>

using namespace std;
using namespace arma;

class IsingModel {
public:
    int L;
    bool spin_file = false; // we set the boulerian spin_file to false as default. That way, if we wish to log instances of spin for each simulation, that can be done by manually setting this to true
    bool flip_file = false; // we set the boulerian flip_file to false as default. That way, if we wish to log instances of amount of flipped spins for each simulation, that can be done by manually setting this to true
    // definition of a couple variables that will be needed to run simulation 
    double E = 0; //energy
    double M = 0; //magnetization
    double T = 0; //temperature
    int flipcount = 0; //variable that will be changed to count amount of accepted flips for our systen
    string filename; // stem of the name of a file we want to log information in 
    int data_lines;
    // vectors that will contain updated information about the state of the system 
    vec boltzmann; // calculated boltzmann factors
    vec expect; // calculated expectation values
    Mat<int> spins; //containing the state of spins in the system, up (1) or down (-1)

    // statements needed for random number-generation
    mt19937 generator; // using mt19937 for random numbers
    uniform_real_distribution<double> dist; // is going to contain random numbers in the intervall [0,1]
    uniform_int_distribution<int> L_dist; // is going to contain random numbers in the intervall [0,L]

    IsingModel(int size, string save_name, int amount_of_data); 
    // we define a couple functions needed to run simulations : 
    void initialize(); // a function that will initialize spins randomly
    void initialize_ordered();  // a function that will initialize spins in the same direction
    void metropolis(); // implementation of the metropolis algorithm
    void calculate_boltzmann(double T); // a function that only calculates boltzmann values for each instance of energies at a temperature T
    void expval(); // calculates expectation values at the currant instance of the system
    void energycalc(); // calculates the energy in the system
    void magcalc(); // calculates the magnetization in the system
    void spinfile(ofstream& file_spins); //writes the state of spins into a file
    void expfile(int cycle, ofstream& filename); // writes expectation values into a file
    void simulate_log(double T, int nr_cycles); // is called to run a simulation that makes calculations and logs neccecary data in .txt files
    void simulate(double T, int nr_cycles); // is called to run a simulation without writing data to files
    void reset(); //resets the entire expectation values
};

#endif