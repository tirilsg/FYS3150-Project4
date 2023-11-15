# Project4 FYS3150 : A 2D Ising Model

This repository contains c++ code constructing a model for a 2D Ising model, and using this model to run different simulations exploring the model and the physics behind it. The reopository is devided into maps `src`, `include`, `data`, `python` and `report`. Outside of these maps, in our main repository, the c++ files that run different kinds of simulations are stored. The model itself is defined in two seperate files, `IsingModel.cpp` and `IsingModel.hpp`, stored in maps `src` (for source files .cpp) and  `include` (for header files .hpp) dependent on their filetype. When running any simulations provided by the code in this repository, these maps have to be linked together when compiling.

----------------------
## The class definition of the 2D Ising Model
The Ising Model is declared in the header file `IsingModel.hpp`, and further defined in our source file `IsingModel.cpp`. The purpouse of these definitions is to model the system, and then bring it as close to equilibrium as possible, which we do by developing methods based on theory presented in the report for our project.

### `IsingModel.hpp`:
Contains a class `IsingModel` that defines and models the behaviour of particles organized in a lattice structure of size L at temperature T.
The header file defines a good amount of methods:
 * `initialize()` and `initialize_ordered()` that will be used to define the direction of the intial spins in our system.
 * `calculate_boltzmann(double T)`, `expval()`, `energycalc()` and `magcalc()` that will be used to calculate boltzmann constants, expected values of energy E, magnetization M, $E^2$, $M^2$ and absoulte magnetization, as well as a short definition of the formulas for energy and magnetization calculations as function of spin.
 * `metropolis()` , `equilibrate(double T, int nr_cycles)`, `simulate_log(double T, int nr_cycles)`, `simulate(double T, int nr_cycles)`, `reset()` which are functions that will make use of the previous functions to flip a single spin at a time, and evaluate whether the system is closer to equilibrium or not, for an arbitrary amount of cycles, and either logs the data by `simulate_log()` or don't by `simulate()`, and then resets the system's expected values by `reset()`.
 * `spinfile(ofstream& file_spins)` and  `expfile(int cycle, ofstream& filename)` are methods that are used to write data into .txt files. 

### `IsingModel.cpp`:
Defines the methods and functions for the IsingModel class presented above




----------------------
## An overview of simulation-files in the main directory of the repository:

We make use of this class when performing all of the following simualtions. 

### `2x2lattice.cpp`:
### `20x20lattice.cpp`:
### `20x20tchange.cpp`:
### `timerpar.cpp`:
### `timerreg.cpp`:
### `varyingsize.cpp`:






------------------------
notes:
### `particle_class.h`:
Contains a class `Particle` that is used to define a particle by mass, charge, a position vector and a velocity vector. This class also contains a function `escape_test()` that checks whether the particle exists within the bounds of the penning trap. 

### `penningtrap_class.h`:
Contains a definition of the class `PenningTrap`, that defines the environment within a penning trap, and uses this to define methods `evolve_forward_Euler(dt,t)` and `evolve_RK4(dt,t)` that estimates particle's movement within the trap, with the method `evolve_RK4(dt,t)` set as the default trajectory-estimation-method. The class contains definitions of Booleans that defines whether the system is modelled as time-dependent (which is set to `false` by default), and whether the system estimates particle-trajectories by taking into consideration forces between particles that exist within the trap (also `false` by default). Also contains a method `particle_add()` that adds a particle into our system, and a function `count_particles()` that counts the amount of particles within the system.

### `filltrap.cpp`:
Contains a single function `fillPenningTrapWithRandomParticles(trap, number_of_particles)` that adds a number of Calcium ions into a system PenningTrap trap. 

###`logging.cpp`:
Contains a single function `saveDataToTxt(filename, data)` that takes an arbitrary filename, as well as a set of data, and writes the data into the file. 

### `twoqmotion.cpp`:
Contains code for a function `twoq(article1, particle2, max_t, iterations,filename, particle_int)` that adds two particles into a penning trap, and simulates their movement within the trap in an interval of time, and saves the data into a file by calling  `saveDataToTxt(filename, data)`. 

### `zdirmotion.cpp`:
Contains the implementations of a function `singleqmotion(particle, max_t, iterations, filename, useRK4)` that simulates the movement of a single particle within the environment of a penning trap. This function takes a Boolean useRK4 as argument, that dictates whether we use the method `evolve_RK4(dt,t)` to estimate the particle trajectory or not. The data is saved in a file.


------------------------

The simulations themselves is done by calling the functions implemented in the files above, which we simply link to each of our main programs. We decided to split the program code into two separate files, dependent on time-dependency, since the program `time_dependent.cpp` takes a long time to run.

### `time_independent.cpp`:
All the simulations done in this program is independent of time.
This file contains code that estimates the movement of a single particle 1 in a penning trap, two particles 1 and 2 in a penning trap with and without particle interactions, and simulations for usage of both `evolve_forward_Euler(dt,t)` and `evolve_RK4(dt,t)` for a single particle 1. 


### `time_dependent.cpp`:
All the simulations done in this program is dependent of time.
This program defines a function `simulateAndLogData(trap, filename, w_v, fs, max_t, iterations)` that fills a penning trap with particles, and estimates the trajectories of each of these particles by `evolve_RK4(dt,t)` and estimates how many particles are still trapped, for a vector containing different values for amplitudes and frequencies $w_v$ and stores the data for frequencies, and amount of trapped particles in a file. This function is called in our `main()`, and is used to create multiple simulations for a time-dependent system.

--------------------

### Linking and compiling of our project files:
To run the program `time_independent.cpp`:
```sh
g++ time_independent.cpp src/*.cpp -I include -o time_independent -O2 -llapack -lblas -larmadillo
```
```sh
./time_independent
```
To run the program `time_dependent.cpp`:
```sh
g++ time_dependent.cpp src/*.cpp -I include -o time_dependent -O2 -llapack -lblas -larmadillo
```
```sh
./time_dependent
```

----------------

The simulations ran by the programs `time_independent.cpp` and `time_dependent.cpp` produce sets of data, which can be imported and visualized by running the python programs with the similar names: 

### `time_independent.py`:
By running our program `time_independent.cpp`, we will obtain a number of files, which we can interpret and visualize by running the code in this python file. The program imports these files, and creates plots showing the particles movement in the xyz-plane, as well as relevant phase space plots, and an error analysis. 


### `time_dependent.py`:
By running our program `time_dependent.cpp`, we will obtain a number of files, which we can interpret and visualize by running the code in this python file. The program imports these files, and returns plots visualizing the fraction of trapped particles as a function of the frequency $w_v$.

