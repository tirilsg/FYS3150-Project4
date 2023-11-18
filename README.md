# Project4 FYS3150 : A 2D Ising Model

This repository contains c++ code constructing a model for a 2D Ising model, and using this model to run different simulations exploring the model and the physics behind it. The reopository is devided into maps `src`, `include`, `data`, `python` and `report`. Outside of these maps, in our main repository, the c++ files that run different kinds of simulations are stored. The model itself is defined in two seperate files, `IsingModel.cpp` and `IsingModel.hpp`, stored in maps `src` (for source files .cpp) and  `include` (for header files .hpp) dependent on their filetype. When running any simulations provided by the code in this repository, these maps have to be linked together when compiling.

----------------------
## The Class Definition of the 2D Ising Model
The Ising Model is declared in the header file `IsingModel.hpp`, and further defined in our source file `IsingModel.cpp`. The purpouse of these definitions is to model the system, and then bring it as close to equilibrium as possible, which we do by developing methods based on theory presented in the report for our project.

### `IsingModel.hpp`:
Contains a class `IsingModel` that defines and models the behaviour of particles organized in a lattice structure of size L at temperature T.
The header file defines a good amount of methods:
 * `initialize()` and `initialize_ordered()` that will be used to define the direction of the intial spins in our system.
 * `calculate_boltzmann(double T)`, `expval()`, `energycalc()` and `magcalc()` that will be used to calculate boltzmann constants, expected values of energy E, magnetization M, $E^2$, $M^2$ and absoulte magnetization, as well as a short definition of the formulas for energy and magnetization calculations as function of spin.
 * `metropolis()` , `simulate_log(double T, int nr_cycles)`, `simulate(double T, int nr_cycles)`, `reset()` which are functions that will make use of the previous functions to flip a single spin at a time, and evaluate whether the system is closer to equilibrium or not, for an arbitrary amount of cycles, and either logs the data by `simulate_log()` or don't by `simulate()`, and then resets the system's expected values by `reset()`.
 * `spinfile(ofstream& file_spins)` and  `expfile(int cycle, ofstream& filename)` are methods that are used to write data into .txt files. 

### `IsingModel.cpp`:
Defines the methods and functions for the IsingModel class presented above: 
* `initialize()` randomly initializes a lattice of size L x L with spins in conditions either +1 $\uparrow$ or -1 $\downarrow$
* `initialize_ordered()` initializes a lattice of size L x L with ordered spins in +1 $\uparrow$ direction
* `metropolis()` implements the metropolis algorithm explained in the report, to flip single spins and update energy and magnetization of the system if this operation closer to equilibrium, and counts the amount of times this algorithm has resulted in an accepted flip
* `calculate_boltzmann(double T)` fills a previously filled with 0's vector with correct values for the boltzmann factor at the temerature T
* `expval()` calculates expectation values defined by functions presented in report
* `energycalc()` and `magcalc()` calculates the total energy and magnetization of the system from equations presented in the report
* `simulate_log(double T, int nr_cycles)` and `simulate(double T, int nr_cycles)` both simulate the behaviour within the system of our lattice by calling the metropolis algorithm a numnber of `nr_cycles` times, or Monte Carlo Cycles, and calculating the expectation values each time `metropolis()` is called.
* `reset()` fills the vector containing expectation values with zeros


----------------------
## An Overview of Simulation-files in the Main Directory of the Repository:

We make use of this class when performing all of the following simualtions. All data collected by these simulations is stored as .txt files in the directory '/data'. 

### `2x2lattice.cpp`:
Simulates and logs the behaviour of a 2 x 2 sized lattice for monte carlo cycles $10, 10^2 10^3, 10^4, 10^5$, in a temperature intervall $T\epsilon [0.5,4]$ with temperaturesteps dT=0.1. Data is stored with the base-name "2x2lattice" in our data-folder, with names _exp.txt, _flips.txt and _en.txt explaining their contents expectation values, amount of accepted flips, and energy. 

### `20x20lattice.cpp`:
Simulates the behaviour of four different systems $20\times 20$ lattices - two at temperature T=1 with different methods of initializing spins `initialize()` and `initialize_ordered()`, and the other two at temperature T=2.4 also with different methods of initializing spins. An amount $10^6$ cycles is used in this simulation. The data of ordered spins is stored in files with base "20x20ordered", and the data for the randomly intialized spins are stored in files with base "20x20random".

### `20x20tchange.cpp`:
Simulates the behaviour of a single $20\times 20$ lattice, by 10^6$ cycles over an interval of temperatures $T\epsilon [0,6]$ with temperaturesteps dT=0.1, and logs the data in .txt files. The data is saved in a .txt file with the base of the name of the file "20x20tchange". 


### `timerpar.cpp`:
Runs a simulation making use of a parallelization of our code, to log how long it takes simulating an environment of $10^5$ monte carlo cycles, for an arbitrary temperature-interval, for different size lattices. The data is saved in a .txt file with the base of the name of the file "tmr5par". 

### `timerreg.cpp`:
Runs a simulation making use of a non-parallized version of our code `timerpar.cpp`,and logs the timing data in a file "tmr5reg". 

### `varyingsize.cpp`:
Makes use of the parallization method to simulate the behaviour of lattices of different sizes , and logging exoected values for an interval of temperature $T\epsilon [2,2.5]$ with temperature-steps dT=0.01 with an usage of by 10^5$ cycles.

### `varyingsizezoom.cpp`:
Includes the same simulation as `varyingsize.cpp`, but in a smaller temperature intervall $T\epsilon [2.2,2.4]$ with temperature-steps dT=0.005.

--------------------

## Linking and Compiling of Our Project Files:
Each simulation is ran by running the following commands in the terminal, while beinf located in the correct directory:

```sh
g++ filename.cpp src/*.cpp -I include -o filename -O2 -llapack -lblas -larmadillo -fopenmp

```
```sh
./filename
```

----------------------

## Visualizing Data By Running Python Programs:
When ran, the simulations produce sets of data, which can be imported and visualized by running the python programs with the similar names, stored in the map 'python' in our repository: python 2x

### `2x2.py`:
Imports the data files created by the simulation `2x2lattice.cpp`, and creates plots comparing the modelled expectation values to the analytic expressions we obtained in our report. The program also creates a scatter plot showing how the accuracy improves with an increase in monte carlo cycles, as well as a plot of the relative error at $10^5$ monte carlo cycles. The program prints these errors in the expectation of energy and magnetization, as well as the amount of monte carlo cycles to the terminal so that they can be discussed.

### `20x20.py`:
Imports the data files created by the simulations `20x20lattice.cpp` and `20x20tchange.cpp`, and uses it to visualize the expected energies and magnetization of each of the four environments of $20\times 20$ lattices as functions of monte carlo cycles. The program also creates plots using histograms for the different energy states, as well as a plot visalizing how the accepted flips of spins develops as temperature changes.  

### `timer.py`:
Imports the data obtained by running the simulations `timerpar.cpp` and `timerreg.cpp`, and visualizes the difference in run-time when using parallized code, versus not doing so in a plot. 


### `sizevariation.py`:
Imports and loads data collected by compiling and running simulation defined in `varyingsize.cpp` and `varyingsizezoom.cpp`, and uses the data to create plots of expected values as a function of temperature for different size lattices. Further, the program creats a plot highlighting the heat capacity, and uses this to further model the critical temperature. The critical temperatures modelled from this data is printed in the terminal. 

