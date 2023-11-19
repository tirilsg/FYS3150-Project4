#python program that imports and processes data collected by compiling and running simulation defined in timerpar.cpp and timerreg.cpp
#collects simulation data for time it takes for simulations by paralellizes and non parallized code for  varying size lattices over a intervall of temperatures
#uses this data to show difference in efficiency by using and not using parallelized code to sun simulations
#each plot is saved in directory "../report/figures/", which is consistent with compilation of latex file to pdf for the report
#we start by importing needed modules
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

plt.style.use("ggplot")
fig, axes = plt.subplots(1, 1, figsize=(8, 6))
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77", 'lightcoral', "slategray", 'lightgreen', 'cadetblue', "#117733", "#661100", "gold", "forestgreen"]
reg5 = np.loadtxt(f"../data/tmr5reg.txt") #loads and reads the data for the non parallized simulation , timerreg.cpp
par5 = np.loadtxt(f"../data/tmr5par.txt")  #loads and reads the data for the parallized simulation , imerpar.cpp
x = np.linspace(10, 100, 100000) # defines an intervall between 10 and 100, which is in agreeance with the intervall used to run simulations 
regt_5 = UnivariateSpline(reg5[:, 0], reg5[:, 1]) # creates a time(size) function that is adapted to the data for the non parallized simulation
part_5 = UnivariateSpline(par5[:, 0], par5[:, 1]) # creates a time(size) function that is adapted to the data for the parallized simulation
axes.scatter(par5[:, 0], par5[:, 1], color=colors[1]) # plots data as a scatter-plot
axes.scatter(reg5[:, 0], reg5[:, 1], color=colors[0]) # plots data as a scatter-plot
axes.plot(x, part_5(x), color=colors[1], label="Parallelized model") # plots the function adapted to the data 
axes.plot(x, regt_5(x), color=colors[0], label="Regular model") # plots the function adapted to the data 
axes.set_title(r"Runtimes for simulations as functions of lattice size, $10^5$ MC Cycles")
axes.set_xlabel("L")
axes.set_ylabel("Time Duration [s]")
axes.legend()
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../report/figures/duration.pdf") # saves plot
plt.show()
# prits lattice size, as well as time simulation takes for regular and parallized coded simulation to run
print(r"$10^5$ Monte Carlo Cycles")
print("L:       Regular time [s]:       Parallized time [s]:")
for i in range(len(reg5[:,0])):
    print(f"{reg5[i,0]}       {reg5[i,1]}              {par5[i,1]}")