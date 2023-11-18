#python program that imports and processes data collected by compiling and running simulation defined in 20x20lattice.cpp, and 20x20tchange.cpp
#collects simulation data for 20x20 lattice: we simulate the behaviour of four different environments of the same size
#for the temperature T=1 and T=2.4 with spins initialized in the same, and random directions. Each plot is saved in directory "../report/figures/", which is consistent with
#compilation of latex file to pdf for the report
#we start by importing needed modules
import numpy as np
import matplotlib.pyplot as plt
#import scipy.interpolate as pol
import pandas as pd
import seaborn as sns
import os

#===========================================================================================================
#--------------------- ENERGY AND MAGNETIZATION AS A FUNCTION OF MONTE CARLO CYCLES  ----------------------#
#===========================================================================================================
L = 20 #size defned
N = L * L
#we import the files created by the simulation code in 20x20lattice.cpp: 
T1_order = np.loadtxt(f"../data/20x20ordered_1.000_exp.txt")
T1_random = np.loadtxt(f"../data/20x20random_1.000_exp.txt")
T24_order = np.loadtxt(f"../data/20x20ordered_2.400_exp.txt")
T24_random = np.loadtxt(f"../data/20x20random_2.400_exp.txt")

#simple plot, expressing how estimated expectation energy and magnetization changes with amounts of monte carlo cycles performed by the simulation
plt.style.use("ggplot")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].set_title(r"$\langle E \rangle$ as a Function of Monte Carlo Cycles")
axes[1].set_title(r" $\langle |M| \rangle $ as a Function of Monte Carlo Cycles")
# any_data[:, 5] is an array containing the sata for cycles performed at the moment the data was recorded by our simulation
axes[0].plot(T1_order[:, 5], T1_order[:, 0] / N, label="T=1, ordered initial spins", color="#88CEEB") 
axes[0].plot(T1_random[:, 5], T1_random[:, 0] / N, label="T=1, random initial spins", color="#CC6677")
axes[0].plot(T24_order[:, 5], T24_order[:, 0] / N, label="T=2.4, ordered initial spins", color="#44AA99")
axes[0].plot(T24_random[:, 5], T24_random[:, 0] / N, label="T=2.4, random initial spins", color="#DDCC77")
axes[1].plot(T1_order[:, 5], T1_order[:, 4] / N, label="T=1, ordered initial spins", color="#88CEEB")
axes[1].plot(T1_random[:, 5], T1_random[:, 4] / N, label="T=1, random initial spins", color="#CC6677")
axes[1].plot(T24_order[:, 5], T24_order[:, 4] / N, label="T=2.4, ordered initial spins", color="#44AA99")
axes[1].plot(T24_random[:, 5], T24_random[:, 4] / N, label="T=2.4, random initial spins", color="#DDCC77")
axes[0].set_xlabel("Amount of Monte Carlo Cycles")
axes[1].set_xlabel("Amount of Monte Carlo Cycles")
axes[0].set_ylabel(r"$\langle E \rangle [J]$")
axes[1].set_ylabel(r"$\langle |M| \rangle $")
axes[0].set_xscale("log") #logaritmic scaling 
axes[1].set_xscale("log")
axes[0].legend()
axes[1].legend()
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../report/figures/20x20randord.pdf")
plt.show()



#===========================================================================================================
#--------------------- 10^6 MONTE CARLO CYCLES, NUMBER OF OCCURANCES OF A GIVEN ENERGY --------------------#
#===========================================================================================================
#Here we create histogram plots for energy in the systems of T=1 and T=2.4 with random initializations
#These plots will vislaize the probability for an energy E to be measured
T1_random_ = np.loadtxt(f"../data/20x20random_1.000_en.txt") #energy for T=1
T24_random_ = np.loadtxt(f"../data/20x20random_2.400_en.txt")  #energy for T=2.4
cycles=10**6
plt.style.use("ggplot")
#creates a plot showing the amounts of instances is measured of certain energies, as fractions of 1
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].set_title(r"Probability for $\langle E \rangle$ at T=1, bins = 10")
axes[1].set_title(r"Probability for $\langle E \rangle$ at T=2.4, bins = 60")
axes[0].hist(T1_random_/(N), label=rf"T=1.0, Average Energy={str(round(np.mean(np.array(T1_random_/(N))),2))}", histtype='bar', color= "#44AA99", bins = 20 , weights = np.ones_like(T1_random_/N)/len(T1_random_/N))
axes[1].hist(T24_random_/(N), label=rf"T=2.4, Average Energy={str(round(np.mean(np.array(T24_random_/(N))),2))}", histtype='bar', color= "#44AA99", bins = 170 , weights = np.ones_like(T24_random_/N)/len(T24_random_/N))
axes[0].set_xlabel(r"$\langle E \rangle [J]$")
axes[1].set_xlabel(r"$\langle E \rangle [J]$")
axes[0].set_ylabel("Probability in fractions of 1")
axes[1].set_ylabel("Probability in fractions of 1")
axes[0].legend()
axes[1].legend()
plt.savefig("../report/figures/20x20probweighted.pdf")
plt.show()
variance1 =(T1_random[-1,1]-(T1_random[-1,0])**2)/N
print("Variance in energy for T = 1:", variance1)
print("Standard deviation in energy for T = 1:", np.sqrt(variance1))
variance2 =(T24_random[-1,1]-(T24_random[-1,0])**2)/N
print("Variance in energy for T = 2.4:", variance2)
print("Standard deviation in energy for T = 2.4:", np.sqrt(variance2))


plt.style.use("ggplot")
fig, axes = plt.subplots(1, 1, figsize=(7, 4))
#creates a plot showing the amounts of instances is measured of certain energies, as percent with logarithmic y-axis
plt.title(r"Probability distribution for $\langle E \rangle$")
axes.hist(T1_random_/(N), label=r"Probability for $\langle E \rangle$ at T=1, bins = 120", color= "#DDCC77", bins = 120 , density = True)
axes.hist(T24_random_/(N), label=r"Probability for $\langle E \rangle$ at T=2.4 , bins = 60", color= "#44AA99", bins = 60, density = True)
axes.set_xlabel(r"$\langle E \rangle [J]$")
axes.set_ylabel("Percentage [%]")
axes.legend()
axes.set_yscale("log")
plt.savefig("../report/figures/20x20prob.pdf")
plt.show()


#===========================================================================================================
#----------------------- ACCEPTED FLIPS AS A FUNCTION OF T AT 10^6 MONTE CARLO CYCLES ---------------------#
#===========================================================================================================
#Here we import the data created by the simulation code in 20x20tchange.cpp, and plot amount of accepted flips after all cycles 10^6 have been performed as a function of the temperature
plt.style.use("ggplot") 
cycles=10**6
fig, axes = plt.subplots(1, 1, figsize=(6, 5))
t_array=np.arange(0 , 6.0 , 0.1)
flip_array=np.zeros(len(t_array))
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77" ,"slategray", 'lightcoral', 'lightgreen', 'cadetblue', "#117733", "#661100", "gold","forestgreen"]
for i in range(len(t_array)):
    flips = np.loadtxt(f"../data/20x20tchange_{t_array[i]:.3f}_flips.txt") #loads file and reads data
    flip_array_element=(flips[:,0][-1])
    flip_array[i]=(flip_array_element/(cycles*L*L)*100) # amount of flips devided by amout of spins in environment, times 100 to get the plot in %
#we find the indexes where the function for accepted flips of spins is the steepest
grad = np.gradient(flip_array)
index1 = np.where(t_array>2.5)[0][0]
index2 = np.where(t_array<2)[0][-1]
maxgrad = np.amax(abs(grad[index2:index1]))
index = np.where(abs(grad) == maxgrad)[0][0]
axes.plot(t_array[index] ,flip_array[index], "X", color="#CC6677", label=f"Critical point [{t_array[index]:.3f}, {flip_array[index]:.3f}]") # plots the point in which the model is steepest
axes.scatter(t_array, flip_array, color="#44AA99")
axes.set_xlabel(r"Temperature T $[J/k_B]$")
axes.set_xlabel(r"Temperature T $[J/k_B]$")
axes.set_ylabel("Accepted flips of spins, in percent [%]")
axes.legend()
plt.savefig("../report/figures/20x20acceptflips_t_EXTRAASF.pdf")
plt.show()