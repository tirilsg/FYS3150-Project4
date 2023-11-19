#python program that imports and processes data collected by compiling and running simulation defined in 2x2lattice.cpp
#collects simulation data for 2x2 lattice over a temperature intervall 0.5-4.0 with dT 0.1, for varying amounts of monte carlo cycles 
#we start by importing needed modules
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from math import*
from numpy import polyfit


#===========================================================================================================
#------------------ COMPARISON BETWEEN MODELLED AND ANAYTICALLY OBTAINED EXPECTATION VALUES----------------#
#----------------------- FOR A SUFICCIENTLY LARGE AMOUNT OF MONTE CARLO CYCLES (10^5) ---------------------#
#===========================================================================================================
N=4 ; k=1 ; J=1 #defines variables
t_array = np.arange(0.5, 4.00, 0.1) #array containing temperatures
#define arrays containing 0 for both expectation values modelled and analytically obtained
E = np.zeros(len(t_array)) ; EE = np.zeros(len(t_array)) 
M = np.zeros(len(t_array)) ; MM = np.zeros(len(t_array)) 
exp_E = np.zeros(len(t_array)) ; exp_EE = np.zeros(len(t_array))
exp_M = np.zeros(len(t_array)) ; exp_MM = np.zeros(len(t_array))
exp_cp= np.zeros(len(t_array)) ; exp_chi = np.zeros(len(t_array))
for i in range(len(t_array)): #looping over every temperature in our array
    beta=1/t_array[i] # calculates beta for this instance
    # estimations for the analytical values for a 2x2 lattice, explained and derived in the report 
    Z=4*(3+cosh(8*beta*J)) # function 17 in the report
    exp_E[i] = -32*J*np.sinh(8*beta*J)/Z  # function 18 in the report
    exp_EE[i] = 4*64*J*J*np.cosh(8*beta*J)/Z  # function 19 in the report
    exp_M[i] = 8*(2 + np.exp(8*beta*J))/Z  # function 20 in the report
    exp_MM[i] = 32*(1 + np.exp(beta*8*J))/Z   # function 21 in the report
    exp_cp[i] = 1/(k*(t_array[i])*(t_array[i]))*(exp_EE[i]-exp_E[i]**2)  # function 22 in the report
    exp_chi[i] = 1/(k*(t_array[i]))*(exp_MM[i]-exp_M[i]**2)  # function 23 in the report
    data = np.loadtxt(f"../data/2x2lattice5_{t_array[i]:.3f}_exp.txt") #opens the data file containing simulation data for 10^5 monte carlo cycles
    E_= data[:,0]; EE_= data[:,1]; M_= data[:,2]; MM_= data[:,3]; abM_= data[:,4] ; cycles= data[:,5] #extracts the data we need
    #fills the arrays with the modelled estimation values, that corresponts to the last value in each data-set, since this is when the system has reached equilibrium and we're sufficiently close to the critical point
    E[i] = E_[-1] 
    M[i] = abM_[-1]
    EE[i] = EE_[-1]
    MM[i] = MM_[-1]
#defining plot, where each of the expectation values are visualized and compared (modelled and analytic) in different subplots
plt.style.use("ggplot") 
fig, axes = plt.subplots(2, 3, figsize=(10, 6))
#fig.suptitle(r"Comparison of Modelled and Analytic Expectationvalues, for $10^5$ Monte Carlo Cycles")
axes[0, 0].scatter(t_array, E/N, label="Modelled E", zorder=1, color='#1A7DA8')
axes[0, 0].plot(t_array, exp_E/N, label="Expected E" )
axes[0, 0].set_title("Expected E")
axes[1, 0].scatter(t_array, EE/(N**2), label=r"Modelled $E^2$", zorder=1, color='#1A7DA8')
axes[1, 0].plot(t_array, exp_EE/(N**2), label=r"Expected $E^2$" )
axes[1, 0].set_title(r"Expected $E^2$")
axes[0, 1].scatter(t_array, M/N, label="Modelled M", zorder=1, color='#1A7DA8')
axes[0, 1].plot(t_array, exp_M/N, label="Expected M" )
axes[0, 1].set_title("Expected M")
axes[1, 1].scatter(t_array, MM/(N**2), label=r"Modelled $M^2$", zorder=1, color='#1A7DA8')
axes[1, 1].plot(t_array, exp_MM/(N**2), label=r"Expected $M^2$" )
axes[1, 1].set_title(r"Expected $M^2$")
axes[0, 2].scatter(t_array, (EE-E**2)/(N*t_array*t_array), zorder=1, label=r"Modelled $C_V$", color='#1A7DA8')
axes[0, 2].plot(t_array, exp_cp/N, label=r"Expected $C_V$" )
axes[0, 2].set_title(r"Expected $C_V$ [$k_B$]", fontsize=14)
axes[1, 2].scatter(t_array, (MM-M**2)/(N*t_array), zorder=1, label=r"Modelled $\chi$", color='#1A7DA8')
axes[1, 2].plot(t_array, exp_chi/N, label=r"Expected $\chi$" )
axes[1, 2].set_title(r"Expected $\chi$ $[1/J]$", fontsize=14)
axes[0, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 0].set_ylabel(r"$\langle E \rangle [J]$")
axes[1, 0].set_ylabel(rf"$\langle E^2 \rangle [J]$")
axes[0, 1].set_ylabel(r"$\langle |M| \rangle $")
axes[1, 1].set_ylabel(r"$\langle M^2 \rangle $")
axes[0, 2].set_ylabel(r"$\langle C_V \rangle $ [$k_B$]")
axes[1, 2].set_ylabel(r"$\langle \chi\rangle $ [$1/J$]")
for ax in axes.flat:
    ax.legend()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("../report/figures/expected_values.pdf") #the figure is saved in a directory ../report/figures/, which is the directory deferenced when compiling the latex file our report is written in 
plt.show()

#===========================================================================================================
#------------------ PLOTS OF ERROR, AND ERROR INCORPERATED IN ANALYTIC/NUMERIC COMPARISON -----------------#
#===========================================================================================================
#calculates difference between the modelled and analytic values for expectation values, and saves these as seperate arrays
error_E = np.abs(E - exp_E)
error_EE = np.abs(EE - exp_EE)
error_M = np.abs(M - exp_M)
error_MM = np.abs(MM - exp_MM)
error_cp = np.abs((EE-E**2)/(t_array*t_array) - exp_cp)
error_chi = np.abs((MM-M**2)/(t_array) - exp_chi)

#This part of the code is pretty unneccecary, but plots the errors as function of temperature if one wants to study this
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle(r"Absolute Errors between Modelled and Analytic Expectationvalues, $10^5$ Monte Carlo Cycles")
axes[0, 0].plot(t_array, error_E)
axes[0, 0].scatter(t_array, error_E, label=f"Average error {np.mean(error_E):.3f}", color='#1A7DA8')
axes[0, 0].set_title("Absolute Error for E")
axes[0, 0].legend()
axes[0, 1].plot(t_array, error_EE)
axes[0, 1].scatter(t_array, error_EE, label=f"Average error {np.mean(error_EE):.3f}", color='#1A7DA8')
axes[0, 1].set_title(r"Absolute Error for $E^2$")
axes[0, 1].legend()
axes[1, 0].plot(t_array, error_M)
axes[1, 0].scatter(t_array, error_M,label=f"Average error {np.mean(error_M):.3f}",  color='#1A7DA8')
axes[1, 0].set_title("Absolute Error for M")
axes[1, 0].legend()
axes[1, 1].plot(t_array, error_MM)
axes[1, 1].scatter(t_array, error_MM, label=f"Average error {np.mean(error_MM):.3f}", color='#1A7DA8')
axes[1, 1].set_title(r"Absolute Error for $M^2$")
axes[1, 1].legend()
axes[0, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 0].set_ylabel(r"$\langle E \rangle [J]$")
axes[1, 0].set_ylabel(rf"$\langle E^2 \rangle [J]$")
axes[0, 1].set_ylabel(r"$\langle |M| \rangle $")
axes[1, 1].set_ylabel(r"$\langle M^2 \rangle $")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("../report/figures/only_errors.pdf")
plt.show()

#creates the plot as for "../report/figures/expected_values.pdf", but with included errorbars for each modelled point 
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 3, figsize=(10, 6))
fig.suptitle(r"Comparison of Modelled and Analytic Values, $10^5$ MC Cycles, with errorbars", fontsize=16)
axes[0, 0].errorbar(t_array, E / N, yerr=error_E / N, fmt='o', label="Modelled E", color='#1A7DA8')
axes[0, 0].plot(t_array, exp_E / N, label="Expected E")
axes[0, 0].set_title("Expected E")
axes[0, 0].legend()
axes[1, 0].errorbar(t_array, EE / (N ** 2), yerr=error_EE / (N ** 2), fmt='o', label=r"Modelled $E^2$", color='#1A7DA8')
axes[1, 0].plot(t_array, exp_EE / (N ** 2), label=r"Expected $E^2$")
axes[1, 0].set_title(r"Expected $E^2$")
axes[1, 0].legend()
axes[0, 1].errorbar(t_array, M / N, yerr=error_M / N, fmt='o', label="Numerical M", color='#1A7DA8')
axes[0, 1].plot(t_array, exp_M / N, label="Modelled M")
axes[0, 1].set_title("Expected M")
axes[0, 1].legend()
axes[1, 1].errorbar(t_array, MM / (N ** 2), yerr=error_MM / (N ** 2), fmt='o',label=r"Modelled $M^2$", color='#1A7DA8')
axes[1, 1].plot(t_array, exp_MM / (N ** 2), label=r"Expected $M^2$")
axes[1, 1].set_title(r"Expected $M^2$")
axes[1, 1].legend()
axes[0, 2].errorbar(t_array, (EE - E ** 2) / (N * t_array * t_array), yerr=error_cp / N, fmt='o', label=r"Modelled $C_V$", color='#1A7DA8' )
axes[0, 2].plot(t_array, exp_cp / N, label=r"Expected $C_V$")
axes[0, 2].set_title(r"$Expected C_V$ [$k_B$]", fontsize=14)
axes[0, 2].legend()
axes[1, 2].errorbar(t_array, (MM - M ** 2) / (N * t_array), yerr=error_chi / N, fmt='o', label=r"Modelled $\chi$", color='#1A7DA8' )
axes[1, 2].plot(t_array, exp_chi / N, label=r"Expacted $\chi$")
axes[1, 2].set_title(r"$Expected \chi$ $[J^{-1}]$", fontsize=14)
axes[1, 2].legend()
axes[0, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 0].set_ylabel(r"$\langle E \rangle [J]$")
axes[1, 0].set_ylabel(rf"$\langle E^2 \rangle  [J]$")
axes[0, 1].set_ylabel(r"$\langle |M| \rangle $")
axes[1, 1].set_ylabel(r"$\langle M^2 \rangle $")
axes[0, 2].set_ylabel(r"$\langle C_V \rangle $ [$k_B$]")
axes[1, 2].set_ylabel(r"$\langle \chi \rangle $ [$1/J$]")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("../report/figures/expected_values_errors.pdf")
plt.show()


#===========================================================================================================
#------------------------ COMPARISON OF AMOUNTS OF MONTE CARLO CYCLES, ALL VARIABLES-----------------------#
#===========================================================================================================

# This code-snippet imports every file created by the simulation 2x2lattice.cpp, and plots, as previously done, all our modelled expectation values in the same plots
#by doing so we obtain a visualization of how the accuracy of the simulation increases with increasing amounts of monte carlo cycles performed
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77" , 'lightcoral', "slategray", 'lightgreen', 'cadetblue', "#117733", "#661100", "gold","forestgreen"]
cycle_ = [10, 100, 1000, 10000, 100000] 
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
number = 0
while number <=4: #looping throuh values 0-4 that will be used to open and read data from each amount of monte carlo cycles performed
    for b in cycle_: # same loop as previously performed
        t_array = np.arange(0.5, 4.00, 0.1)
        E = np.zeros(len(t_array)) ; EE = np.zeros(len(t_array))
        M = np.zeros(len(t_array)) ; MM = np.zeros(len(t_array))
        exp_E = np.zeros(len(t_array)) ; exp_EE = np.zeros(len(t_array))
        exp_M = np.zeros(len(t_array)) ; exp_MM = np.zeros(len(t_array))
        exp_cp = np.zeros(len(t_array)) ; exp_chi = np.zeros(len(t_array))
        for i in range(len(t_array)):
            beta=1/t_array[i]
            Z=4*(3+cosh(8*beta*J))
            exp_E[i] = -32*J*np.sinh(8*beta*J)/Z
            exp_EE[i] = 4*64*J*J*np.cosh(8*beta*J)/Z
            exp_M[i] = 8*(2 + np.exp(8*beta*J))/Z
            exp_MM[i] = 32*(1 + np.exp(beta*8*J))/Z 
            exp_cp[i] = 1/(k*(t_array[i])*(t_array[i]))*(exp_EE[i]-exp_E[i]**2)
            exp_chi[i] = 1/(k*(t_array[i]))*(exp_MM[i]-exp_M[i]**2)
            data = np.loadtxt(f"../data/2x2lattice{number+1}_{t_array[i]:.3f}_exp.txt") #opens the file with 10^(number+1) Mc cycles performed at temperature t_array[i]
            E_= data[:,0]; EE_= data[:,1]; M_= data[:,2]; MM_= data[:,3]; abM_= data[:,4] ; cycles= data[:,5]
            E[i] = E_[-1]
            M[i] = abM_[-1]
            EE[i] = EE_[-1]
            MM[i] = MM_[-1]
        #plots data as scatter-plot for each instance of MC cycles
        axes[0, 0].scatter(t_array, E/N, label=f"{b} Monte Carlo cycles", zorder=1, color=colors[number])
        axes[1, 0].scatter(t_array, EE/(N**2), label=f"{b} Monte Carlo cycles", zorder=1, color=colors[number])
        axes[0, 1].scatter(t_array, M/N, label=f"{b} Monte Carlo cycles", zorder=1, color=colors[number])
        axes[1, 1].scatter(t_array, MM/(N**2), label=f"{b} Monte Carlo cycles", zorder=1, color=colors[number])
        axes[0, 2].scatter(t_array, (EE-E**2)/(N*t_array*t_array), zorder=1, label=f"{b} Monte Carlo cycles", color=colors[number])
        axes[1, 2].scatter(t_array, (MM-M**2)/(N*t_array), zorder=1, label=f"{b} Monte Carlo cycles", color=colors[number])
        number+=1
#also plots the analytically obtained expectation values to compare
axes[0, 0].set_title("Expected E")
axes[1, 0].set_title(r"Expected $E^2$")
axes[0, 1].set_title("Expected M")
axes[1, 1].set_title(r"Expected $M^2$")
axes[1, 2].set_title(r"$Expected \chi$ $[J^{-1}]$", fontsize=14)
axes[0, 2].set_title(r"$Expected C_V$ [$k_B$]", fontsize=14)
axes[0, 0].plot(t_array, exp_E/N, label="Expected E")
axes[1, 0].plot(t_array, exp_EE/(N**2), label=r"Expected $E^2$")
axes[0, 1].plot(t_array, exp_M/N, label="Expected M")
axes[1, 1].plot(t_array, exp_MM/(N**2), label=r"Expected $M^2$")
axes[0, 2].plot(t_array, exp_cp/N, label=r"Expected $C_v$")
axes[1, 2].plot(t_array, exp_chi/N, label=r"Expected $\chi$")
axes[0, 0].legend(fontsize='small')
axes[1, 0].legend(fontsize='small')
axes[0, 1].legend(fontsize='small')
axes[1, 1].legend(fontsize='small')
axes[0, 2].legend(fontsize='small')
axes[1, 2].legend(fontsize='small')
axes[0, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1, 2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0, 0].set_ylabel(r"$\langle E \rangle [J]$")
axes[1, 0].set_ylabel(rf"$\langle E^2 \rangle [J]$")
axes[0, 1].set_ylabel(r"$\langle |M| \rangle $")
axes[1, 1].set_ylabel(r"$\langle M^2 \rangle$")
axes[0, 2].set_ylabel(r"$\langle C_V \rangle$ [$k_B$]")
axes[1, 2].set_ylabel(r"$\langle \chi \rangle$ [$1/J$]")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("../report/figures/mcc_comp_all.pdf")
plt.show()

#===========================================================================================================
#------------------------------- COMPARISON OF AMOUNTS OF MONTE CARLO CYCLES ------------------------------#
#===========================================================================================================
#creates the same plot as previously, but zooms into only three interesting developments for energy, magnetization and heat capacity as function of temperature, and MC cycles
print("MC Cycles:       Average relative error E:     Average relative error abs(M):")
plt.style.use("ggplot")
fig, axes = plt.subplots(1, 3, figsize=(12, 5))
number = 0
while number <=4:
    for b in cycle_:
        t_array = np.arange(0.5, 4.00, 0.1)
        E = np.zeros(len(t_array)) ; EE = np.zeros(len(t_array))
        M = np.zeros(len(t_array)) ; MM = np.zeros(len(t_array))
        exp_E = np.zeros(len(t_array)) ; exp_EE = np.zeros(len(t_array))
        exp_M = np.zeros(len(t_array)) ; exp_MM = np.zeros(len(t_array))
        exp_cp = np.zeros(len(t_array)) ; exp_chi = np.zeros(len(t_array))
        for i in range(len(t_array)):
            beta=1/t_array[i]
            Z=4*(3+cosh(8*beta*J))
            exp_E[i] = -32*J*np.sinh(8*beta*J)/Z
            exp_EE[i] = 4*64*J*J*np.cosh(8*beta*J)/Z
            exp_M[i] = 8*(2 + np.exp(8*beta*J))/Z
            exp_MM[i] = 32*(1 + np.exp(beta*8*J))/Z 
            exp_cp[i] = 1/(k*(t_array[i])*(t_array[i]))*(exp_EE[i]-exp_E[i]**2)
            exp_chi[i] = 1/(k*(t_array[i]))*(exp_MM[i]-exp_M[i]**2)
            data = np.loadtxt(f"../data/2x2lattice{number+1}_{t_array[i]:.3f}_exp.txt")
            E_= data[:,0]; EE_= data[:,1]; M_= data[:,2]; MM_= data[:,3]; abM_= data[:,4] ; cycles= data[:,5]
            E[i] = E_[-1]
            M[i] = abM_[-1]
            EE[i] = EE_[-1]
            MM[i] = MM_[-1]
        axes[0].scatter(t_array, E/N, label=f"{b} Monte Carlo cycles", zorder=1, color=colors[number])
        axes[1].scatter(t_array, M/N, label=f"{b} Monte Carlo cycles", zorder=1, color=colors[number])
        axes[2].scatter(t_array, (EE-E**2)/(N*t_array*t_array), zorder=1, label=f"{b} Monte Carlo cycles", color=colors[number])
        #prints the errors between energy and magnetization calculations into the terminal, as well as the cycle number, to use and present in report 
        print(f"{cycle_[number]}               {np.mean(np.abs(E - exp_E))}             {np.mean(np.abs(M - exp_M))}") 
        number+=1
axes[0].set_title("Expected E")
axes[1].set_title("Expected M")
axes[2].set_title(r"$Expected C_V$ [$k_B$]", fontsize=14)
axes[0].plot(t_array, exp_E/N, label="Expected E" )
axes[1].plot(t_array, exp_M/N, label="Expected M" )
axes[2].plot(t_array, exp_cp/N, label=r"Expected $C_V$" )
axes[0].legend()
axes[1].legend()
axes[2].legend()
axes[0].set_xlabel(r"Temperature T $[J/k_B]$")
axes[1].set_xlabel(r"Temperature T $[J/k_B]$")
axes[2].set_xlabel(r"Temperature T $[J/k_B]$")
axes[0].set_ylabel(r"$\langle E \rangle [J]$")
axes[1].set_ylabel(r"$\langle |M| \rangle $")
axes[2].set_ylabel(r"$\langle C_V \rangle$ [$k_B$]")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("../report/figures/mcc_comp.pdf")
plt.show()
