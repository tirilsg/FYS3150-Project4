#python program that imports and processes data collected by compiling and running simulation defined in varyingsize.cpp, and varyingsizezoom.cpp
#Collects simulation data for varying size lattices in the interval 40-120 with dL 20, and varying temperatures from 2.0 to 2.5 with dT=0.01
#uses data to create plot of expected values as a function of temperature, and uses the heat capacity model to further model the critical temperature
#Each plot is saved in directory "../report/figures/", which is consistent with compilation of latex file to pdf for the report
#if the simulation is to be ran for more lattice-sizes, just append these to the size-list
#we start by importing needed modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import scipy.interpolate as pol
from scipy import stats


#===========================================================================================================
#------------------------- 10^5 MOONTE CARLO CYCLES, ZOOMED INTO INTERVAL 2.0 - 2.5 -----------------------#
#===========================================================================================================
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 2, figsize=(10, 6))
size = [40 , 60 , 80 , 100, 120] #different size lattices
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77" , 'lightcoral', "slategray", 'lightgreen', 'cadetblue', "#117733", "#661100", "gold","forestgreen"]
t_array=np.arange(2.0 , 2.5 , 0.01) #temperature array
n=0
cv_zoom=[] ; tc_values = []
for L in size: #loops through each lattice-size we modelled
    N=L*L #defines N as amount of particles in the lattice
    #defines empty arrays that will contain modelled expectation values E, E^2, M, M^2, CV
    E = np.zeros(len(t_array)) ; EE = np.zeros(len(t_array)) 
    M = np.zeros(len(t_array)) ; MM = np.zeros(len(t_array)) ; CV= np.zeros(len(t_array))
    zoom_cv=[] #defines a list for cv for temperatures within a special intervall
    for i in range(len(t_array)): #loops through each temperature the simulation is ran for 
        data = np.loadtxt(f"../data/mcc5{L}_{t_array[i]:.3f}_exp.txt") #opens and reads the file where expectation values were stored
        E_= data[:,0]; EE_= data[:,1]; M_= data[:,2]; MM_= data[:,3]; abM_= data[:,4] ; cycles= data[:,5] #extracts columns of data
        #defines expectation values as the last calculated value by our simulation, since this is as close to critical temperature that we got
        E[i] = E_[-1]
        M[i] = abM_[-1]
        EE[i] = EE_[-1]
        MM[i] = MM_[-1]
        CV[i] = (EE[i]-E[i]**2)/(N*t_array[i]*t_array[i]) #defines heat capacity as by formula (7) in the report
        if t_array[i]>=2.20 and t_array[i]<=2.35: #we want to create a seperate zoomed in plot, which is done by testing if the temperature is within a special intervall
            zoom_cv.append((EE[i]-E[i]**2)/(N*t_array[i]*t_array[i]))      #if it is, the heat capacity value is added to the list  
        i+=1
    cv_zoom.append(zoom_cv) # the entire list of cv values in our special intervall is added to another list
    #we create an function that is adapted as well as possible to the modelled data for cv by using pol.UnivariateSpline
    cs = pol.UnivariateSpline(t_array,CV) 
    TS=np.linspace(2.0, 2.5, 10000) #array containing 10000 points in our temperature intervall, to make function more accurate
    tc_values.append(TS[np.argmax(cs(TS))]) #appends the temperature that has the same index as the modelled function for cv's apex to a list. this is where the model's critical point is
    #for each size lattice, the data for expectation values E, M, heat capacity and susceptibility is plotted as scatter-plots in their own subplots, with labels explaining what size lattice the data is extracted from
    axes[0 , 0].scatter(t_array[1:], E[1:]/N, label=f"L={L}", color=colors[n])
    axes[1 , 0].scatter(t_array[1:], M[1:]/N, label=f"L={L}", color=colors[n])
    axes[0 , 1].scatter(t_array[1:], (EE[1:]-E[1:]**2)/(N*t_array[1:]*t_array[1:]), label=f"L={L}", color=colors[n])
    axes[1 , 1].scatter(t_array[1:], (MM[1:]-M[1:]**2)/(N*t_array[1:]), label=f"L={L}", color=colors[n])
    axes[0 , 0].set_title(r"Expected Energy")
    axes[1 , 0].set_title(r"Expected Magnetization")
    axes[0 , 1].set_title(r"Expected Heat Capacity")
    axes[1 , 1].set_title(r"Expected Suceptibility")
    axes[0 , 0].set_xlabel(r"T $[J/k_B]$")
    axes[1 , 0].set_xlabel(r"T $[J/k_B]$")
    axes[0 , 1].set_xlabel(r"T $[J/k_B]$")
    axes[1 , 1].set_xlabel(r"T $[J/k_B]$")
    axes[0 , 0].set_ylabel(r"E [J]")
    axes[1 , 0].set_ylabel(r"$\langle |M| \rangle $")
    axes[0 , 1].set_ylabel(r"$C_V [k_B]$")
    axes[1 , 1].set_ylabel(r"$\chi [1/J]$")
    axes[0 , 0].legend()
    axes[1 , 0].legend()
    axes[0 , 1].legend()
    axes[1 , 1].legend()
    n+=1
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../report/figures/sizevariations.pdf") #saves the figure as pdf
plt.show()



plt.style.use("ggplot")
fig, axes = plt.subplots(1, 4, figsize=(12, 3))
size = [40 , 60 , 80 , 100, 120] #different size lattices
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77" , 'lightcoral', "slategray", 'lightgreen', 'cadetblue', "#117733", "#661100", "gold","forestgreen"]
t_array=np.arange(2.205 , 2.4 , 0.005) #temperature array
n=0
cv_zoom=[] ; tc_values = []
for L in size: #loops through each lattice-size we modelled
    N=L*L #defines N as amount of particles in the lattice
    #defines empty arrays that will contain modelled expectation values E, E^2, M, M^2, CV
    E = np.zeros(len(t_array)) ; EE = np.zeros(len(t_array)) 
    M = np.zeros(len(t_array)) ; MM = np.zeros(len(t_array)) ; CV= np.zeros(len(t_array))
    zoom_cv=[] #defines a list for cv for temperatures within a special intervall
    for i in range(len(t_array)): #loops through each temperature the simulation is ran for 
        data = np.loadtxt(f"../data/mcc5zoom{L}_{t_array[i]:.3f}_exp.txt") #opens and reads the file where expectation values were stored, for the zoomed in interval 
        E_= data[:,0]; EE_= data[:,1]; M_= data[:,2]; MM_= data[:,3]; abM_= data[:,4] ; cycles= data[:,5] #extracts columns of data
        #defines expectation values as the last calculated value by our simulation, since this is as close to critical temperature that we got
        E[i] = E_[-1]
        M[i] = abM_[-1]
        EE[i] = EE_[-1]
        MM[i] = MM_[-1]
        CV[i] = (EE[i]-E[i]**2)/(N*t_array[i]*t_array[i]) #defines heat capacity as by formula (7) in the report
        if t_array[i]>=2.20 and t_array[i]<=2.35: #we want to create a seperate zoomed in plot, which is done by testing if the temperature is within a special intervall
            zoom_cv.append((EE[i]-E[i]**2)/(N*t_array[i]*t_array[i]))      #if it is, the heat capacity value is added to the list  
        i+=1
    cv_zoom.append(zoom_cv) # the entire list of cv values in our special intervall is added to another list
    #we create an function that is adapted as well as possible to the modelled data for cv by using pol.UnivariateSpline
    cs = pol.UnivariateSpline(t_array,CV) 
    TS=np.linspace(2.0, 2.5, 10000) #array containing 10000 points in our temperature intervall, to make function more accurate
    tc_values.append(TS[np.argmax(cs(TS))]) #appends the temperature that has the same index as the modelled function for cv's apex to a list. this is where the model's critical point is
    #for each size lattice, the data for expectation values E, M, heat capacity and susceptibility is plotted as scatter-plots in their own subplots, with labels explaining what size lattice the data is extracted from
    axes[0].scatter(t_array[1:], E[1:]/N, label=f"L={L}", color=colors[n])
    axes[1].scatter(t_array[1:], M[1:]/N, label=f"L={L}", color=colors[n])
    axes[2].scatter(t_array[1:], (EE[1:]-E[1:]**2)/(N*t_array[1:]*t_array[1:]), label=f"L={L}", color=colors[n])
    axes[3].scatter(t_array[1:], (MM[1:]-M[1:]**2)/(N*t_array[1:]), label=f"L={L}", color=colors[n])
    axes[0].set_title(r"Expected Energy")
    axes[1].set_title(r"Expected Magnetization")
    axes[2].set_title(r"Expected Heat Capacity")
    axes[3].set_title(r"Expected Suceptibility")
    axes[0].set_xlabel(r"T $[J/k_B]$")
    axes[1].set_xlabel(r"T $[J/k_B]$")
    axes[2].set_xlabel(r"T $[J/k_B]$")
    axes[3].set_xlabel(r"T $[J/k_B]$")
    axes[0].set_ylabel(r"E [J]")
    axes[1].set_ylabel(r"$\langle |M| \rangle $")
    axes[2].set_ylabel(r"$C_V [k_B]$")
    axes[3].set_ylabel(r"$\chi [1/J]$")
    axes[0].legend()
    n+=1
plt.tight_layout()
plt.savefig("../report/figures/sizevariationszoom.pdf") #saves the figure as pdf
plt.show()

#creates a zoomed in plot of the heat capacity data extracted previously
cv_zoom = np.array(cv_zoom) #turns the list containing lists into an array 
plt.style.use("ggplot")
fig, axes = plt.subplots(1, 1, figsize=(8, 6)) #initializes plot
tc_values=[]
for i in range(len(cv_zoom)): #loops through each list within cv_zoom array
    n=np.array(cv_zoom[i]) #turns the list in place i within the array into an array
    t_zoom=np.linspace(2.20, 2.35, len(n)-1) #defiens a temperature-intervall that is zoomed into
    cs = pol.UnivariateSpline(t_zoom,n[1:]) #creates a function that is adapted to the cv-points
    TS=np.linspace(2.20, 2.35, 10000) #defines the temperatureintervall for the model-function 
    list=cs(TS)
    cvmax=0; index=0
    for o in range(len(TS)):
        if list[o]>cvmax: 
            cvmax=list[o]
            index=o
    tc_values.append(TS[index]) #adds the critical point of the adapted function to a list of critical temperatures
    axes.plot(TS, list, color=colors[i]) #plots the function dapted to the data points 
    axes.plot(TS[index], list[index], "X", color="#000000") #plots the critical points
    axes.scatter(t_zoom, n[1:], label=f"L={size[i]}", color=colors[i]) #plots the data as scattered single points "o"
LS=(1/np.array([40 , 60 , 80 , 100, 120])) #sorts the size-array after redefining each point 1/size
tc_values=(np.array(tc_values)) #sorts the critical temperature values 
slope, intercept, r, p, se = stats.linregress(LS, np.array(tc_values)) #creates a linear regression for the critical temperature-values
#x=np.linspace(0, LS[-1], 10000) #defines an x-intervall from 0 to the max 1/size value in the LS array 
TC_value = 2.269185 #expected value for critical temperaure
rel_error = (intercept-TC_value)/(TC_value) #calculates error in each point
axes.set_title(r"Expected Heat Capacity, zoomed")
axes.set_xlabel(r"T $[J/k_B]$")  
axes.set_ylabel(r"$C_V [k_B]$")
axes.legend()
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../report/figures/heat_capacity_zoom.pdf")
plt.show()


#===========================================================================================================
#-------------------------------- TC ESTIMATION FROM 10^5 MONTE CARLO CYCLES ------------------------------#
#===========================================================================================================
#creates a plot visulaizing the linear regression for the critical temperatures
plt.style.use("ggplot")
xs=[]
fig, axes = plt.subplots(1, 1, figsize=(8, 6))
for i in range(len(size)):
    xs.append(1/size[i])
x=np.linspace(0,max(xs), 10000)
slope, intercept, r, p, se = stats.linregress(LS, np.array(tc_values)) #creates a linear regression for the critical temperature-values
axes.plot(x, intercept+slope*x, label=rf"Fitted line, $T_C=$ {intercept:.3f} $\pm$ relative error = {abs(rel_error):.3f}")
axes.plot(LS, np.array(tc_values), "o")
axes.set_xlabel("1/L")
axes.set_ylabel(r"Estmation of $T_C$, T $[J/k_B]$")
axes.set_title(r"Linear regression for estimation of $T_C$")
axes.legend()
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../report/figures/linear.pdf")
plt.show()
print("L:")
print(size)
print(r"T_C:")
print(np.array(tc_values))



#===========================================================================================================
#------------------------- 10^5 MOONTE CARLO CYCLES, ZOOMED INTO INTERVAL 2.0 - 2.5 -----------------------#
#------------------------- 10^5 MOONTE CARLO CYCLES, ZOOMED INTO INTERVAL 2.0 - 2.5 -----------------------#
#===========================================================================================================
plt.style.use("ggplot")
fig, axes = plt.subplots(2, 2, figsize=(10, 6))
size = [40 , 60 , 80 , 100, 120] #different size lattices
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77" , 'lightcoral', "slategray", 'lightgreen', 'cadetblue', "#117733", "#661100", "gold","forestgreen"]
t_array=np.arange(2.205 , 2.4 , 0.005) #temperatures
n=0 
cv_zoom=[] ; tc_values = []
for L in size: #loops through each lattice-size we modelled
    N=L*L #defines N as amount of particles in the lattice
    #defines empty arrays that will contain modelled expectation values E, E^2, M, M^2, CV
    E = np.zeros(len(t_array)) ; EE = np.zeros(len(t_array)) 
    M = np.zeros(len(t_array)) ; MM = np.zeros(len(t_array)) ; CV= np.zeros(len(t_array))
    zoom_cv=[] #defines a list for cv for temperatures within a special intervall
    for i in range(len(t_array)): #loops through each temperature the simulation is ran for 
        data = np.loadtxt(f"../data/mcc5zoom{L}_{t_array[i]:.3f}_exp.txt") #opens and reads the file where expectation values were stored
        E_= data[:,0]; EE_= data[:,1]; M_= data[:,2]; MM_= data[:,3]; abM_= data[:,4] ; cycles= data[:,5] #extracts columns of data
        #defines expectation values as the last calculated value by our simulation, since this is as close to critical temperature that we got
        E[i] = E_[-1]
        M[i] = abM_[-1]
        EE[i] = EE_[-1]
        MM[i] = MM_[-1]
        CV[i] = (EE[i]-E[i]**2)/(N*t_array[i]*t_array[i]) #defines heat capacity as by formula (7) in the report
        if t_array[i]>=2.20 and t_array[i]<=2.35: #we want to create a seperate zoomed in plot, which is done by testing if the temperature is within a special intervall
            zoom_cv.append((EE[i]-E[i]**2)/(N*t_array[i]*t_array[i]))      #if it is, the heat capacity value is added to the list  
        i+=1
    cv_zoom.append(zoom_cv) # the entire list of cv values in our special intervall is added to another list
    #we create an function that is adapted as well as possible to the modelled data for cv by using pol.UnivariateSpline
    cs = pol.UnivariateSpline(t_array,CV) 
    TS=np.linspace(2.0, 2.5, 10000) #array containing 10000 points in our temperature intervall, to make function more accurate
    tc_values.append(TS[np.argmax(cs(TS))]) #appends the temperature that has the same index as the modelled function for cv's apex to a list. this is where the model's critical point is
    #for each size lattice, the data for expectation values E, M, heat capacity and susceptibility is plotted as scatter-plots in their own subplots, with labels explaining what size lattice the data is extracted from
    axes[0 , 0].scatter(t_array[1:], E[1:]/N, label=f"L={L}", color=colors[n])
    axes[1 , 0].scatter(t_array[1:], M[1:]/N, label=f"L={L}", color=colors[n])
    axes[0 , 1].scatter(t_array[1:], (EE[1:]-E[1:]**2)/(N*t_array[1:]*t_array[1:]), label=f"L={L}", color=colors[n])
    axes[1 , 1].scatter(t_array[1:], (MM[1:]-M[1:]**2)/(N*t_array[1:]), label=f"L={L}", color=colors[n])
    axes[0 , 0].set_title(r"Expected Energy")
    axes[1 , 0].set_title(r"Expected Magnetization")
    axes[0 , 1].set_title(r"Expected Heat Capacity")
    axes[1 , 1].set_title(r"Expected Suceptibility")
    axes[0 , 0].set_xlabel(r"T $[J/k_B]$")
    axes[1 , 0].set_xlabel(r"T $[J/k_B]$")
    axes[0 , 1].set_xlabel(r"T $[J/k_B]$")
    axes[1 , 1].set_xlabel(r"T $[J/k_B]$")
    axes[0 , 0].set_ylabel(r"E [J]")
    axes[1 , 0].set_ylabel(r"$\langle |M| \rangle $")
    axes[0 , 1].set_ylabel(r"$C_V [k_B]$")
    axes[1 , 1].set_ylabel(r"$\chi [1/J]$")
    axes[0 , 0].legend()
    axes[1 , 0].legend()
    axes[0 , 1].legend()
    axes[1 , 1].legend()
    n+=1
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../report/figures/sizevariationszoom4.pdf") #saves the figure as pdf
plt.show()