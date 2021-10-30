from sys import path
import yaml
import glob
import os
import numpy as np
from utils.functions import stiffnessmatrix as sm
import matplotlib.pyplot as plt


# Make sure your working directory is same as the scripts locations
path = os.path.dirname(__file__)

# Changing the CWD to the current working directory
os.chdir(path)

# loading the YAML files in the the data directory
files = glob.glob(pathname="../data/*.yml")
temp = open(files[0],"r")

# YAML loads the files a dictionary
config = yaml.full_load(temp)

# assigning values to the variables. Reading from the input.yml file
su0 , k , gamma, e50, d0, t, l , E, yld ,nx, tol, yn, pn, J, BC, delta_force, head_cond, tip_cond = list(config.values())

# Calculation of derived quantities of the piles
I = np.pi/64 * (d0**4 - (d0-2*t)**4)            # moment of inertia for circular hollow pile
S = 2 * I/d0                                    # section modulus
EI = E*I/(12**2)                                  # Changing the unit to lbs.ft2

# reference deflection
y50 = 2.5 * e50 * d0 / 12                       # Unit from here onwards is ft.

# Calcualtion of Pmax in the the P-y curve equation
dx = l/nx
num = nx + 1
x = np.linspace(0,l,num)
xi = []
su = []
Np = []
Pmax = []

for i in np.linspace(1,num,num):
    i = int(i)                                  # indexing can only be done with int not float
    xi.append((i - 1)*dx)
    su.append(su0 + k * xi[i-1])
    Np.append(3 + (gamma-62.4)/su[i-1]*xi[i-1] + J/(d0/12)*xi[i-1])  # all units in ft
    Pmax.append(Np[i-1] * su[i-1])

# limiting the maximum value of Np as 9
for idx, val in enumerate(Np):
    if val < 9:
        pass
    else:
        Np[idx] = 9

# Initializing intial displacement field with triangular displacement field with max of 4 in of displacement
y = np.linspace(0.1,0.0,num)

# initializing the P values for a given set of displacements based on the given API P-y curves
P = np.zeros(num)

# initializing ks - slope of P-y curve values for each and very node.
ks = np.zeros(num)

# difference between predicted from the inversion and the P-y curve
y_diff = 1

# Running while loop until the difference between y that statifies the boundary condition and the P-y curve is met
loop_count = 0 

while y_diff > tol:

    # calculating P and ks values for each and every node given the displacement field
    for idx,val in enumerate(y/y50):
        if val < 8:
            P[idx] = 0.5 * np.abs(val)**(1/3) * Pmax[idx]
            ks[idx] = P[idx]/np.abs(val*y50)                                # Make sure that Ks is positive when iterating
        else:
            P[idx] = Pmax[idx]
            ks[idx] = P[idx]/(val*y50)


    [k, R] = sm(BC,num,dx,ks,delta_force, head_cond, tip_cond,EI)

    y_solved = np.dot(np.linalg.inv(k),R)
    y_diff = np.mean(np.abs(y-y_solved))
    y = y_solved
    loop_count+=1

plt.plot(y,-1*np.linspace(1,num,num))
plt.savefig("Temp.jpg")
print("This is the end")


