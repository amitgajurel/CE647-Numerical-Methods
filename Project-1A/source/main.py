from sys import path
import yaml
import glob
import os
import numpy as np

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
su0 , k , gamma, e50, d0, t, l , E, yld ,nx, tol, yn, pn, J, BC = list(config.values())

# Calculation of derived quantities of the piles
I = np.pi/64 * (d0**4 - (d0-2*t)**4)            # moment of inertia for circular hollow pile
S = 2 * I/d0                                    # section modulus
EI = E*I

# reference deflection
yc = 2.5 * e50 * d0

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
    Np.append(3 + (gamma-9.81)/su[i-1]*xi[i-1] + J/d0*xi[i-1])
    Pmax.append(Np[i-1] * su[i-1])

# limiting the maximum value of Np as 9
for idx, val in enumerate(Np):
    if val < 9:
        pass
    else:
        Np[idx] = 9

cvc = 2


