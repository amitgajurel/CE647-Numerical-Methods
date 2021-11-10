from sys import path
import yaml
import glob
import os
import numpy as np
from utils.functions import stiffnessmatrix as sm
from utils.functions import custom_plot 
import matplotlib.pyplot as plt
import csv


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
I = np.pi/64 * (d0**4 - (d0-2*t)**4)                # moment of inertia for circular hollow pile
S = 2 * I/d0                                        # section modulus
EI = E*I/(12**2)                                    # Changing the unit to lbs.ft2

# reference deflection
y50 = 2.5 * e50 * d0 / 12                           # Unit from here onwards is ft.

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
y = np.linspace(0.5,0.3,num)

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
        kmax = pn[1]/yn[1] * Pmax[idx]/y50 *10                                               # Limiting the value of maximum stiffness
        if val < 8:
            P[idx] = 0.5 * np.abs(val)**(1/3) * Pmax[idx]
            temp_ks = P[idx]/np.abs(np.float128(val*y50))                                # Make sure that Ks is positive when iterating
            if temp_ks < kmax:
                ks[idx] = temp_ks
            else:
                ks[idx] = kmax        
        else:
            P[idx] = Pmax[idx]
            temp_ks = P[idx]/np.abs(np.float128(val*y50))         
            if temp_ks < kmax:
                ks[idx] = temp_ks
            else:
                ks[idx] = kmax 


    [k, R] = sm(BC,num,dx,ks,delta_force, head_cond, tip_cond,EI)

    y_solved = np.float128(np.dot(np.linalg.inv(k),R))
    y_diff = np.max(np.abs(y-y_solved))
    y = y_solved
    loop_count+=1

### Post Processing for plotting and exporting the results

# Depth Vs Displacement
depth = np.linspace(0,l,num)                                             # Creating a plotting value along the depth
displacement = [np.float64(x) for x in y]*12                             # Changing the list of one element arrary to float and multiplying to get disp. in in.                                                    # Changing the ft to in. to display the displacement
plt = custom_plot(x=y*12,y=depth,xlabel="Displacement(in.)",
                  ylabel="Depth (ft)",title="Displacement Vs Depth ")   
plt.savefig('../output/DisplvsDepth.jpg')                               # Saving figure outputs

# Writing a CSV file that was used to plot the above figure
with open('../output/DisplvsDepth.csv','w') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    writer.writerow(['Depth (Ft)', 'Displacement (in)'])
    writer.writerows(zip(depth,displacement))



BM = np.zeros(num)
SF = np.zeros(num)
SR = np.zeros(num)

# Depth, Shear Force, Soil resistance Vs Bending moment
for idx,val in enumerate(np.arange(0,num-4)):
    # using forward difference for 0 to num-4. i.e. leaving last 4 nodes
    BM[val] = EI * (y[val]-2*y[val+1]+y[val+2])/dx**2
    SF[val] = EI * (-y[val]+3*y[val+1]-3*y[val+2]+y[val+3])/dx**3
    SR[val] = EI * (y[val]-4*y[val+1]+6*y[val+2]-4*y[val+3]+y[val+4])/dx**4

for idx,val in enumerate(np.arange(num-4,num)):
    # using backward difference for num-4 to num. i.e. including last 4 nodes
    BM[val] = EI * (y[val]-2*y[val-1]+y[val-2])/dx**2
    SF[val] = EI * (y[val]-3*y[val-1]+3*y[val-2]-y[val-3])/dx**3
    SR[val] = EI * (y[val]-4*y[val-1]+6*y[val-2]-4*y[val-3]+y[val-4])/dx**4

depth = np.linspace(0,l,num)                                             # Creating a plotting value along the depth
plt = custom_plot(x=BM,y=depth,xlabel="Moment(ft-lb)",
                  ylabel="Depth (ft)",title="Moment Vs Depth ")   
plt.savefig('../output/MomentvsDepth.jpg')                               # Saving figure outputs

# Writing a CSV file that was used to plot the above figure
with open('../output/MomentvsDepth.csv','w') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    writer.writerow(['Depth (Ft)', 'Moment (Ft-lb)'])
    writer.writerows(zip(depth,BM))


plt = custom_plot(x=SF,y=depth,xlabel="Shear(lbs.)",
                  ylabel="Depth (ft)",title="Shear Vs Depth ")   
plt.savefig('../output/ShearvsDepth.jpg')                               # Saving figure outputs

# Writing a CSV file that was used to plot the above figure
with open('../output/ShearvsDepth.csv','w') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    writer.writerow(['Depth (Ft)', 'Shear (lbs.)'])
    writer.writerows(zip(depth,SF))


plt = custom_plot(x=SR*-1,y=depth,xlabel="Soil Resistance(lbs/ft.ft)",
                  ylabel="Depth (ft)",title="Soil Resistance Vs Depth ")   
plt.savefig('../output/SRvsDepth.jpg')                               # Saving figure outputs

# Writing a CSV file that was used to plot the above figure
with open('../output/SRvsDepth.csv','w') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    writer.writerow(['Depth (Ft)', 'Shear (lbs.)'])
    writer.writerows(zip(depth,SR*-1))

print("This is the end")



