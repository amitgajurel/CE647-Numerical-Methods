from sys import path
import yaml
import glob
import os
import numpy as np
from utils.functions import stiffnessmatrix as sm
from utils.functions import custom_plot 
import matplotlib.pyplot as plt
import csv
from scipy.interpolate import interp1d


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
su0 , k , gamma, e50, alpha,Nc,\
    d0, t, l , E, yld, \
        nx, tol,\
        Qn, Z1n, tn, Z2n,\
             BC, delta_force, head_cond, tip_cond, plug_cond = list(config.values())

# Calculation of derived quantities of the piles
A = np.pi/4 * (d0**2 - (d0-2*t)**2)                 # Crosssectional Area of the pile in in^2
A0 = np.pi/4 * (d0**2)                              # Crossection area of the pile for plugged case in in^2
EA = E*A                                            # Output units in pounds               


# Initializing the value of Tmax in the the t-z curve equation along all the nodes
dx = l/nx                                           # Total number of elements
num = nx + 1                                        # Total number of nodes
x = np.linspace(0,l,num)                            # Divides the "l" into "num" nodes starting from 0
su = []
tmax = []

for idx,val in enumerate(x):                        # populating each node with su and tmax value
    su.append(su0 + k * val)
    tmax.append(alpha * su[idx] * np.pi * d0/12)    # all units in ft


# Initializing intial displacement field with triangular displacement field with max of 4 in of displacement
deltamin=0.002                                  # Initial minimum displacement at the bottom of the pile
deltamax = 0.4                                 # Initial maximum displacement at the top of the pile
z = -1*np.linspace(deltamax,deltamin,num)

# initializing the t values for calculating given set of displacements based on the given API t-z curves
t = np.zeros(num)

# initializing ks - slope of t-z curve values for each and every node.
ks = np.zeros(num)

# initialzing the Q value for q-z curve for each and every node
kp = np.zeros(num)

# difference between predicted from the inversion and the t-z curve
z_diff = 1

# Running while loop until the difference between y that statifies the boundary condition and the P-y curve is met
loop_count = 0 

while z_diff > tol:
    # Working with the API t-z curves to fill up the matrix of ks
    intp_tn = interp1d(x=Z2n,y=tn,kind='linear')(np.abs(z/(d0/12)))             # value of interpolated tn i.e. y-axis of t-z API curve
    
    for idx,val in enumerate(tmax):
        if np.abs(z[idx])>0.00000001:
            tmob = intp_tn[idx] * tmax[idx] * np.sign(z[idx])                       # using sign to preserve the direction
            ks[idx] =  np.abs(tmob * np.pi * d0/12) / np.abs(z[idx])                # Stiffness in terms of Force/Area
        else:
           ks[idx] = 0
    
    # Working witht the API q-z curve
    intp_qn = interp1d(x=Z1n, y=Qn,kind='linear')(np.abs(z[num-1]/(d0/12)))
    qmax = su[num-1] * Nc * A0
    qmob = intp_qn * qmax * np.sign(z[num-1])                                   # ensuring direction is preserved for sign
    
    [k, R] = sm(num,dx,ks,z,qmob,delta_force, EA, BC,tip_cond, plug_cond)

    z_solved = np.dot(np.linalg.inv(k),R)
    z_diff = np.max(np.abs(z-z_solved))/(d0/12)
    z = z_solved
    loop_count+=1
    if loop_count > 1000:
        break

### Post Processing for plotting and exporting the results

# Depth Vs Displacement
depth = np.linspace(0,l,num)                                             # Creating a plotting value along the depth
displacement = [x for x in np.abs(z)]*12                             # Changing the list of one element arrary to float and multiplying to get disp. in in.                                                    # Changing the ft to in. to display the displacement
plt = custom_plot(x=z*-12,y=depth,xlabel="Displacement(in.)",
                  ylabel="Depth (ft)",title="Displacement Vs Depth ")   
plt.savefig('../output/DisplvsDepth.jpg')                               # Saving figure outputs

# Writing a CSV file that was used to plot the above figure
with open('../output/DisplvsDepth.csv','w') as csvfile:
    writer = csv.writer(csvfile,delimiter=',')
    writer.writerow(['Depth (Ft)', 'Displacement (in)'])
    writer.writerows(zip(depth,displacement))

z = np.concatenate(z,axis = 0)                        
SR = np.multiply(ks, np.abs(z))     # Skin Friction at every point

# Creating a plotting value along the depth
plt = custom_plot(x=SR,y=depth,xlabel="Mobilized Skin Friction (lb)",
                  ylabel="Depth (ft)",title="Moment Vs Depth ")   
plt.savefig('../output/MomentvsDepth.jpg')                               # Saving figure outputs

print("This is the end")



