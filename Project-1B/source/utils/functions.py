# this script contains necessary function for the main scripts

def stiffnessmatrix(num,dx ,ks,z,qmob,delta_Force ,EA ,BC = 0,tip_cond=0, plug_cond = 0):
    """
    INPUTS

    BC              : Boundary Condition If BC = 0 (disp. controlled) or 1 (force controlled)
    num             : Number of nodes i.e. number of segments + 1
    dx              : Size of each element
    ks              : Spring stiffness for side friction at all the internal points/nodes along the pile
    qmob            : Mobilited qz at the end bearing for the last node of the pile
    delta_force     : if BC = 0 provide the delta (in.) and if BC = 1 provide force (lbs.)
    tip_cond        : is equal to 0 for displacement possible at bottom ; 1 for zero displacement i.e. pile in rock
    plug_cond       : is equal to 0 for plugged pile
    z               : axial displacement profile
    EA              : Axial Stiffness
    
    RETURNS
    [k, R]          : Where k -> stiffness matrix and R -> Load Vector

    """
    import numpy as np

    # Creating a numpy arrary of zeros with size num x num
    k = np.zeros((num,num))
    R = np.zeros((num,1))

    # Central Nodes
    for idx,val in enumerate(np.arange(1,num-1)):
        # Using 2nd order finite difference approximation for axial load tranfer equation i.e. EA d2z/dx2 -ksz=0 
        k[val,val-1] = 1
        k[val,val] = -1 * (2 + (ks[val] * dx**2)/ EA)                # Spring axial stiffness term
        k[val,val+1] = 1

    # First node is used for force or displacement control BC at the tip of the pile
    if BC == 0:
        k[0,] = 0                                           # zero all
        k[0,0] = 1                                          # for displacement control
        R[0] = -1 * delta_Force/12                             # Changing the unit from input in. to ft.   
    elif BC == 1:
        # TODO: Need to check
        k[0,] = 0
        k[1,1] = 1
        k[1,0] = -1
        R[1] = -1 * delta_Force * dx / EA                        # Changing force to displacement
        
    # last row of node is used for type of boundary condition at the bottom of the pile
    
    if plug_cond == 0:
        # last row of node is used for type of boundary condition at the bottom of the pile
        # imposing plugged at the bottom i.e. there is resistance from the soil at the bottom of the pile
        k[num-1,num-1] = 1
        k[num-1,num-2] = -1
        R[num-1] = -qmob * dx / EA
        
    elif plug_cond == 1:
        # last row of node is used for type of boundary condition at the bottom of the pile
        # imposing no force at the boundary
        # TODO: Check
        pass
 

    return([k,R])

def custom_plot(x,y,xlabel,ylabel,title):
    """The function takes makes a pretty plot with labels for piles

    Args:
        x (array/list): List or array to plot in X-axis
        y (array/list): List or array to plot in Y-axis
        xlabel (string): String for Xlabel
        ylabel (string): Strimg fpr Ylabel
        title (string) : String for title
    Output:
        Plot object
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sb
    import pandas as pd
    import numpy as np

    fig, axs = plt.subplots(dpi=600)
    fig.set_size_inches(6,6)
    fig.tight_layout(pad=4)

    axs.xaxis.set_tick_params(which='major', size=5, width=0.7, direction='in', top='on')
    axs.xaxis.set_tick_params(which='minor', size=2.5, width=0.35, direction='in', top='on') 

    axs.yaxis.set_tick_params(which='major', size=5, width=0.7, direction='in', right='on')
    axs.yaxis.set_tick_params(which='minor', size=2.5, width=0.35, direction='in', right='on')

    axs.plot(x,y,'--')
    axs.invert_yaxis()

    axs.set_title(title,fontsize=10,fontweight='bold')
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.grid(linewidth='0.25')
    
    plt.ioff()
    return(fig)
