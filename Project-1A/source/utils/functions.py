# this script contains necessary function for the main scripts

def stiffnessmatrix(BC = 0,num,dx,ks,delta_Force,head_cond,EI):
    """
    BC              : Boundary Condition If BC = 0 (disp. controlled) or 1 (force controlled)
    num             : Number of nodes i.e. number of segments + 1
    dx              : Size of each element
    ks              : Spring stiffness
    delta_force     : if BC = 0 provide the delta (in.) and if BC = 1 provide force (lbs.)
    head_condition  : is equal to 0 for free head and 1 for fixed head
    EI              : Flexural Stiffness
    
    """
    import numpy as np

    # Creating a numpy arrary of zeros with size num x num
    k = np.zeros((num,num))
    R = np.zeros((num,1))

    # Central Nodes
    for idx,val in enumerate(np.arange(2,num-2,)):
        # Using 4th order beam bending problem with central difference 
        k[val,val-2] = 1
        k[val,val-1] = -4
        k[val,val] = # Spring stiffness term
        k[val,val+1] = -4
        k[val,val+2] = 1

    # First node is used for force or displacement control BC
    if BC == 0:
        k[0,] = 0               # zero all
        k[0,0] = 1              # for displacement control
        R[0,0] = delta_Force    
    elif BC == 1:
        # Third order forward difference of shear force
        k[0,] = 0
        k[0,0] = -1             # -ve sign is to account for the direction of force and direction of displacement
        k[0,1] = 3
        k[0,2] = -3
        k[0,3] = 1
        R[0,0] = delta_Force * dx**3 / EI
    
    # Second row of nodes is used for type of boundary condition
    if head_cond == 0:
        # Free head -> Moment = 0 Using 2nd order forward difference 
        k[1,0] = 1
        k[1,1] = -2
        k[1,2] = 1
        R[1,0] = 0          # Free head has moment zero 
    if head_cond == 1:
        # Fixed head -> Slope = 0 Using 1st order forward difference
        k[1,0] = 1
        k[1,1] = -1
        R[1,0] = 0          # Fixed head has slope 0




    return(k)

    