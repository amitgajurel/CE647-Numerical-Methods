# this script contains necessary function for the main scripts

def stiffnessmatrix(BC,num,dx,ks,delta_Force,head_cond,tip_cond,EI):
    """
    INPUTS

    BC              : Boundary Condition If BC = 0 (disp. controlled) or 1 (force controlled)
    num             : Number of nodes i.e. number of segments + 1
    dx              : Size of each element
    ks              : Spring stiffness for all the internal points/nodes along the pile
    delta_force     : if BC = 0 provide the delta (in.) and if BC = 1 provide force (lbs.)
    head_cond       : is equal to 0 for free head and 1 for fixed head
    tip_cond        : is equal to 0 for zero moment & zero shear at bottom i.e. floating pile; 1 for zero displacement and zero slope i.e. pile in rock
    EI              : Flexural Stiffness
    
    RETURNS
    [k, R]          : Where k -> stiffness matrix and R -> Load Vector

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
        k[val,val] = 6 + ks[val] * dx**4 / EI               # Spring stiffness term
        k[val,val+1] = -4
        k[val,val+2] = 1

    # First node is used for force or displacement control BC at the tip of the pile
    if BC == 0:
        k[0,] = 0                                           # zero all
        k[0,0] = 1                                          # for displacement control
        R[0,0] = delta_Force    
    elif BC == 1:
        # Third order forward difference of shear force
        k[0,] = 0
        k[0,0] = -1                                         # -ve sign is to account for the direction of force and direction of displacement
        k[0,1] = 3
        k[0,2] = -3
        k[0,3] = 1
        R[0,0] = delta_Force * dx**3 / EI
    
    # Second row of nodes is used for type of boundary condition at the top of the pile
    if head_cond == 0:
        # Free head -> Moment = 0 Using 2nd order forward difference 
        k[1,0] = 1
        k[1,1] = -2
        k[1,2] = 1
        R[1,0] = 0          # Free head has moment zero 
    elif head_cond == 1:
        # Fixed head -> Slope = 0 Using 1st order forward difference
        k[1,0] = 1
        k[1,1] = -1
        R[1,0] = 0          # Fixed head has slope 0

    # last two row of node is used for type of boundary condition at the bottom of the pile
    if tip_cond == 0:
        # last row of node is used for type of boundary condition at the bottom of the pile
        # imposing moment equal to zero at the tip/bottom of the pile
        k[num-1,num-3] = 1
        k[num-1,num-2] = -2
        k[num-1,num-1] = 1
        
        # second to the last row of node is used for type of boundary condition at the bottom of the pile
        # imposing zero shear at the bottom
        k[num-2,num-4] = -1
        k[num-2,num-3] = 3
        k[num-2,num-2] = -3
        k[num-2,num-1] = 1
    elif tip_cond == 1:
        # last row of node is used for type of boundary condition at the bottom of the pile
        # imposing displacement equal to zero at the tip/bottom of the pile
        k[num-1,num-1] = 1

        # second to last tow of node is used for slope is zero at the bottom
        k[num-2,num-4] = 1
        k[num-2,num-3] = -2
        k[num-2,num-2] = 1

    return([k,R])

    