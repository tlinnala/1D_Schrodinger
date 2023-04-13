import numpy as np

# Return the cosine potential
def cosine(y:list,n_atoms:int):
    return np.cos(2*np.pi*y*(n_atoms-1))/2+1/2

# Returns the position of the i:th nucleus in an array of n nuclei
# while accounting for the surfaces
def pos(i:int,n_atoms:int,s:float):
    return i*(1-2*s)/(n_atoms-1)+s

# Returns the KP potential
def KP(y:list,n_atoms:int,d:float):
    N=len(y)
    U=np.zeros(N)
    for i in range(N):
        for j in range(n_atoms):
            if abs(pos(j,n_atoms,0)-y[i])<=d/2:
                U[i]+=1
    if max(U)>1:
        raise Exception("Potentials overlap, check input values")
    return U

# Returns the KP potential with surfaces
def KP_surfaces(y:list,n_atoms:int,d:float,s:float,U_vac:float):
    N=len(y)
    U=np.zeros(N)
    l=(1-2*s-(n_atoms-1)*d)/(n_atoms-1)+d/2
    for i in range(N):
        for j in range(n_atoms):
            if abs(pos(j,n_atoms,s)-y[i])<=d/2:
                U[i]+=1
        if y[i]<pos(0,n_atoms,s)-l:
            U[i]+=U_vac
        elif y[i]>pos(n_atoms-1,n_atoms,s)+l:
            U[i]+=U_vac
    if max(U)>1:
        raise Exception("Potentials overlap, check input values")
    return U