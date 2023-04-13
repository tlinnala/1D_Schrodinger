import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.constants import hbar,m_e
import matplotlib.pyplot as plt
import potential_functions as pot

N=1000                # Number of discretization points
n_energies=25         # Number of energies to print
n_atoms=9             # Number of nuclei
l=2*10**-10           # Period (m)
L=(n_atoms-1)*l       # Number of periods in interval
h=1/(N-1)             # Step size
a=1/(h**2)            # Parameter a
b=-1/(2*h**2)         # Parameter b
la=5000               # Lambda parameter
d=0.04                # Range of of the KP potential
s=0.2                 # Length of the surface
U_vac=0.25            # Surface potential

# Here k is used to convert energy units to eV
k=6.24151*10**18*hbar**2/(m_e*L**2)

# Define the variable range
y=np.linspace(0,1,num=N)

# Compute the values of the potential function using the potential_functions module
# Comment/uncomment the following lines to use different funtions
#V_discrete=pot.cosine(y,n_atoms)
#V_discrete=pot.KP(y,n_atoms,d)
V_discrete=pot.KP_surfaces(y,n_atoms,d,s,U_vac)

# Plot and show the potential
# You should check that the potential looks as intended
plt.plot(y,V_discrete)
plt.title("Potential function")
plt.xlabel("Position")
plt.ylabel("Potential")
#plt.savefig("potential")
plt.show()

# Format lists A and B which contain the diagnoal
# and off-diagonal elements of the matrix
A=a*np.ones(N)+[la*i for i in V_discrete]  # Diagonal elements
B=b*np.ones(N-1)                           # Off-diagonal elements

# Solve the energies and the states
energy,state=eigh_tridiagonal(A,B,select="a")

# Plot the i:th eigenstate T[i]
i=0 # Ground state by default
plt.plot(state.T[0]**2,label="$\psi^2$ for state "+str(i+1))
plt.title("Electron probability amplitude")
plt.xlabel("Position ($h$)")
plt.ylabel("Probability amplitude")
plt.legend()
#plt.savefig("probability_amplitudes")
plt.show()

# Plot the eigenvalues and print five lowest eigenvalues
plt.scatter(np.arange(0,n_energies,1)+1,k*energy[0:n_energies],s=20,marker='.')
print("Lowest energies (eV): ",k*energy[0:5])

plt.title("Energies of an electron")
plt.xlabel("$n$")
plt.ylabel("Energy (eV)")
#plt.savefig("energies")
plt.show()