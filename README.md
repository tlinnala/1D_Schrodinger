# A Solution to the one-dimensional and time-independent Schrödinger equation
This simple program can be used to solve the one-dimensional time-independent Schrödinger equation using a matrix method. The method is based on the fact that the differential equation can be discretised. The main program file eigval_prob contains the actual solution and functions used to format the potential are defined in the potential_functions module.

One can use the program to compute the electronic band structures of one-dimensional solids using different potentials. The relevant parameters to change are N (number of discretisation steps), n_energies (number of energies to print), n_atoms (number of atoms in the lattice), l (distance between adjacent atoms) and lambda (describes the strength of the potential) as well as the potential function used. Three potentials are defined in the module:

cosine(y,n_atoms): returns the cosine potential,

KP(y,n_atoms,d): returns the KP potential without surfaces,

KP_surfaces(y,n_atoms,d,s,U_vac): returns the KP potential with surfaces.

Try out different values of d, s and U_vac for different KP potentials. Be careful when adjusting these parameters, since it is possible to, for example, choose d such that the potentials of different nuclei overlap. The main program plots the resulting potentials, so check that it looks as intended. You can also define your own potential function.
