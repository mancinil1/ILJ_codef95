Description

The codes have been developed in order to calculate the interaction potential considering different models, such as the Harmonic, Lennard-Jones (LJ), Improved Lennard-Jones (ILJ) and Morse. The program starts from the boundaries, then integrates using an arbitrarily small quantity of integration step to the matching point. The initial parameters and boundary conditions are given by the input files. The trial energy is set to a value slightly larger than the well depth of the corresponding potential. First part, the program calculates the radial Schrödinger equation to get guess energies using trial energy from input file, by performing two integrations using the Numerov method. Second part, the program will evaluate the current and the previous energies using the shooting method. This method searches the energy, between those two energies, that gives first derivative difference (error) less than the tolerance value that we give in the input file. The program will reintegrate the radial Schrödinger equation using Numerov subroutine and this procedure will be repeated until it gets the best energy. This procedure is built to solve even quantum number states, the criteria in the shooting method is changed in reverse way if we deal with odd quantum numbers. The last part, after we find the best energy (En) (it means that the radial Schrödinger equations have been solved), the entire obtained wavefunctions will be normalized and saved in some output files.

How to cite 

Y. B. Apriliyanto, A. Lombardi, L. Mancini, F. Pirani, N. Faginas-Lago (Indonesia Defense University, Università degli Studi di Perugia, 2024)

https://github.com/mancinil1/ILJ_codef95
