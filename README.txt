README File for POD-Schrodinger Equation Solver
2D POD Quantum Dot Solver 
Authors: Martin Veresko and Min-Cheng Cheng

Please cite the following papers:
* M.C. Cheng, A reduced-order representation of the Schrödinger equation, AIP Advances 6, 095121 (2016) https://doi.org/10.1063/1.4963835 
* M. Veresko & M.C. Cheng, Physics-informed Reduced-Order Learning from the First Principles for Simulation of Quantum Nanostructures, (2023). arXiv preprint arXiv:2302.00100.

The following code will solve the Time Independent Schrodinger equation
using POD. Make sure in your directory in which the code is installed you have all files located https://github.com/CompResearchLab/POD-Schrodinger-Equation-Solver. If you are missing a file or directory from this site the code will not run. The code is divided into three parts: Schro_Solver_QDs.m, POD_Library_Creation.m
and POD_Schro_Solver.m. This first part, Schro_Solver_QDs.m, performs direct numerical simulations (DNSs) to train/collect 
wave function data. The second part, POD_Library_Creation.m generates the POD library. Finally,
POD_Schro_Solver.m will solve the POD Hamiltonian and will plot the error with 
regards to direct numerical simulation.

To run the entire code first you must gather wave function data to train the POD modes.
Hence, you must run Schro_Solver_QDs.m. Next you need to generate the POD modes and the library. Therefore, 
you need to run POD_Library_Creation.m. Finally, to solve the POD solution for a
particular case, you need to run POD_Schro_Solver.m. POD_Schro_Solver.m will put all graphics and results in the 
folder titled Results. Every time you run a new case the data in this folder will be deleted; hence, make sure
you save all data you want to keep from a previous test in another folder. POD_Schro_Solver.m will produce contour
plots, profile plots, a least square error plot, a table of the state energies and a recording of the computational 
time saved via the POD method.

Procedure to run code from start to finish:
1) Run Schro_Solver_QDs.m
2) Run POD_Library_Creation.m
3) Run POD_Schro_Solver.m (You can continue to run this code after training the modes to check different cases)

Usage: The user should only modify Training_Parameters.m and POD_Solver_Parameters.m. Training_Parameters.m 
allows the user to adjust the structure of interest and training parameters.
The file POD_Solver_Parameters.m allows the user to adjust POD-space parameters and plotting
guidelines. Information pertaining to each input is specified as a comment.

Note:
Do not delete the Library folder. Training data, stored modes and POD matrices are grabbed from here.
Sometimes the Schro_Solver_QDs.m program will print everything to the display. Simply running clear will correct this issue.




