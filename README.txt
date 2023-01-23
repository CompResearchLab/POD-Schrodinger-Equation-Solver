README File for POD-Schrodinger Equation Solver
2D POD Quantum Dot Solver 
Authors: Martin Veresko and Ming-Cheng Cheng

The following code will solve the Time Independent Schrodinger equation
using POD. The code is divided into three parts: Schro_Solver_QDs.m, POD_Library_Creation.m
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

References:

Jiang, Lin, et al. "An effective physics simulation methodology based on a 
data-driven learning algorithm." Proceedings of the Platform for Advanced 
Scientific Computing Conference. 2022.

Veresko, Martin, and Ming-Cheng Cheng. "An Effective Simulation Methodology 
of Quantum Nanostructures based on Model Order Reduction." 2021 International 
Conference on Simulation of Semiconductor Processes and Devices (SISPAD). IEEE, 2021.

Cheng, Ming-C. "A reduced-order representation of the Schrödinger equation." 
AIP Advances 6.9 (2016): 095121.

Cheng, Ming-C. "Quantum element method for quantum eigenvalue problems derived 
from projection-based model order reduction." AIP Advances 10.11 (2020): 115305.



