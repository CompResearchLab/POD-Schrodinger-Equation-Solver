%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PODSolverInput%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------Electric Fields to use---------------------%
% These are the two components of the electric field to verify the POD
% model with. 
XField=50; %[kV/cm] Field in the X direction 
YField=-10; % kV/cm Field in the Y direction
%------------------------------States-------------------------------% 
StatesToBeCollectedForPODComparision=8; %This is the number of states to verify the POD model with
%------------------------------Hamiltonian Set up-------------------------%
TotalModesInHamiltonian=30;% Total Modes in POD Hamiltonian Matrix.
% The Hamiltonian matrix is solved once and then solutions using N number
% of the modes are plotted and used to determine the LS error. Additionally,
% this number needs to be larger than the total number of states requested.
NumberofModesToPlot=20; %Total Number of Modes to consider in the LS Error plot.
% This number must be less than TotalModesInHamiltonian.
%----------------------------Plotting Commands----------------------------%
%UsedModes reveals the Mode solution to plot for each state. The number of
%rows needs to equal StatesToBeCollectedForPODComparision. Additionally, the number of
%columns needs to be equal for each state. To specify that no more Mode
%solutions are needed just input zeros.
UsedModes=[1,2,6,0;   1,4,8,0;   1,2,6,9;   1,4,8,0;  1,6,13,0;   1,6,9,0;   1,7,13,0;   1,8,13,0];

HeightAboveWave=.57 %[eV] % This is a plotting parameter to increase the y-axis of the profile plots.
% If WFs are being cut off on the top of the figure, increase this
% parameter.




