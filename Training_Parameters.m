%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% README: This is the input file used to train the POD model. 
%Do not delete any variables, only adjust their value.
%-----------------------------Effective Mass Data%%-----------------------%
% These two parameters allow one to adjust the effective mass of the
% material inside and outside of each quantum well. 
%Effective Mass outside Wells
Eff1=.067;
%Effective Mass inside Wells
Eff2=.023;
%----------------------------Grid Size-----------------------------------%
% This parameter here will adjust the grid size. The Total size of the 
%nanostructure needs to be divisible by this grid size.
dxDistance=.1;%[nm] 
%---------------------------QD Substructure-------------------------------%
% The next inputs are the dimensions of the nanostructure. 
% The nanostructure is created by building blocks of a total height and
% width. Inside, located at the center of each block is a quantum well with
% a width and a height. After this block has been repeated for the
% specified number of rows and columns a boarder is added to the structure.
% Example:
% Basic building block. 

%-------------   ^ Total Height 
%-------------   |
%----ooooo----   | ^
%----ooooo----   | | Well Height 
%----ooooo----   | v 
%-------------   |
%-------------   v
%    <----> Well Width
%<------------> Total width

% Repeat structure for 2 columns and 3 rows
%--------------------------     
%--------------------------   
%----ooooo--------ooooo----    
%----ooooo--------ooooo----    
%----ooooo--------ooooo----  
%--------------------------   
%-------------------------- 
%--------------------------     
%--------------------------   
%----ooooo--------ooooo----    
%----ooooo--------ooooo----    
%----ooooo--------ooooo----  
%--------------------------   
%--------------------------
%--------------------------     
%--------------------------   
%----ooooo--------ooooo----    
%----ooooo--------ooooo----    
%----ooooo--------ooooo----  
%--------------------------   
%--------------------------
% Add boarder now:
%------------------------------
%------------------------------
%------------------------------     
%------------------------------   
%------ooooo--------ooooo------    
%------ooooo--------ooooo------    
%------ooooo--------ooooo------  
%------------------------------   
%------------------------------ 
%------------------------------     
%------------------------------   
%------ooooo--------ooooo------    
%------ooooo--------ooooo------    
%------ooooo--------ooooo------  
%------------------------------   
%------------------------------
%------------------------------     
%------------------------------   
%------ooooo--------ooooo------    
%------ooooo--------ooooo------    
%------ooooo--------ooooo------
%------------------------------  
%------------------------------

Depth=0.544;  %[eV] % Conduction band offset 
Borderthickness=2; %[nm] This is the thickness of the boarder added to the nanostructure
TotalHeight=5; %[nm] Total Height of the nanostructure building block.
TotalWidth=5;%[nm] Total Width of the nanostructure building block.
WellHeight=4;%[nm] Well Height inside the nanostructure building block.
WellWidth=4;%[nm] Well Width inside the nanostructure building block.
Rows=4; % Total number of rows to repeat the pattern.
Columns=4; % Total number of columns to repeat the pattern.
%---------------------------Test Electric fields and States collected-----%
% To train the POD modes, the nanostructure is subjected with electric
% fields. Here we train with single component electric fields. 
MaxElectricFieldPerDirection=35; %[kV/cm] Maximum magnitude of the training electric fields. 
NumberOfSamplesPerdirection=4; %number of fields per direction. (xhat, yhat,-xhat, -yhat).
NumberSavedStates=6;
%----------------------------POD library Creator---------------%
% Total modes for the POD kinetic energy and potential energy matrices to
% be stored in the libary. This number should be bigger that the number of
% modes one needs to use during application. 
LibraryModes=30;% Total Modes in POD Hamiltonian Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


