%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%DNS Schrodinger Equation Solver%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solves the Schrodinger Equation via a finite difference scheme for specified 
%amount of test electric. This program uses both Drichelet and Newman
%boundary conditions. 
%-------------Clear Command Windom/Workspace and Close Figures------------%
clc;
clear;
close all;
%----------------------------Read Input File-------------------------------%
Training_Parameters;
%Variables Include:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eff1:Effective mass outside of well
%Eff2:Effective mass inside the well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dxDistance: Grid spacing in nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Depth: Conduction Band Offset of QDs in eV
%Borderthickness: Thickness of the boarder surrounding the lattice
%structure in nm
%TotalHeight: Total height of a single quantum well subunit in nm (inside and
%outside)
%TotalWidth: Total width of a single quantum well subunit in nm (inside and
%outside)
%WellHeight: Total height of the inside of the well (must be smaller than
%the total height)
%WellWidth: Total width of the inside of the well (must be smaller than the
%total width of the well)
%Rows: Times the quantum well structure will be repeated along the y-direction
%Columns: Times the quantum well structure will be repeated along the y-direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MaxElectricFieldPerDirection: Maximum Magnitude of Tests kV/cm
%NumberOfSamplesPerdirection: %number of samples per direction
%NumberSavedStates:Total Number of Saved States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%StatesToBeCollectedForPODComparision: Number of States To be Collected To
%Test POD model
%-----------------------Constants-----------------------------------------%
h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
q_el = 1.60218e-19; % elementary charge, in coulombs
Mass=9.10938215e-31; % Rest mass of particle in kg
%------------------------------Creation of Test Potentials----------------%
MaxElectricFieldPerDirection=MaxElectricFieldPerDirection*100000; % convert to eV/m
%Determine electric field magnitudes
PotentialIncrement=linspace(0,MaxElectricFieldPerDirection,NumberOfSamplesPerdirection+1);
%Determine number of total samples used
SampleNumber=(NumberOfSamplesPerdirection)*4+1; % +1 is from the unbiased test
%PotentialCoefficients: Stores every test potential used. First Column is
%the x-component and the second column is the y-component
PotentialCoefficients=zeros(SampleNumber,2);
PotentialCoefficients(2:NumberOfSamplesPerdirection+1,2)=PotentialIncrement(2:end)';
PotentialCoefficients((NumberOfSamplesPerdirection+2):(2*NumberOfSamplesPerdirection+1),2)=-PotentialIncrement(2:end)';
PotentialCoefficients(2*NumberOfSamplesPerdirection+2:3*NumberOfSamplesPerdirection+1,1)=PotentialIncrement(2:end)';
PotentialCoefficients(3*NumberOfSamplesPerdirection+2:end,1)=-PotentialIncrement(2:end)';
%------------------------Convert Input spacial variables to meters--------%
dx=dxDistance*(10^(-9));
dy=dxDistance*(10^(-9));
TotalHeight=TotalHeight*10^(-9);
TotalWidth=TotalWidth*10^(-9);
WellHeight=WellHeight*10^(-9);
WellWidth=WellWidth*10^(-9);  
Borderthickness=Borderthickness*(10)^(-9);
%Break out of the program if spatial grid is not completely divisible by dx.
if mod(2*Borderthickness+TotalHeight,dx)~=0||mod(2*Borderthickness+TotalWidth,dx)~=0
    disp(" Lattice parameters must be divisble by dx and dy");
    return;
end
%---------------------------Create Single QD substructure-----------------%
conditionx=(TotalWidth-WellWidth)/2;
conditiony=(TotalHeight-WellHeight)/2;
[XCELL,YCELL]=meshgrid(0:dx:TotalWidth,0:dy:TotalHeight);
UCell=Depth*~((conditiony<YCELL&YCELL<conditiony+WellHeight)&(conditionx<XCELL&XCELL<conditionx+WellWidth));
%--------------------------Form Rows and Columns of QDs-------------------%
U=repmat(UCell(2:end,2:end),Rows,Columns);% Repeat QD for Rows and Columns
%The zero x-coordinates and y-coordinates need to be remove before
U=[Depth*ones(size(U,1),1),U]; %Put potential located at x=0 back
U=[Depth*ones(1,size(U,2));U]; %put back zero at y=0 back
%-----------------------------Add Boarders--------------------------------%
XBorderthickness=round(Borderthickness/dx);
YBorderthickness=round(Borderthickness/dy);
U=[Depth*ones(size(U,1),XBorderthickness),U];
U=[U,Depth*ones(size(U,1),XBorderthickness)];
U=[Depth*ones(YBorderthickness,size(U,2));U];
U=[U;Depth*ones(YBorderthickness,size(U,2))];
%-----------------------------Create Grid--------------------------------%
y=[-Borderthickness:dy:TotalHeight*Rows+Borderthickness];
x=[-Borderthickness:dx:TotalWidth*Columns+Borderthickness];
y=y-y(1);
x=x-x(1);
[X,Y]=meshgrid(x,y);
%Get number of rows and columns of grid with boundaries removed. 
[rw,col]=size(U(3:end-2,3:end-2));
%Number of Rows and Columns with boundaries
ROWGRID=rw+4;
COLGRID=col+4;
%------------------------Determine Effective Mass for each Coordinate------%
EffectiveMassMatrix=zeros(size(U,1),size(U,2));
%Assign an effective Mass to each coordinate of the grid. If the potential
%at a grid position is equal to the Depth set the effective mass to
%Eff1;elsewise set the effective mass to Eff2. 
for q=1:1:size(U,2)
    for k=1:1:size(U,1)
               if(U(k,q)==Depth) 
                    EffectiveMassMatrix(k,q)=Eff1; %
                else
                    EffectiveMassMatrix(k,q)=Eff2;
                end 
    end
end
%--------------------------Average Effective Mass-------------------------%
RemovedEffectiveMass=EffectiveMassMatrix(3:end-2,3:end-2);% Remove 2 boundaries. 
EffectiveMassLeftAverage=Eff1*ones(rw,col);
EffectiveMassRightAverage=Eff1*ones(rw,col);
EffectiveMassUpperAverage=Eff1*ones(rw,col);
EffectiveMassLowerAverage=Eff1*ones(rw,col);
EffectiveMassCenter=zeros(rw,col);
%Get Right Averaged Effective Mass
for z=1:rw
   for q=1:col-1
      EffectiveMassRightAverage(z,q)=(RemovedEffectiveMass(z,q)+RemovedEffectiveMass(z,q+1))/2;
   end
end
%Get Left Averaged Effective Mass
for z=1:rw
   for q=2:col
      EffectiveMassLeftAverage(z,q)=(RemovedEffectiveMass(z,q)+RemovedEffectiveMass(z,q-1))/2;
   end
end
%Get Upper Averaged Effective Mass
for z=1:col
   for q=2:rw
      EffectiveMassUpperAverage(q,z)= (RemovedEffectiveMass(q-1,z)+RemovedEffectiveMass(q,z))/2;
   end
end
%Get Lower Averaged Effective Mass
for z=1:col
   for q=1:rw-1
      EffectiveMassLowerAverage(q,z)= (RemovedEffectiveMass(q+1,z)+RemovedEffectiveMass(q,z))/2;
   end
end
%Get Center Averaged Effective Mass
for z=1:col
   for q=1:rw
      EffectiveMassCenter(q,z)=1/EffectiveMassLowerAverage(q,z)...
          +1/EffectiveMassUpperAverage(q,z)+1/EffectiveMassLeftAverage(q,z)...
          +1/EffectiveMassRightAverage(q,z);
   end
end
%---------------------Create Hamiltonian Matrix Without Potential---------%
%Invert Average Effective Masses
%Invert Average Effective Masses
t=(h_bar)^2/(2*Mass*(dx)^2)/(q_el); % h^2/2ma in eV
EffectiveMassRightAverage= EffectiveMassRightAverage.^(-1);
EffectiveMassLeftAverage=EffectiveMassLeftAverage.^(-1);
EffectiveMassUpperAverage=EffectiveMassUpperAverage.^(-1);
EffectiveMassLowerAverage=EffectiveMassLowerAverage.^(-1);
%Set Matrices to zero when wave function is zero. 
EffectiveMassRightAverage(:,col)=0;
EffectiveMassLeftAverage(:,1)=0;
EffectiveMassUpperAverage(1,:)=0;
EffectiveMassLowerAverage(rw,:)=0;
%Convert the EffectiveMass Matrix to vectors
Centervector=ToVector(EffectiveMassCenter,rw,col);
B=sparse([1:rw*col],[1:rw*col],t*Centervector,rw*col,rw*col);
RightVector=ToVector(EffectiveMassRightAverage(1:end,1:col-1),rw,col-1);
B=B+sparse(1:rw*(col-1),rw+1:rw*col,-t*RightVector,rw*col,rw*col);
LeftVector=ToVector(EffectiveMassLeftAverage(:,2:end),rw,col-1);
B=B+sparse(rw+1:rw*col,1:rw*(col-1),-t*LeftVector,rw*col,rw*col);
UpperVector=ToVector(EffectiveMassUpperAverage(2:end,:),rw-1,col);
RowIndices=2:rw*col;
ColIndices=1:rw*col-1;
RowIndices(rw:rw:end)=[];
ColIndices(rw:rw:end)=[];
B=B+sparse(RowIndices,ColIndices,-t*UpperVector,rw*col,rw*col);
LowerVector=ToVector(EffectiveMassLowerAverage(1:end-1,:),rw-1,col);
RowIndices=1:rw*col-1;
ColIndices=2:rw*col;
RowIndices(rw:rw:end)=[];
ColIndices(rw:rw:end)=[];
B=B+sparse(RowIndices,ColIndices,-t*LowerVector,rw*col,rw*col);
clear Centervector RightVector LeftVector UpperVector LowerVector EffectiveMassRightAverage...
EffectiveMassLeftAverage EffectiveMassUpperAverage EffectiveMassLowerAverage EffectiveMassCenter
%-------------------------Get WF and Energies for each Potential----------%
%It: Iteration Number; Loop Stops when the Schrodinger equation for every 
%Potential Has been solved
for It=1:SampleNumber
    %Get Total Potential For each iteration
    UwPot=U+PotentialCoefficients(It,1)*X+PotentialCoefficients(It,2)*Y;
    RemoveBoundaries=UwPot(3:end-2,3:end-2);
    %Solve Schrodinger Equation and get Wave Functions and Energies
    [WaveVectors,WaveEnergy]=eigs(B+sparse(1:rw*col,1:rw*col,ToVector(RemoveBoundaries,rw,col),rw*col,rw*col),NumberSavedStates,'sm');
    clear UwPot
    clear RemoveBoundaries
    %Sort energies and wave functions
    WaveEnergy=diag(WaveEnergy);
    [~,Index]=sort(WaveEnergy);
    WaveVectors=WaveVectors(:,Index);
    %Put back boundaries and normalize WF
    NewV=zeros(ROWGRID*COLGRID,NumberSavedStates);
    for i=1:1:NumberSavedStates
         Z=[zeros(2,COLGRID);...
           zeros(rw,2),Matrix(WaveVectors(:,i),rw,col),zeros(rw,2);...
           zeros(2,COLGRID)];
         Z=normalize(Z,x,y)*Z; 
        NewV(:,i)=ToVector(Z,ROWGRID,COLGRID);   
    end
    SampleWFs(:,:,It)=NewV;
end
%-----------------------Save Data-----------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save all Test Samples, Number of States Saved Per sample, Total number of
%Samples
save('Library/SavedSamples','SampleWFs','NumberSavedStates','SampleNumber');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save Grid information
save('Library/Grid','X','Y','dx','dy','ROWGRID','COLGRID','x','y');
%Save the unbiased potential of the nanostructure
save('Library/TrainingPotentialandMass','U','EffectiveMassMatrix')
%-------------------Function To convert Matrix TO a vector----------------%
% Function to turn a Matrix to a vector
function Vector=ToVector(Matrix,row,col)
    Vector=zeros(row*col,1);
    I=1;
       for k=1:col;
           for q=1:row;
            Vector(I)=Matrix(q,k);
            I=I+1;
           end
       end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%