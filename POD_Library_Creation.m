%%%%%%%%%%%%%%%%%%%%%%%%%POD Mode Generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate POD modes from test samples
%-----------------Clear Workspace,Command Window and Close Figures--------%
clc
clear all;
close all;
%-------------------------------Load Data from DNS program----------------%
Training_Parameters;
load('Library/SavedSamples');
%Variables Include:
%SampleWF: Includes WFs from sample Electric feilds
%NumberSavedStates:Number of states saved from each sample
load('Library/Grid');
%Variables Include
%X: contains x cordinates of mesh
%Y: contains y cordinates of mesh
%rw: contains the number of rows of the mesh without boarder of two layers
%col contains the number of columns of the mesh without boarder of two
%layers
%dx:spacing along x-direction
%dy: spacing along y-direction (same as dx)
%ROWGRID: total number of rows of the grid
%COLGRID: total number of columns of the grid
%x: contains x axis grid points
%y: contains y axis grid points
load('Library/TrainingPotentialandMass')
%--------------------Get Total number Of snapshots------------------------%
Nstate=NumberSavedStates;
Nsamp=SampleNumber;
Snapshots=Nsamp*Nstate; %total number of snapshots
%--------------Condense all samples into a two dimension matrix-----------%
WaveF = SampleWFs(:,1:Nstate,1);
for j = 2:Nsamp
    WaveF = [WaveF,SampleWFs(:,1:Nstate,j)];
end
%---------------Normalize all Wave Functions------------------------------%
for i=1:1:size(WaveF,2)
    WaveF(:,i)=normalize(Matrix(WaveF(:,i),ROWGRID,COLGRID),x,y)*WaveF(:,i);
end
%------------Set the polarity of all Wave Functions to be the same--------%
for i=1:1:size(WaveF,2)
   condition=Matrix(WaveF(:,i),ROWGRID,COLGRID);
   if condition(1,1)<0
       WaveF(:,i)=WaveF(:,i)*-1;
   end
end
%------------------------------Create A matrix----------------------------%
Amat_a=zeros(Snapshots,Snapshots);
 for i = 1:Snapshots;
        for j = 1:Snapshots;
            Product=Matrix(WaveF(:,i).*WaveF(:,j),ROWGRID,COLGRID);
            Amat_a(i,j) = trapz(y,trapz(x,Product,2))/Snapshots; 
        end
 end
%----------------------------Solve Eigenvalue Problem---------------------%
[AMatrixVector,UnsortedLamda] =eig(Amat_a);
LamdaUnsortedVector=diag(UnsortedLamda);
[SortedLambda,Index]=sort(LamdaUnsortedVector,'descend');
AmatrixSortedVector=zeros(size(AMatrixVector,1),size(AMatrixVector,2));
for i=1:1:length(SortedLambda)
    AmatrixSortedVector(:,i)=AMatrixVector(:,Index(i));
end
%----------------------Plot Eigen Energy against Mode---------------------%
ModeNumbers=1:1:length(SortedLambda);
NumberofModes=length(ModeNumbers);
Lambda=figure(1);
semilogy(ModeNumbers,SortedLambda,'k-','Linewidth',2.5);
title('Eigenvalues')
xlabel("Mode")
ylabel("Eigenvalue")
ax=gca;
set(ax,'FontSize',12);
%---------------------Create Theoretical Error Plot-----------------------%
errs_lambda_a=zeros(1,NumberofModes-1);
tot_lam_a=sum(abs(SortedLambda));
for m=1:NumberofModes-1
    errs_lambda_a(m)=sqrt(sum(abs(SortedLambda(m+1:NumberofModes)))/tot_lam_a);
end
TheoreticalErrorPlot=figure(2);
semilogy(ModeNumbers(1:end-1), errs_lambda_a*100,'r-o');
title("Theoretical Error");
xlabel("Modes");
ylabel("Errors (%)");
%------------------------------Create Modes-------------------------------%
Modes=zeros(size(WaveF,1),NumberofModes);
for m=1:NumberofModes
   for j=1:1:NumberofModes
       Modes(:,m)=WaveF(:,j)*AmatrixSortedVector(j,m)+Modes(:,m);
   end
   Modes(:,m)=Modes(:,m)/(NumberofModes)/SortedLambda(m);
end
%------------------------------Normalize Each Mode------------------------%
for i=1:1:NumberofModes
   Modes(:,i)=normalize(Matrix(Modes(:,i),ROWGRID,COLGRID),x,y).*Modes(:,i);
end
for i=1:1:NumberofModes
    for j=1:1:NumberofModes
        Orthognaolcheck(i,j)=trapz(y,trapz(x,Matrix(Modes(:,i).*Modes(:,j),ROWGRID,COLGRID),2));
    end
end
%------------------------Create the matrix libary for POD-----------------%
h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
m0 = 9.10938215e-31; % free electron mass in kg
q_el = 1.60218e-19; % elementary charge, in coulombs
% Create the Matrix library for POD
%Create kinetic energy Matrix
T=zeros(LibraryModes);
for i=1:1:LibraryModes
    [FXi,FYi] = gradient(Matrix(Modes(:,i),ROWGRID,COLGRID),dx,dy);
    for j=1:1:LibraryModes
        [FXj,FYj] = gradient(Matrix(Modes(:,j),ROWGRID,COLGRID),dx,dy);
        Integrand=((FXi.*FXj)+(FYi.*FYj)).*h_bar.^2/(2*m0)./(EffectiveMassMatrix);
        T(i,j)=trapz(y,trapz(x,Integrand,2))/(q_el);
    end
end

%One does not know the electric field which will be applied during
%application. However,we already know the underlying potential due to the
%structure. The total potential is UTotal=Ustruct+Ex*x+Ey*y. Hence, we can
%precompute the integrals. The integral of  mode_i*mode_j*x and
%mode_i*mode_j*y will be multiplied by Ex and Ey respectively. 
UStruct=zeros(LibraryModes); % 
for i=1:1:LibraryModes;
  for j=1:1:LibraryModes;
      IntegralProduct=(Matrix(Modes(:,i),ROWGRID,COLGRID).*Matrix(Modes(:,j),ROWGRID,COLGRID)).*U;
      UStruct(i,j)=trapz(y,trapz(x,IntegralProduct,2));
  end    
end
% Compute matrix of mode_i*mode_j*x
UX=zeros(LibraryModes); % 
for i=1:1:LibraryModes;
  for j=1:1:LibraryModes;
      IntegralProduct=(Matrix(Modes(:,i),ROWGRID,COLGRID).*Matrix(Modes(:,j),ROWGRID,COLGRID)).*X;
      UX(i,j)=trapz(y,trapz(x,IntegralProduct,2));
  end    
end
% Compute matrix of mode_i*mode_j*y
UY=zeros(LibraryModes); % 
for i=1:1:LibraryModes;
  for j=1:1:LibraryModes;
      IntegralProduct=(Matrix(Modes(:,i),ROWGRID,COLGRID).*Matrix(Modes(:,j),ROWGRID,COLGRID)).*Y;
      UY(i,j)=trapz(y,trapz(x,IntegralProduct,2));
  end    
end
%------------------------------Save data----------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the POD Matrix Library
save('Library/MatrixLibrary','T','UStruct','UX',"UY")
%Save the Modes and the number of total modes
save('Library/PODMODES','Modes','NumberofModes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save Theoretical Error plot
save('Library/TheoreticalError','errs_lambda_a');
%-----------------------------Save Figures--------------------------------%
%Save the Eigenvalue Energy vs modes graph
savefig(Lambda,'Library/Lambda');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the EigenvalueError vs modes graph
savefig(TheoreticalErrorPlot,'Library/TheoreticalErrorPlot');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




