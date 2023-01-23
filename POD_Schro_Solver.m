%%%%%%%%%%%%%%Generate Wave Functions Using Global 1 Model%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code To generate wave functions and plots using POD model
%----------------Clear Command Window, Workspace and close Figures--------%
clear;
clc;
close all;
%Solve DNS directly for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Library/PODMODES');
%Variables Include:
%Modes: Contains POD modes
%NumberofModes: Contains Total Number of Modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Library/MatrixLibrary')
%Variables Include: 
%T Kinetic energy matrix
%UStruct, UX and UY Components of potential energy matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Library/TheoreticalError')
%Variables Include:
%errs_lambda_a: Contains Theoretical Error
%-------------------------------Run Input 2-------------------------------% 
POD_Solver_Parameters;
DNS_Compare; % Run DNSCOMPARE to create DNS data to compare against POD solution
%Variables Include: 
%TotalModesInHamiltonian:Modes to Construct POD Hamiltionian
%UsedModes: Number of Modes Used for Plotting WF
%NumberofModesToPlot:Number of Modes Used for LS plot
%-----------------------------Constants-----------------------------------%
PotentialUsed=SavePotential(:,:,end);
[Shifting,Index]=min(PotentialUsed,[],'all','linear');
h_bar = 1.054571726e-34; % Reduced Planck Constant in J*s
m0 = 9.10938215e-31; % free electron mass in kg
q_el = 1.60218e-19; % elementary charge, in coulombs
%------------------------------RenameDNS Data-----------------------------% 
DNS=SavedWFDNSTest; 
EEDNS=SavedEnergiesDNSTest;
NumberofStates=StatesToBeCollectedForPODComparision;
%------------------------CREATE Hamiltonian MATRIX------------------------% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create and solve Hamiltonian Matrix
H=UStruct(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian)+...
  T(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian)+...
  XField*100000*UX(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian)+...
  YField*100000*UY(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian);
[a,EE]=eig(H);
[EE,INDEX]=sort(diag(EE));
a=a(:,INDEX);
%--------------------------Plot each WF-----------------------------------%
Figures=gobjects(NumberofStates*2,1);
Linestyle= {'r-','b:','g--','M-.'};
Linesize=[2,4,2.25,2];
for k=1:NumberofStates;
    %initial conditions for each loop:
    %Set iteration number to 1
    %Grab the desired number of used modes for this particular state
    StateModes=UsedModes(k,:);
    StateModes(StateModes==0)=[];
    iterationnumber=1;
    %Normalize the DNS wavefunction created from DNS code 
    DNS(:,k)=normalize(Matrix(DNS(:,k),ROWGRID,COLGRID),x/(10^(-9)),y/(10^(-9)))*DNS(:,k);
   %---------------------------Select Figure and create Title-------------%
    Figures(k)=figure(k);
    set(gcf,'color','w');   
    TitleString=strcat(" |\Psi|^2 of State ",string(k), " in x and y Direction");
    xtitle=strcat(" |\Psi|^2 of State ",string(k), " in x Direction");
    ytitle=strcat(" |\Psi|^2 of State ",string(k), " in y Direction");
    %-------------------------Determine Scaling Factor--------------------%
    [Max,Index]=max(Matrix(DNS(:,k).^2,ROWGRID,COLGRID),[],'all','linear');
    [RowOfIntrest,ColOfIntrest]=ind2sub([ROWGRID,COLGRID],Index);
    Wavefactor=Depth/(2*Max);
    %Scale DNS acording to scaling factor
    WaveStateDNS=Matrix(Wavefactor*DNS(:,k).^2+EEDNS(k)-Shifting,ROWGRID,COLGRID);
    %-----------------------Plot Potential and DNS for x Subplot----------%
    subplot(1,2,1);
    hold on;
    plot(x/(10^(-9)),PotentialUsed(RowOfIntrest,:)-Shifting,'k','LineWidth',1.5,'HandleVisibility','off');
    plot(x/(10^(-9)),WaveStateDNS(RowOfIntrest,:),'k-','LineWidth',2.75,'DisplayName','DNS');
    title(xtitle)
    xlabel('X Position (nm)','FontSize',12);
    ylabel('Energy Band (eV)','FontSize',12);
    box on
    %----------------------Plot Potential and DNS for y Subplot-----------%
    subplot(1,2,2);
    hold on;
    %plot potential and DNS
    plot(y/(10^(-9)),PotentialUsed(:,ColOfIntrest)-Shifting,'k','LineWidth',1.5,'HandleVisibility','off');
    plot(y/(10^(-9)),WaveStateDNS(:,ColOfIntrest),'k-','LineWidth',2.75,'DisplayName',"DNS");
    title(ytitle)
    xlabel('Y Position (nm)','FontSize',12);
    ylabel('Energy Band (eV)','FontSize',12);
    box on
    sgtitle(TitleString);
    %----------------------------Add POD plot to every figure-------------%
        for m=StateModes;
            
            %-------------------------Create POD WF and Normalize---------%
            WavePOD=Modes(:,1:m)*a(1:m,k);
           
            WavePOD=normalize(Matrix(WavePOD,ROWGRID,COLGRID),x/(10^(-9)),y/(10^(-9)))*WavePOD;
            WaveStatePOD=Matrix(Wavefactor*WavePOD.^2+EEDNS(k)-Shifting,ROWGRID,COLGRID);
            %--------------------Add POD Plot to x subplot----------------%
            subplot(1,2,1);
            hold on
            %plot POD wave function
            plot(x/(10^(-9)),WaveStatePOD(RowOfIntrest,:),Linestyle{iterationnumber},'LineWidth',Linesize(iterationnumber),'DisplayName',strcat("POD:",string(m)));
           %-------------------Add legend During Last Iteration-----------%
            if iterationnumber==length(StateModes)
                axis([x(1)/(10^(-9)),x(end)/(10^(-9)),0,Depth/2+HeightAboveWave]);
            Leg=legend('location','southeast');
            Leg.FontSize=10
            legend boxon;
            ax=gca
            pbaspect(ax,[0.7598257575,1,0.7598257575])
            end
            %----------------------Add POD Plot to y subplot--------------%
            subplot(1,2,2);
            hold on;
            plot(y/(10^(-9)),WaveStatePOD(:,ColOfIntrest),Linestyle{iterationnumber},'LineWidth',Linesize(iterationnumber),'MarkerSize',7,'DisplayName',strcat("POD:",string(m)));
            %------------------Add Legend During Last Iteration-----------%
            if iterationnumber==length(StateModes)
            axis([y(1)/10^(-9),y(end)/(10^(-9)),0,Depth/2+HeightAboveWave])
            Leg=legend('location','southeast');
            Leg.FontSize=10
            legend boxon;
            ax=gca
            pbaspect(ax,[0.7598257575,1,0.7598257575])
            end
            %Increase Iteration Number By 1
            iterationnumber=iterationnumber+1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create Contours%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Figures(NumberofStates+k)=figure(NumberofStates+k);
    set(gcf,'color','w');
    WaveStatePOD=Matrix(WavePOD.^2+EEDNS(k)-Shifting,ROWGRID,COLGRID);
    WaveStateDNS=Matrix(DNS(:,k).^2+EEDNS(k)-Shifting,ROWGRID,COLGRID);
    %-------------------------Contour of POD------------------------------%
    S1=subplot(2,1,1);
    hold on
    contour(X./(10^(-9)),Y./(10^(-9)),WaveStatePOD,'Fill','on')
    %Create Lines to show where WF plots where taken from
    plot(x(ColOfIntrest)*ones(length(y),1)/(10^(-9)),y/(10^(-9)),'r--','LineWidth',1.5);
    plot(x/(10^(-9)),y(RowOfIntrest)*ones(length(x),1)/(10^(-9)),'r--','LineWidth',1.5);
    % Set the axis to square
    set(gca,'Layer','top');
    set(gca,'Box','on');
    axis('equal')
    title("POD: State "+string(k)+" using "+string(StateModes(end))+" Modes")
    %-------------------Contour of DNS------------------------------------%
    S2=subplot(2,1,2);
    contour(X/(10^-9),Y/(10^(-9)),WaveStateDNS,'Fill','on');
    title("DNS: State "+string(k))
    hold on;
    %Plot axis to know where to look
    plot(x(ColOfIntrest)*ones(length(y),1)/(10^(-9)),y/(10^(-9)),'r--','LineWidth',1.5);
    plot(x/(10^(-9)),y(RowOfIntrest)*ones(length(x),1)./(10^(-9)),'r--','LineWidth',1.5);
    %set axis to square
    set(gca,'Layer','top');
    set(gca,'Box','on');
    axis('equal')
  
end
%-------------------------------Percent Difference in Energy--------------%
PercentDifference=abs(EE(1:NumberofStates)-Shifting-(EEDNS(1:NumberofStates)-Shifting))./((EE(1:NumberofStates)-Shifting+EEDNS(1:NumberofStates)-Shifting)./2).*100
EnergyCompare=[EE(1:NumberofStates)-Shifting,EEDNS(1:NumberofStates)-Shifting,PercentDifference];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Create LS Error Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNS=SavedWFDNSTest; %Get data from DNS
WaveFunctionerror=zeros(NumberofStates,NumberofModesToPlot);%Stores Error per mode
AverageError=zeros(NumberofModesToPlot,1);
Numerators=zeros(NumberofStates,NumberofModesToPlot);%Stores Numerator for finding LS
%---------------------------------Create Hamiltonian Again----------------%
tic % Start Timer
H=UStruct(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian)+...
  T(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian)+...
  XField*100000*UX(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian)+...
  YField*100000*UY(1:TotalModesInHamiltonian,1:TotalModesInHamiltonian);
[a,EE]=eig(H);
[EE,INDEX]=sort(diag(EE));
a=a(:,INDEX);
HamTime=toc; %End Timer

%-------------------For each # of Modes create POD WF and find LS error---%
for m=1:1:NumberofModesToPlot
    tic
    WavesPOD=Modes(:,1:m)*a(1:m,1:NumberofStates);
    CreationTime=toc; % determine time to create the WF for each state
    TotalTime(m)=HamTime+CreationTime; % calculate total time
    for k=1:1:NumberofStates
    WavesPOD(:,k)=normalize(Matrix(WavesPOD(:,k),ROWGRID,COLGRID),x,y)*WavesPOD(:,k);
    end
    for k=1:1:NumberofStates
        POSPODMATRIX=Matrix(WavesPOD(:,k),ROWGRID,COLGRID);
        NEGPODMATRIX=Matrix(-1*WavesPOD(:,k),ROWGRID,COLGRID);
        DNSMATRIX=Matrix(DNS(:,k),ROWGRID,COLGRID);
        IntegrandPOS=(POSPODMATRIX-DNSMATRIX).^2;
        IntegrandNeg=(NEGPODMATRIX-DNSMATRIX).^2;
        POSSUM=trapz(y,trapz(x,IntegrandPOS,2));
        NEGSUM=trapz(y,trapz(x,IntegrandNeg,2));
        if POSSUM>NEGSUM
            NUMERATOR=NEGSUM;
        else
            NUMERATOR=POSSUM;
        end
        Numerators(k,m)=NUMERATOR;
        DENOM=trapz(y,trapz(x,DNSMATRIX.^2,2));
        WaveFunctionerror(k,m)=sqrt(NUMERATOR)/sqrt(DENOM);
    end
end
%---------------------------Calculate the average Error-------------------%
sumErrorsquared=sum(Numerators(1:NumberSavedStates,:),1);
for i=1:length(AverageError)
   AverageError(i)=sqrt((sumErrorsquared(i)/NumberSavedStates));
end
%------------------------------Create LS Error Plot-----------------------%
LSEPlot=figure(2*NumberofStates+1)
hold on
%Plot LS error for each state
for k=1:NumberofStates   
semilogy(1:1:NumberofModesToPlot,WaveFunctionerror(k,:)*100,'-s','LineWidth',1.5);
 legendInfo{k} = strcat("State: ",string(k));
end

semilogy(1:1:NumberofModesToPlot,AverageError*100,'-- m','LineWidth',1.5);
%Plot theoretical lambda
semilogy(1:1:NumberofModesToPlot,errs_lambda_a(1:1:NumberofModesToPlot)*100,'-s k','LineWidth',1.5)
legendInfo{NumberofStates+1} = "Average Error";
legendInfo{NumberofStates+2}="Theoretical";
set(gcf,'color','w')
ylabel("LS Error (%)")
xlabel("Number of Modes")
ax = gca;
ax.YAxis.Scale = 'log';
axis([-inf,inf,4*10^(-2),200])
set(ax,'FontSize',12);
legend(legendInfo);
%--------------------------Print Energies To Command Window---------------%
format long
disp(EnergyCompare)
%-------------------------------Create Table to Save Energies-------------%
EnergyTable=table(EnergyCompare(:,1),EnergyCompare(:,2),EnergyCompare(:,3));
EnergyTable.Properties.VariableNames([1 2 3]) = {'POD_Energy' 'DSN_Energy','Percent_Difference'}
Lambda=openfig('Library/Lambda');
ax=gca;
TheoreticalErrorPlot=openfig('Library/TheoreticalErrorPlot');
%-----------------------------Delete and remake Results Directory---------%
%Delete old Results Folder
Status=rmdir('Results','s');
%Make results folder
mkdir Results
%----------------------------Save Data to Results Folder------------------%
%Save Energy Table in Results folder 
writetable(EnergyTable,'Results/EnergyTable.csv')
%Save wavefunction figures(side view and contour plots in results folder
for i=1:1:NumberofStates
   saveas(Figures(i),strcat("Results/ProfileState",string(i),".fig"))
end
% Write a file outputting the times
fileID = fopen('./Results/TimeFile.txt','w');
fprintf(fileID,'This file records computation time results for POD and DNS.\n');
fprintf(fileID,'DNS time: %fs \n',DNSEndTime);
fprintf(fileID,'POD Time Results: \n');
fprintf(fileID,'Hamiltonian time creation (%i) Modes : %fs \n',TotalModesInHamiltonian,HamTime);
fprintf(fileID,"Total time for solution composed of M modes: \n");
for i=1:NumberofModesToPlot
    Improvement=DNSEndTime/TotalTime(i);
    fprintf(fileID," POD computation time and improvement for  %i Modes: %fs, %f times quicker \n",i,TotalTime(i),Improvement);
end
fclose(fileID);
%Save contours to results
for i=NumberofStates+1:1:2*NumberofStates
   saveas(Figures(i),strcat("Results/ContourState",string(i-NumberofStates),".fig"))
end
saveas(LSEPlot,"Results/LSE.fig")
saveas(Lambda,"Results/ModeEigenvalue.fig")
saveas(TheoreticalErrorPlot,"Results/TheoreticalErrorPlot.fig")

