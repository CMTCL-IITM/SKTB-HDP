% Slater-Koster tight-binding model for mixed cation halide double perovskite (Cs2AgInxBi1-xCl6).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace
clear
clc
close all
tic
%% Tight-binding parameters
% Select the compound
selectcompd = 2;  % 1 for In0.25Bi0.75, 2 for In0.5Bi0.5, and 3 for In0.75Bi0.25
%%  AgInCl	AgBiCl for In0.25Bi0.75
if selectcompd == 1
    EsA =	[6.36	6.36];
    Eeg	=   [-0.65	-0.65];
    Et2g=	[-1.5	-1.5];
    Vdds=	[-0.014	-0.014];
    Vddp=	[0.007	0.007];
    Vddd=	[0.0	0.0];
    VssA=	[0.0	0.0];
    VsdA=	[0.0	0.0];
    
    EsB	=   [3.27	-1.35];
    EpB	=   [7.36	3.75];
    VssB= 	[0.0	-0.03];
    VspB =	[0.018	0.03];
    VppsB =	[-0.092	-0.054];
    VpppB =	[0.0	0.0];
    
    VsAsB=	[0.375	0.27];
    VspAB=	[0.48	0.565];
    VsdAB=	[0.45	0.223];
    Vpds=	[0.63	0.6];
    Vpdp=	[0.24	0.0];
    
    sA	=   [0.09	0.09];
    sB	=   [0.16	0.54];
end
%%	AgInCl	AgBiCl for In0.5Bi0.5
if selectcompd == 2
    EsA =	[6.39	6.39];
    Eeg	=   [-0.761	-0.761];
    Et2g=	[-1.58	-1.58];
    Vdds=	[-0.014	-0.014];
    Vddp=	[0.007	0.007];
    Vddd=	[0.0	0.0];
    VssA=	[0.0	0.0];
    VsdA=	[0.0	0.0];
    
    EsB	=   [3.21	-1.14];
    EpB	=   [7.10	3.70];
    VssB= 	[0.0	-0.03];
    VspB =	[0.018	0.03];
    VppsB =	[-0.092	-0.054];
    VpppB =	[0.0	0.0];
    
    VsAsB=	[0.375	0.27];
    VspAB=	[0.48	0.565];
    VsdAB=	[0.45	0.223];
    Vpds=	[0.63	0.6];
    Vpdp=	[0.24	0.0];
    
    sA	=   [0.09	0.09];
    sB	=   [0.16	0.54];
end
%%  	AgInCl	AgBiCl for In0.75Bi0.25
if selectcompd == 3
    EsA =	[6.75	6.75];
    Eeg	=   [-0.361	-0.361];
    Et2g =	[-1.14	-1.14];
    Vdds =	[-0.014	-0.014];
    Vddp =	[0.007	0.007];
    Vddd =	[0.0	0.0];
    VssA =	[0.0	0.0];
    VsdA =	[0.0	0.0];
    
    EsB	=   [3.65	-0.66];
    EpB	=   [7.32	4.08];
    VssB = 	[0.0	-0.03];
    VspB =	[0.018	0.03];
    VppsB =	[-0.092	-0.054];
    VpppB =	[0.0	0.0];
    
    VsAsB =	[0.375	0.27];
    VspAB =	[0.48	0.565];
    VsdAB =	[0.45	0.223];
    Vpds =	[0.63	0.6];
    Vpdp =	[0.24	0.0];
    
    sA	=   [0.09	0.09];
    sB	=   [0.16	0.54];
end
%%
if selectcompd == 1
    dft = load('AgInBi_0.25.dat');
elseif selectcompd == 2
    dft = load('AgInBi_0.5.dat');
elseif selectcompd == 3
    dft = load('AgInBi_0.75.dat');
end
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
kpath=load('klist_sc.dat');
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt = 203;
nbnd =  999;
% position coordinates of 2nd-nearest and 4th-nearest neighbor B and B' atoms in crystal unit
% 4th-nearest neighbor
posAA = [1  1  0;
        -1 -1  0;
        -1  1  0;
        1 -1  0;
        0  1  1;
        0 -1 -1;
        0  1 -1;
        0 -1  1;
        1  0  1;
        -1  0 -1
        -1  0  1;
        1  0 -1];

posAB = [1  0  0;
        -1  0  0;
        0 -1  0;
        0  1  0;
        0  0  1;
        0  0 -1];
% Desigh the Hamiltonian and calculate the eigenvalues along the high-symmetry kpath.
% hopping parameters
for p = 1:2
    t1(p,:) = [Eeg(p),Et2g(p),Vdds(p),Vddp(p),Vddd(p),VssA(p),EsA(p),VsdA(p)];   % A-A
    t2(p,:) = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];                 % B-B
    t3(p,:) = [VsAsB(p),VspAB(p),VsdAB(p),Vpds(p),Vpdp(p)];                      % A-B
end
% diagonolize the Hamiltonian
for i= 1:length(kpath)
    kx = kpath(i,1)*2*pi;
    ky = kpath(i,2)*2*pi;
    kz = kpath(i,3)*2*pi;
    % get the Hamiltonian matrix
    H = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd);
    % get the eigenvalues
    v = eig(H);
    % store the eigenvalues
    eig_val(i,:) = v';
end
%% plot the DFT and TB band structure
plot(dft(1:nkpt,1),reshape(dft(:,2),[nkpt,nbnd]),'r','LineWidth',2)
hold on
plot(dft(1:nkpt,1),eig_val,'b','LineWidth',2);
hold off
%% axis properties
xlim([0,dft(nkpt,1)]);
ylim([-1,5]);
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 0.001;
ax.TickDir = 'in';
ax.TickLength = [0.001 0.001];
ax.FontSize = 22;

toc
