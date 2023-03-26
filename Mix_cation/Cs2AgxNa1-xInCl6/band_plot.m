% Slater-Koster tight-binding model for mixed cation halide double perovskite (Cs2AgxNa1-xInCl6).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
% Read the comments carefully before run the code.
% Change the setup in "full_Hamiltonian_pd_sc.m" for each x value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
tic
%% Tight-binding parameters
selectcompd = 2;
%% CsAg0.75Na0.25InCl
if selectcompd == 1
    % CsAgInCl
    EsA =	6.93;
    Ez2	=   -0.06;
    Et2g=	-0.8;
    Vdds=	-0.01;
    Vddp=	0.007;
    Vddd=	0.0;
    VssA=	0.0;
    VsdA=	0.0;
    
    EsB	=   4.24;
    EpB	=    7.98;
    VssB= 	-0.0;
    VspB =	0.018;
    VppsB =	-0.092;
    VpppB =	0.00;
    
    VsAsB=	0.375;
    VspAB=	0.48;
    VsdAB=	0.45;
    Vpds=	0.63;
    Vpdp=	0.24;
    
    sA	=   0.09;
    sB	=   0.16;
    
    %%  CsNaInCl
    
    EsNa	=   12.4;
    EpNa	=   13.5;
    VssNa   =  -0.02;
    VspNa   =	0.018;
    VppsNa  =  -0.282;
    VpppNa =	0.05;
    
    VsAsNa  = -0.0;
    VspANa  =  0.0;
    VsdANa  =  0.04;
    VpdsANa =  0.16;
    VpdpANa =  0.16;
    
    VssBNa  = -0.24;
    VspBNa  =  0.16;
    VppsBNa =  -0.;
    VpppBNa =   0.;
end
%%  CsAg0.5Na0.5InCl
if selectcompd == 2
    EsA	= 6.6;
    Ez2	 = -0.06;
    Et2g =	-0.84;
    Vdds =	-0.01;
    Vddp =	0.007;
    Vddd =	0;
    VssA =	0;
    VsdA=	0;
    
    EsB	= 4.08;
    EpB	= 7.76;
    VssB =	-0.012;
    VspB =	0.018;
    VppsB =	-0.092;
    VpppB =	0;
    
    VsAsB=	0.375;
    VspAB=	0.48;
    VsdAB=	0.45;
    Vpds=	0.63;
    Vpdp=	0.24;
    
    sA	= 0.09;
    sB	= 0.16;
    
    EsNa	=	12.4;
    EpNa	=	13.5;
    VssNa =		-0.02;
    VspNa =		0.018;
    VppsNa =	-0.282;
    VpppNa =	0.05;
    
    VsAsNa  = -0.0;
    VspANa  =  0.0;
    VsdANa  =  0.04;
    VpdsANa =  0.16;
    VpdpANa =  0.16;
    
    
    VssBNa	=	-0.24;
    VspBNa	=	0.16;
    VppsBNa	=	0;
    VpppBNa	=	0;
end
%%  CsAg0.25Na0.75InCl
if selectcompd == 3
    EsA	= 7.23;
    Ez2	 = 0.195;
    Et2g =	-0.56;
    Vdds =	-0.01;
    Vddp =	0.007;
    Vddd =	0;
    VssA =	0;
    VsdA=	0;
    
    EsB	= 4.44;
    EpB	= 8;
    VssB =	-0.012;
    VspB =	0.018;
    VppsB =	-0.092;
    VpppB =	0;
    
    VsAsB=	0.375;
    VspAB=	0.48;
    VsdAB=	0.45;
    Vpds=	0.63;
    Vpdp=	0.24;
    
    sA	= 0.09;
    sB	= 0.16;
    
    EsNa	=	12.4;
    EpNa	=	13.5;
    VssNa =		-0.02;
    VspNa =		0.018;
    VppsNa =   -0.282;
    VpppNa =	0.05;
    
    VsAsNa  = -0.0;
    VspANa  =  0.0;
    VsdANa  =  0.04;
    VpdsANa =  0.16;
    VpdpANa =  0.16;
    
    
    VssBNa	=		-0.24;
    VspBNa	=		0.16;
    VppsBNa	=		0;
    VpppBNa	=		0;
end
%%
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
if selectcompd == 2
    kpath=load('klist_sc_0.5.dat'); % change for Ag0.5Na0.5
else
    kpath=load('klist_sc.dat');
end
% Select the compound
if selectcompd == 1
    dft = load('CsAg0.75Na0.25InCl.dat');
elseif selectcompd == 2
    dft = load('CsAg0.5Na0.5InCl.dat');
elseif selectcompd == 3
    dft = load('CsAg0.25Na0.75InCl.dat');
end
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
if selectcompd == 2
    nkpt = 199; % for Ag0.5Na0.5: nkpt = 199
else
    nkpt = 203;
end

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

% Design the Hamiltonian and calculate the eigenvalues along the high-symmetry kpath.
% hopping parameters
t1 = [Ez2,Et2g,Vdds,Vddp,Vddd,VssA,EsA,VsdA];   % A-A
t2 = [VssB,VspB,VppsB,VpppB,EpB,EsB];          % B-B
t3 = [VsAsB,VspAB,VsdAB,Vpds,Vpdp];                 % A-B
tNa1 = [VssNa,VspNa,VppsNa,VpppNa,EpNa,EsNa];
tNaIn = [VssBNa,VspBNa,VppsBNa,VpppBNa,0,0];
tNaAg = [VsAsNa,VspANa,VsdANa,VpdsANa,VpdpANa];
% diagonolize the Hamiltonian

for i= 1:length(kpath)
    kx = kpath(i,3)*2*pi;
    ky = kpath(i,2)*2*pi;
    kz = kpath(i,1)*2*pi;
    % get the Hamiltonian matrix
    H = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd);
    % get the eigenvalues
    v = eig(H);
    % store the eigenvalues
    eig_val(i,:) = v';
end

%% plot the DFT and TB band structure
plot(dft(1:nkpt,1),reshape(dft(:,2),[nkpt,nbnd]),'r','LineWidth',1.5)
hold on
plot(dft(1:nkpt,1),eig_val,'b','LineWidth',1.5);
hold off
%% axis properties
xlim([0,dft(nkpt,1)]);
ylim([-1,8]);
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 0.001;
ax.TickDir = 'in';
ax.TickLength = [0.001 0.001];
ax.FontSize = 22;