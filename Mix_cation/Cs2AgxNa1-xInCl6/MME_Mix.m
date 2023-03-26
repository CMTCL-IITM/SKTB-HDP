% Slater-Koster tight-binding model for mixed cation halide double perovskite (Cs2AgxNa1-xInCl6).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
% Read the comments carefully before run the code.
% Change the setup in "full_Hamiltonian_pd_sc.m" for each x value.
%% clear the workspace
close all
clear
clc
tic
%% Tight-binding parameters
selectcompd = 2; % Input
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
    
    EsB	 =   4.24;
    EpB	 =    7.98;
    VssB = 	-0.0;
    VspB =	0.018;
    VppsB =	-0.092;
    VpppB =	0.00;
    
    VsAsB =	0.375;
    VspAB =	0.48;
    VsdAB =	0.45;
    Vpds =	0.63;
    Vpdp =	0.24;
    
    sA	=   0.09;
    sB	=   0.16;
    
    %%  CsNaInCl
    
    EsNa	=   12.4;
    EpNa	=   13.5;
    VssNa   =  -0.02;
    VspNa   =	0.018;
    VppsNa  =  -0.282;
    VpppNa  =	0.05;
    
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
    
    EsA	 = 6.6;
    Ez2	 = -0.06;
    Et2g =	-0.84;
    Vdds =	-0.01;
    Vddp =	0.007;
    Vddd =	0;
    VssA =	0;
    VsdA =	0;
    
    EsB	 =  4.08;
    EpB	 =  7.76;
    VssB =	-0.012;
    VspB =	 0.018;
    VppsB =	-0.092;
    VpppB =	 0;
    
    VsAsB =	0.375;
    VspAB =	0.48;
    VsdAB =	0.45;
    Vpds  =	0.63;
    Vpdp  =	0.24;
    
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
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
if selectcompd == 2
    nkpt = 199; % for Ag0.5Na0.5: nkpt = 199
else
    nkpt = 203;
end

nbnd =  999;
% crystal structure matrix
a = 10.6;     % lattice constant
rvec = [a, 0, 0;
        0, a, 0;
        0, 0, a];
% calculate reciprocal lattice vectors
Vol = dot(rvec(1,:),cross(rvec(2,:),rvec(3,:)));
b1 =  2*pi*cross(rvec(2,:),rvec(3,:))/Vol;
b2 =  2*pi*cross(rvec(3,:),rvec(1,:))/Vol;
b3 =  2*pi*cross(rvec(1,:),rvec(2,:))/Vol;

recip = [b1;b2;b3];
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
t1 = [Ez2,Et2g,Vdds,Vddp,Vddd,VssA,EsA,VsdA];  % A-A
t2 = [VssB,VspB,VppsB,VpppB,EpB,EsB];          % B-B
t3 = [VsAsB,VspAB,VsdAB,Vpds,Vpdp];            % A-B
tNa1 = [VssNa,VspNa,VppsNa,VpppNa,EpNa,EsNa];
tNaIn = [VssBNa,VspBNa,VppsBNa,VpppBNa,0,0];
tNaAg = [VsAsNa,VspANa,VsdANa,VpdsANa,VpdpANa];

%%
momentum = [];
ndiv = 10;
ss = 1;
%%
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
if selectcompd == 2
    kpath=load('klist_sc_0.5.dat'); % change for Ag0.5Na0.5
else
    kpath=load('klist_sc.dat');
end
% diagonolize the Hamiltonian

for istate = 19:20       % VB
    for fstate = 21:22   % CB
        
        for i = 1:length(kpath)
            
            kx = kpath(i,1)*2*pi;
            ky = kpath(i,2)*2*pi;
            kz = kpath(i,3)*2*pi;
            
            H0 = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd);
            % store the eigenvalues
            [v,e] = eig(H0);
            [energy1,perm] = sort(real(diag(e)'));
            energy(i,:) = energy1;
            vectors(:,:) = v(:,perm);
            % Calculation of MME
            %%
            dk = 2*pi/(ndiv-1);
            
            dHx = (full_Hamiltonian_pd_sc(kx+dk,ky,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd)-full_Hamiltonian_pd_sc(kx-dk,ky,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd))/(2*dk);
            dHy = (full_Hamiltonian_pd_sc(kx,ky+dk,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd)-full_Hamiltonian_pd_sc(kx,ky-dk,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd))/(2*dk);
            dHz = (full_Hamiltonian_pd_sc(kx,ky,kz+dk,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd)-full_Hamiltonian_pd_sc(kx,ky,kz-dk,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd))/(2*dk);
            
            px = transpose(conj(v(:,fstate)))*dHx*v(:,istate) ;
            py = transpose(conj(v(:,fstate)))*dHy*v(:,istate) ;
            pz = transpose(conj(v(:,fstate)))*dHz*v(:,istate) ;
            
            momentum(i,ss,:) = [kx,ky,kz,px,py,pz,px*conj(px)+py*conj(py)+pz*conj(pz)];
            
        end
        ss = ss+1;
    end
end
% constants
m0 = 9.11e-31;
hbar = 1.054e-34;
prefact = m0/hbar;
%% k-length calculation
kpts = [];
for j = 1: length(kpath)
    kpts = [kpts;sum(recip.*transpose(kpath(j,1:3)))];
end
s = zeros(1,size(kpts,1));
for i = 1:size(kpts,1)-1
    d = sqrt(sum((kpts(i,:)-kpts(i+1,:)).^2));
    s(i+1) = s(i) + d;
end
% plot the MME and TB band structure
figure(1)
plot(s,energy,'r','LineWidth',2);
figure(2)
plot(s,prefact*sum(momentum(:,:,end),2)/(nkpt*ss),'k','LineWidth',2)
%% axis properties
ax=gca;
xlim([0,s(end)]);
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 22;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];
ax.YLim = [0,2];

toc
