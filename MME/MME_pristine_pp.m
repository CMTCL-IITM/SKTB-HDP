% Momentum matrix element (MME) calculation for halide double perovskite (Cs2BB'X6,
% B = In/Tl, B' = Bi/Sb, and X = Cl).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace
clear
clc
close all
addpath(genpath([fileparts(fileparts(pwd)), filesep]));
%% Tight-binding parameters
%%     CsInBiCl   CsInSbCl   CsTlBiBr   CsTlSbBr  CsTlBiCl   CsTlSbCl
EsA =	[-0.98	-1.16	-1.1	-1.08  -1.02	-1.16];
EpA =	[4.50	4.07	5.24	4.60    6.10	5.6];
VssA = 	[-0.053	-0.008	-0.022	0.00   -0.022	0];
VspA = 	[-0.016	0.036	0.044	0.028   0.044	0.028];
VppsA = [0	    0    	0.04	0.0     0.04	0.0];
VpppA = [0	    0	    0.0	    0.0     0.0	    0.0];

EsB =	[-3.156	-2.676	-1.22	-1.62 -1.54	-1.68];
EpB =	[2.42	1.9	     3.25	 2.8   3.75	 3.3];
VssB = 	[0	    0	     0.005	0      0.00   0];
VspB = 	[0	    0	     0	    0      0	   0];
VppsB = [0	    0	     0	    0      0	   0];
VpppB = [0	    0	     0	    0      0	   0];

VssAB =	[-0.377	-0.303	-0.21	-0.22  -0.235	-0.235];
VspAB =	[0.44	0.51	0.50	0.470  0.470	0.46];
VppsAB =[0.665	0.665	0.675	0.720  0.570	0.65];
VpppAB =[0.09	0.09	0.0899	0.0785 0.0899	0.0785];
% SOC strength
sA =	[0.16	0.16	0.46	0.37 0.46	0.37];
sB =	[0.56	0.17	0.49	0.17 0.49	0.17];
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt =  202;
nbnd =  [322 314  414  406 324 324];
% Select the compound
p = 5;
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
kpath=load('../DFT_data/klist.dat');
% hopping parameters
t1 = [VssA(p),VspA(p),VppsA(p),VpppA(p),EpA(p),EsA(p)];
t2 = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];
t3 = [VssAB(p),VspAB(p),VppsAB(p),VpppAB(p),0,0];
% size of basis set
nbasis = 16;
%% crystal structure matrix
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
posAA = [1 1 0;
        -1 1 0;
        -1 -1 0;
         1 -1 0;
         0 1 1;
         0 -1 -1;
         0 1 -1;
         0 -1 1;
         1 0 1;
        -1 0 1;
        -1 0 -1
         1 0 -1];

posAB = [1 0 0;
        -1 0 0;
         0 -1 0;
         0  1 0;
         0 0 1;
         0 0 -1];
% step size in k-vector
ndiv = 10;
%% Desigh the Hamiltonian and calculate the eigenvalues along the high-symmetry kpath.
momentum = [];
ss = 1;
% initial and final states between which the the transition is happening.
% (check the eig_val data for VBE and CBE).
for istate = 3:4       % VB
    for fstate = 5:6   % CB
        % diagonolize the Hamiltonian
        for i = 1:length(kpath)
            
            kx = kpath(i,1)*2*pi;
            ky = kpath(i,2)*2*pi;
            kz = kpath(i,3)*2*pi;
            % get the Hamiltonian matrix
            H0 = full_Hamiltonian_pp(kx,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p));
            % get the eigenvalues
            [v,e] = eig(H0);
            % store the eigenvalues
            eig_val(i,:) = diag(e)';
            % Calculation of MME
            dk = 2*pi/(ndiv-1);
            dHx = (full_Hamiltonian_pp(kx+dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pp(kx-dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
            dHy = (full_Hamiltonian_pp(kx,ky+dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pp(kx,ky-dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
            dHz = (full_Hamiltonian_pp(kx,ky,kz+dk,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pp(kx,ky,kz-dk,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
            
            px = transpose(conj(v(:,fstate)))*dHx*v(:,istate) ;
            py = transpose(conj(v(:,fstate)))*dHy*v(:,istate) ;
            pz = transpose(conj(v(:,fstate)))*dHz*v(:,istate) ;
            
            momentum(i,ss,:) = [kx,ky,kz,px,py,pz,px*conj(px)+py*conj(py)+pz*conj(pz)];
        end
        ss = ss + 1;
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
plot(s,eig_val,'r','LineWidth',2);
figure(2)
plot(s,prefact*sum(momentum(:,:,end),2)/(nkpt*ss),'k','LineWidth',2)
%% axis properties
ylim([0,6]);
xlim([0 s(end)]);
ax=gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 22;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];

toc
