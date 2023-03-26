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
tic
addpath(genpath([fileparts(fileparts(pwd)), filesep]));
%% Tight-binding parameters
% CsNaInCl, CsNaBiCl, CsKInCl,   CsKBiCl
EsB	  = [12.4, 9,   15,  15];
EpB	  = [13.5, 10,  15,  15];
VssB 	  = [0.0, 0,   0.0,   0];
VspB  	  = [0.0, 0,   0.0,   0];
VppsB 	  = [0.0, 0.0, -0.0,  0];
VpppB     = [0.0, 0.0,  0.0,  0];

EsB1   =  [5.8,   -0.08,    6.21,  0.0];
EpB1   =  [9.64,    5.2,    9.35,  5.12];
VssB1  =  [-0.036, -0.034,   0.01,  0.01];
VspB1  =  [ 0.018,   0.05,   0.05,  0.04];
VppsB1 =  [-0.152,  -0.08,  -0.15, -0.01];
VpppB1 =  [-0.00,   0.02,    0.0,  0.0]; 
 
VssBB   = [-0.0, -0.0,-0.33, -0.4];
VspBB   = [ 0.0,  0.0,  0.0,  0.0];
VppsBB  = [	0.0,  0.0,  0.0,  0.0];
VpppBB  = [	0.0,  0.0,  0.0,  0.0];

ECs     = [-2,   -2,   -4,    -5];
ECp     = [-0.1, -0.95, -0.15,  -1.1];
VppsMH1 = [-1.4, -1.3,  -1.4,  -1.4];
VpppMH1 = [0.7,  0.7,  0.7,   0.7];
VppsMH2 = [1.6,  0.6,  1.4,   0.5];
VpppMH2 = [0.4,  0.4,  0.7,   0.7];
% SOC strength
sB  = [0, 0, 0, 0];
sB1 = [0.0, 0.58, 0.0, 0.0];
%%
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt =  202;
nbnd =  [276, 284, 294, 308];
p = 2;
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
kpath = load('../DFT_data/klist.dat');
% hopping parameters
t1 = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];
t2 = [VssB1(p),VspB1(p),VppsB1(p),VpppB1(p),EpB1(p),EsB1(p)];
t3 = [VssBB(p),VspBB(p),VppsBB(p),VpppBB(p),0,0];
t4 = [ECs(p), ECp(p), VppsMH1(p), VpppMH1(p), VppsMH2(p), VpppMH2(p)];
% size of basis set
nbasis = 64;
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
for istate = 49:50       % VB
    for fstate = 51:52   % CB
        % diagonolize the Hamiltonian
        for i = 1:length(kpath)
            
            kx = kpath(i,1)*2*pi;
            ky = kpath(i,2)*2*pi;
            kz = kpath(i,3)*2*pi;
            % get the Hamiltonian matrix
            H0 = full_Hamiltonian_BpXp(kx,ky,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p));
            % get the eigenvalues
            [v,e] = eig(H0);
            % store the eigenvalues
            eig_val(i,:) = diag(e)';
            % Calculation of MME
            dk = 2*pi/(ndiv-1);
            dHx = (full_Hamiltonian_BpXp(kx+dk,ky,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p))-full_Hamiltonian_BpXp(kx-dk,ky,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)))/(2*dk);
            dHy = (full_Hamiltonian_BpXp(kx,ky+dk,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p))-full_Hamiltonian_BpXp(kx,ky-dk,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)))/(2*dk);
            dHz = (full_Hamiltonian_BpXp(kx,ky,kz+dk,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p))-full_Hamiltonian_BpXp(kx,ky,kz-dk,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)))/(2*dk);
            %%
            Hp = [dHx,dHy,dHz];
            
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
