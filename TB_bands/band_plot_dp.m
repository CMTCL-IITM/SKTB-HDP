% Slater-Koster tight-binding model for halide double perovskite (Cs2BB'X6,
% B = Ag, B' = In/Bi/Sb/Tl, and X = Cl).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace
close all
clear
clc
tic
addpath(genpath([fileparts(fileparts(pwd)), filesep]));
%% Tight-binding parameters
%   	AgInCl	AgBiCl AgSbCl  AgTlCl
EsA =	[7.15	7.76    5.86  7.76];
Eeg =   [-0.08	-0.4    -1.1  0.04];
Ex2 =   [-0.08	-0.4    -1.1  0.04];
Et2g=	[-0.82	-0.91   -1.81  -1.16];
Vdds=	[-0.01	-0.01    0    0.016];
Vddp=	[0.007	0.007    0   -0.017];
Vddd=	[0.0	0.0      0     0];
VssA=	[0.0	0.0      0.03   0];
VsdA=	[0.0	0.0      0      0];

EsB =   [4.6	-1.85  -1.16   1.46];
EpB =   [7.98	4.18    3.41   6.81];
EpzB=   [7.98	4.18    3.41   6.81];
VssB= 	[0.0	-0.01   -0.065  0.048];
VspB =	[0.018	0.03    0.02   0.126];
VppsB =	[-0.092	-0.028  -0.076  -0.078];
VpppB =	[0.00	0.00     0.0   0.0];

VsAsB=	[0.375	0.0    0.0   0.532];
VspAB=	[0.48	0.475  0.57  0.54];
VsdAB=	[0.45	0.193   0.25  0.39];
Vpds=	[0.63	0.54    0.63  0.84];
Vpdp=	[0.24	0.0     0.14  0.2];
% SOC strength
sA  =   [0.09	0.09   0.09  0.09];
sB  =   [0.16	0.54   0.17  0.24];
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt =  202;
nbnd =  [ 282 256 292  292];
% Select the compound
p = 1; % Input
%%
if p == 1
    BS=load('../DFT_data/CsAgInCl.dat');
elseif p == 2
    BS=load('../DFT_data/CsAgBiCl.dat');
elseif p == 3
    BS=load('../DFT_data/CsAgSbCl.dat');
elseif p == 4
    BS=load('../DFT_data/CsAgTlCl.dat');
end
% plot the dft band structure
dft = reshape(BS(:,2),[nkpt,nbnd(p)]);
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
kpath=load('../DFT_data/klist.dat');
%% crystal structure matrix
a = 10.6;     % lattice constant
rvec = [a, 0, 0;
    0, a, 0;
    0, 0, a];
%% calculate reciprocal lattice vectors
Vol = dot(rvec(1,:),cross(rvec(2,:),rvec(3,:)));
b1 =  2*pi*cross(rvec(2,:),rvec(3,:))/Vol;
b2 =  2*pi*cross(rvec(3,:),rvec(1,:))/Vol;
b3 =  2*pi*cross(rvec(1,:),rvec(2,:))/Vol;
recip = [b1;b2;b3];
% position coordinates of 2nd-nearest and 4th-nearest neighbor B and B' atoms in crystal unit
% 4th-nearest neighbor
posAA = [1  1  0;
    -1  1  0;
    -1 -1  0;
    1 -1  0;
    0  1  1;
    0 -1 -1;
    0  1 -1;
    0 -1  1;
    1  0  1;
    -1  0  1;
    -1  0 -1;
    1  0 -1];
% 2nd-nearest neighbor
posAB = [1  0  0;
    -1  0  0;
    0 -1  0;
    0  1  0;
    0  0  1;
    0  0 -1];
%% Desigh the Hamiltonian and calculate the eigenvalues along the high-symmetry kpath.
% hopping parameters
t1 = [Eeg(p),Et2g(p),Vdds(p),Vddp(p),Vddd(p),VssA(p),EsA(p),VsdA(p),Eeg(p)];   % A-A
t2 = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p),EpB(p)];                 % B-B
t3 = [VsAsB(p),VspAB(p),VsdAB(p),Vpds(p),Vpdp(p)];                             % A-B
% diagonolize the Hamiltonian
for i= 1:length(kpath)
    kx = kpath(i,1)*2*pi;
    ky = kpath(i,2)*2*pi;
    kz = kpath(i,3)*2*pi;
    % get the Hamiltonian matrix
    H = full_Hamiltonian_pd(kx,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p));
    % get the eigenvalues
    v = eig(H);
    % store the eigenvalues
    eig_val(i,:) = v';
end
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
%% plot the DFT and TB band structure
figure(1)
plot(s,dft,'r','LineWidth',2)
hold on
plot(s,eig_val,'b','LineWidth',2);
hold off
%% axis properties
ax=gca;
ylim([-2,8]);
xlim([0 s(end)]);
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 22;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];

toc
