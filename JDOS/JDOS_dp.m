% Momentum matrix element (MME) and joint densty of states (JDOS) calculation for halide double perovskite (Cs2BB'X6,
% B = Ag, B' = In/Bi/Sb/Tl, and X = Cl).
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
%  	AgInCl	AgBiCl AgSbCl  AgTlCl
EsA =	[7.15	7.56    5.86  7.76];
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

VsAsB=	[0.375	0.27    0.0   0.532];
VspAB=	[0.48	0.4755  0.57  0.54];
VsdAB=	[0.45	0.193   0.25  0.39];
Vpds=	[0.63	0.54    0.63  0.84];
Vpdp=	[0.24	0.0     0.14  0.2];
% SOC strength
sA  =   [0.09	0.09   0.09  0.09];
sB  =   [0.16	0.54   0.17  0.24];
%%
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt =  202;
nbnd =  [ 282 256 292  292];
% Select the compound
p = 4;
% hopping parameters
t1 = [Eeg(p),Et2g(p),Vdds(p),Vddp(p),Vddd(p),VssA(p),EsA(p),VsdA(p),Eeg(p)];   % A-A
t2 = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p),EpB(p)];                 % B-B
t3 = [VsAsB(p),VspAB(p),VsdAB(p),Vpds(p),Vpdp(p)];
% size of basis set
nbasis = 20;
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
% desigh the kmesh
[X,Y,Z] = meshgrid(linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv));
kpath = [X(:),Y(:),Z(:)];
%% Desigh the Hamiltonian and calculate the eigenvalues
sig = 0.1;
momentum = [];
parfor E1 = 1:300
    ss = 1;
    jdos = [];
    eps2 = [];
    E = 0+E1/100;
    % initial and final states between which the the transition is happening.
    for istate = 9:10       % VB
        for fstate = 11:12   % CB
            % diagonolize the Hamiltonian
            for i = 1:length(kpath)
                
                kx = kpath(i,1)*2*pi;
                ky = kpath(i,2)*2*pi;
                kz = kpath(i,3)*2*pi;
                % get the Hamiltonian matrix
                H0 = full_Hamiltonian_pd(kx,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p));
                % get the eigenvalues
                [v,e] = eig(H0);
                % Calculation of MME
                dk = 2*pi/(ndiv-1);
                dHx = (full_Hamiltonian_pd(kx+dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pd(kx-dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
                dHy = (full_Hamiltonian_pd(kx,ky+dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pd(kx,ky-dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
                dHz = (full_Hamiltonian_pd(kx,ky,kz+dk,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pd(kx,ky,kz-dk,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
                
                px = transpose(conj(v(:,fstate)))*dHx*v(:,istate) ;
                py = transpose(conj(v(:,fstate)))*dHy*v(:,istate) ;
                pz = transpose(conj(v(:,fstate)))*dHz*v(:,istate) ;
                
                MME = (px*conj(px)+py*conj(py)+pz*conj(pz));
                % Calculate the joint density of states and epsilon_2(omega)
                jdos(ss) = (2/(sig*sqrt(2*pi))*exp(-((e(fstate,fstate)-e(istate,istate)-E)^2/(2*sig^2))));
                eps2(ss) = MME*(2/(sig*sqrt(2*pi))*exp(-((e(fstate,fstate)-e(istate,istate)-E)^2/(2*sig^2))));
                ss = ss+1;
            end
        end
    end
    % store the output
    sum_jdos(E1,:) = [E,sum(jdos)/ss,sum(eps2)/ss];
end
% plot the E vs JDOS/epsilon_2
figure(1)
plot(sum_jdos(:,1),sum_jdos(:,2),'r','LineWidth',2)
figure(2)
plot(sum_jdos(:,1),sum_jdos(:,3),'b','LineWidth',2)
%% axis properties
ax=gca;

ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 22;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];
ax.YLim = [0,0.3];

toc
