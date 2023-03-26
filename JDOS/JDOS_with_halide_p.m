% Momentum matrix element (MME) calculation for halide double perovskite (Cs2BB'X6,
% B = Na/K, B' = In/Bi, and X = Cl).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
tic
addpath(genpath([fileparts(fileparts(pwd)), filesep]));
%% Tight-binding parameters
%% CsNaInCl, CsNaBiCl, CsKInCl,   CsKBiCl
EsB	  = [12.4, 9,   15,  15];
EpB	  = [13.5, 10,  15,  15];
VssB  = [0.05, 0,   0.0,   0];
VspB  =	[0.05, 0,   0.0,   0];
VppsB =	[0.3, 0.0,  -0.0,   0];
VpppB =	[0.0, 0.0,  0.0,   0];

EsB1   = [5.9,   -0.02,    6.21,  0.032];
EpB1   = [9.64,  5.23,    9.35,  5.12];
VssB1  = [-0.022, -0.01,  0.015,  0.015];
VspB1  = [ 0.018,  0.0,     -0.02,  0.02];
VppsB1 = [-0.152,  -0.03,  -0.18,  0.04];
VpppB1 = [-0.00, -0.00,  0.0,  -0.02];

VssBB   = [-0.24, -0.3, -0.36,  -0.44];
VspBB   = [ 0.16,  0.3,  0.12,  0.2];
VppsBB  = [	0.1,  0.40,  0.2,  0.2];
VpppBB  = [	0.0,  0.0,  0.0,  0.0];
% SOC strength
sB  = [0, 0, 0, 0];
sB1 = [0.16, 0.58, 0.16, 0.6];

% halide p strength
ECs     = [-2,   -2,   -4,    -5];
ECp     = [-0.1, -0.95, -0.15,  -1.1];
VppsMH1 = [-1.4, -1.3,  -1.4,  -1.4];
VpppMH1 = [0.7,  0.7,  0.7,   0.7];
VppsMH2 = [1.6,  0.6,  1.4,   0.5];
VpppMH2 = [0.4,  0.4,  0.7,   0.7];
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt =  202;
nbnd =  [276, 284, 294, 308];
% Select the compound
%%
p = 1;

t1 = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];
t2 = [VssB1(p),VspB1(p),VppsB1(p),VpppB1(p),EpB1(p),EsB1(p)];
t3 = [VssBB(p),VspBB(p),VppsBB(p),VpppBB(p),0,0];

t4 = [ECs(p), ECp(p), VppsMH1(p), VpppMH1(p), VppsMH2(p), VpppMH2(p)];

%%
nbasis = 64;

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
ndiv = 20;
% desigh the kmesh
[X,Y,Z] = meshgrid(linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv));

kpath = [X(:),Y(:),Z(:)];
%% Desigh the Hamiltonian and calculate the eigenvalues
%%
sig = 0.1;

momentum = [];

parfor E1 = 1:400
    ss = 1;
    jdos = [];
    E = 3+E1/100;
    for istate = 47:48       % VB
        for fstate = 49:50   % CB
            % diagonolize the Hamiltonian
            for i = 1:length(kpath)
                kx = kpath(i,1)*2*pi;
                ky = kpath(i,2)*2*pi;
                kz = kpath(i,3)*2*pi;
                % get the Hamiltonian matrix
                H0 = full_Hamiltonian_BpXp(kx,ky,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p));
                % get the eigenvalues
                [v,e] = eig(H0);
                % Calculation of MME
                dk = 2*pi/(ndiv-1);
                
                dHx = (full_Hamiltonian_BpXp(kx+dk,ky,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)) - full_Hamiltonian_BpXp(kx-dk,ky,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)))/(2*dk);
                dHy = (full_Hamiltonian_BpXp(kx,ky+dk,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)) - full_Hamiltonian_BpXp(kx,ky-dk,kz,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)))/(2*dk);
                dHz = (full_Hamiltonian_BpXp(kx,ky,kz+dk,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)) - full_Hamiltonian_BpXp(kx,ky,kz-dk,t1,t2,t3,t4,posAA,posAB,sB(p),sB1(p)))/(2*dk);
                %%
                px = transpose(conj(v(:,fstate)))*dHx*v(:,istate);
                py = transpose(conj(v(:,fstate)))*dHy*v(:,istate);
                pz = transpose(conj(v(:,fstate)))*dHz*v(:,istate);
                % Calculate the joint density of states and epsilon_2(omega)
                MME = (px*conj(px)+py*conj(py)+pz*conj(pz));
                
                jdos(ss) = MME*(2/(sig*sqrt(2*pi))*exp(-((e(fstate,fstate)-e(istate,istate)-E)^2/(2*sig^2))));
                % store the output
                ss = ss+1;
            end
        end
    end
    
    sum_jdos(E1,:) = [E,sum(jdos)/ss];
    
end
% plot the E vs JDOS/epsilon_2
plot(sum_jdos(:,1),sum_jdos(:,2),'b','LineWidth',2)
%%  axis properties
ax=gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 22;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];
ax.YLim = [0,0.1];

toc
