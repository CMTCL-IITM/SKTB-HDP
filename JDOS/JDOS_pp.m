% Momentum matrix element (MME) calculation for halide double perovskite (Cs2BB'X6,
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
p = 6;
% hopping parameters
t1 = [VssA(p),VspA(p),VppsA(p),VpppA(p),EpA(p),EsA(p)];
t2 = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];
t3 = [VssAB(p),VspAB(p),VppsAB(p),VpppAB(p),0,0];
% size of basis set
nbasis = 16;
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
        -1  0 -1
         1  0 -1];

posAB = [1  0  0;
        -1  0  0;
         0 -1  0;
         0  1  0;
         0  0  1;
         0  0 -1];
% step size in k-vector
ndiv = 20;
% desigh the kmesh
[X,Y,Z] = meshgrid(linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv));
kpath = [X(:),Y(:),Z(:)];

%% Desigh the Hamiltonian and calculate the eigenvalues
sig = 0.1;
momentum = [];
parfor E1 = 1:600
    ss = 1;
    jdos = [];
    eps2 = [];
    E = E1/100;
    % initial and final states between which the the transition is happening.       
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
                % Calculation of MME
                %%
                dk = 2*pi/(ndiv-1);
                
                dHx = (full_Hamiltonian_pp(kx+dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pp(kx-dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
                dHy = (full_Hamiltonian_pp(kx,ky+dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pp(kx,ky-dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);
                dHz = (full_Hamiltonian_pp(kx,ky,kz+dk,t1,t2,t3,posAA,posAB,sA(p),sB(p))-full_Hamiltonian_pp(kx,ky,kz-dk,t1,t2,t3,posAA,posAB,sA(p),sB(p)))/(2*dk);

                px = transpose(conj(v(:,fstate)))*dHx*v(:,istate) ;
                py = transpose(conj(v(:,fstate)))*dHy*v(:,istate) ;
                pz = transpose(conj(v(:,fstate)))*dHz*v(:,istate) ;
                % Calculate the joint density of states and epsilon_2(omega)
                MME = px*conj(px)+py*conj(py)+pz*conj(pz);
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
