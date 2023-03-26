% MME for mixed cation halide double perovskite (Cs2AgInxBi1-xCl6).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
% Note: check istate and fstate from eig_val before each MME calculations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear the workspace
clear
clc
%% Tight-binding parameters
selectcompd = 2;

%%  AgInCl	AgBiCl for In0.25Bi0.75
if selectcompd == 1
EsA =	[6.36	6.36];
Eeg	=   [-0.67	-0.67];
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
% load the dft kpoints (kx,ky,kz) in the unit of 1/scale.
kpath=load('klist_sc.dat');
% other variables: nkpt = no. of dft k-points, nbnd = total no. of dft
% bands
nkpt = 203;
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
posAA = [1 1 0;
        -1 -1 0;
        -1 1 0;
        1 -1 0;
        0 1 1;
        0 -1 -1;
        0 1 -1;
        0 -1 1;
        1 0 1;
        -1 0 -1
        -1 0 1;
        1 0 -1];

posAB = [1 0 0;
        -1 0 0;
        0 -1 0;
        0  1 0;
        0 0 1;
        0 0 -1];
% Desigh the Hamiltonian and calculate the eigenvalues along the high-symmetry kpath.
% hopping parameters
energy=[];
for p = 1:2
    t1(p,:) = [Eeg(p),Et2g(p),Vdds(p),Vddp(p),Vddd(p),VssA(p),EsA(p),VsdA(p)];   % A-A
    t2(p,:) = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];                 % B-B
    t3(p,:) = [VsAsB(p),VspAB(p),VsdAB(p),Vpds(p),Vpdp(p)];                      % A-B
end
%% diagonolization of the Hamiltonian
ndiv = 10;
nbasis = 80;
momentum = [];
ss = 1;
% initial and final states between which the the transition is happening.
% change istart and fstart for Cs2AgInCl6 and Cs2AgTlCl6 (check the
% eig_val data for VBE and CBE).
for istate = 43:44       % VB
    for fstate = 45:48   % CB
        % diagonolize the Hamiltonian
        for i = 1:length(kpath)
            kx = kpath(i,1)*2*pi;
            ky = kpath(i,2)*2*pi;
            kz = kpath(i,3)*2*pi;
            % get the Hamiltonian matrix
            H0 = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd);
            % store the eigenvalues
            [v,e] = eig(H0);
            [energy1,perm] = sort(real(diag(e)'));
            energy(i,:) = energy1;
            vectors(:,:) = v(:,perm);
            % Calculation of MME
            dk = 2*pi/(ndiv-1);
            
            dHx = (full_Hamiltonian_pd_sc(kx+dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd)-full_Hamiltonian_pd_sc(kx-dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd))/(2*dk);
            dHy = (full_Hamiltonian_pd_sc(kx,ky+dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd)-full_Hamiltonian_pd_sc(kx,ky-dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd))/(2*dk);
            dHz = (full_Hamiltonian_pd_sc(kx,ky,kz+dk,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd)-full_Hamiltonian_pd_sc(kx,ky,kz-dk,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd))/(2*dk);

            px = transpose(conj(vectors(:,fstate)))*dHx*vectors(:,istate) ;
            py = transpose(conj(vectors(:,fstate)))*dHy*vectors(:,istate) ;
            pz = transpose(conj(vectors(:,fstate)))*dHz*vectors(:,istate) ;
            
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
plot(s,energy,'r','LineWidth',2);
figure(2)
plot(s,prefact*sum(momentum(:,:,end),2)/(nkpt*ss),'k','LineWidth',2)
%% axis properties
ax = gca;
ax.Box = 'on';
ax.TickDir = 'in';
ax.TickLength = [0.0001,0.0001];
ax.LineWidth = 2;
ax.XLim = [0 s(end)];
ax.YLim = [0,2];

toc
