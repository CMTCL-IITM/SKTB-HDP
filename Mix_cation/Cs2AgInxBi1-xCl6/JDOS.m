% JDOS and MME for mixed cation halide double perovskite (Cs2AgInxBi1-xCl6).
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for p = 1:2

        t1(p,:) = [Eeg(p),Et2g(p),Vdds(p),Vddp(p),Vddd(p),VssA(p),EsA(p),VsdA(p)];   % A-A
        t2(p,:) = [VssB(p),VspB(p),VppsB(p),VpppB(p),EpB(p),EsB(p)];                 % B-B
        t3(p,:) = [VsAsB(p),VspAB(p),VsdAB(p),Vpds(p),Vpdp(p)];                      % A-B  
end

%%
nbasis = 80;

%%
nkpt = 203;
nbnd =  999;

%%
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

%%
ndiv = 10;

[X,Y,Z] = meshgrid(linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv),linspace(-pi,pi,ndiv));

kpath = [X(:),Y(:),Z(:)];

%%
sig = 0.1;

momentum = [];

parfor E1 = 1:400
  
    ss = 1;
    jdos = [];
    eps2 = [];
    E = 1+E1/100;
        
    for istate = 43:44       % VB
        for fstate = 45:48   % CB
            
            for i = 1:length(kpath)
                
                kx = kpath(i,1)*2*pi;
                ky = kpath(i,2)*2*pi;
                kz = kpath(i,3)*2*pi;
                
                H0 = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd);
                
                [v,e] = eig(H0);
                
                %%
                dk = 2*pi/(ndiv-1);
                
                dHx = (full_Hamiltonian_pd_sc(kx+dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd)-full_Hamiltonian_pd_sc(kx-dk,ky,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd))/(2*dk);
                dHy = (full_Hamiltonian_pd_sc(kx,ky+dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd)-full_Hamiltonian_pd_sc(kx,ky-dk,kz,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd))/(2*dk);
                dHz = (full_Hamiltonian_pd_sc(kx,ky,kz+dk,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd)-full_Hamiltonian_pd_sc(kx,ky,kz-dk,t1,t2,t3,posAA,posAB,sA(p),sB(p),selectcompd))/(2*dk);
                %%
                Hp = [dHx,dHy,dHz];
                
                %%
                
                intraatomic = 0.00*1j*(e(fstate,fstate)-e(istate,istate))*(transpose(conj(v(:,fstate)))*(ones(nbasis)-eye(nbasis))*v(:,istate));
                
                px = transpose(conj(v(:,fstate)))*dHx*v(:,istate) + intraatomic;
                py = transpose(conj(v(:,fstate)))*dHy*v(:,istate) + intraatomic;
                pz = transpose(conj(v(:,fstate)))*dHz*v(:,istate) + intraatomic;
                
               MME = (px*conj(px)+py*conj(py)+pz*conj(pz));
               
               jdos(ss) = (2/(sig*sqrt(2*pi))*exp(-((e(fstate,fstate)-e(istate,istate)-E)^2/(2*sig^2))));
           
               eps2(ss) = MME*(2/(sig*sqrt(2*pi))*exp(-((e(fstate,fstate)-e(istate,istate)-E)^2/(2*sig^2))));
              
               ss = ss+1;
            
            end
        end
    end
    
    sum_jdos(E1,:) = [E,sum(jdos)/ss,sum(eps2)/ss];
    
end
figure(1)
plot(sum_jdos(:,1),sum_jdos(:,2),'r','LineWidth',3)
figure(2)
plot(sum_jdos(:,1),sum_jdos(:,3),'b','LineWidth',3)
%%
ax=gca;

ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 22;
ax.TickDir = 'in';
ax.TickLength = [0.01 0.01];
ax.YLim = [0,0.1];
ax.XLim = [1,5];

toc
