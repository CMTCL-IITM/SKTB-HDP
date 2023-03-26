function HMC = full_Hamiltonian_pp(kx,ky,kz,t1,t2,t3,t4,posAA,posAB,sA,sB)
% Import kpath, TB parameters, and position of atoms from
% band_plot_pristine.m file. To run this code use:
% H = full_Hamiltonian_pp(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB).
% As an output it gives a 80X80 Hamiltonian (with SOC included)
% of d-p orbitals-based HDP.
%% Design the SK-TB Hamiltonian
% For s-p orbitals interaction (B-B)
    HAA  = Hamiltonian_p(kx,ky,kz,t1,posAA,1/sqrt(2));
    HBB  = Hamiltonian_p(kx,ky,kz,t2,posAA,1/sqrt(2));
% For s-p orbitals interaction (B-B')   
    HAB  = Hamiltonian_p(kx,ky,kz,t3,posAB,1/2);
    HBA = HAB';
% Hamiltonian without SOC
    Htb  = [HAA, zeros(4,4), HAB, zeros(4,4);
            zeros(4,4), HAA, zeros(4,4), HAB;
            HBA, zeros(4,4), HBB, zeros(4,4);
            zeros(4,4), HBA, zeros(4,4), HBB];
% Hamiltonian of SOC  
    [Hpsoc,Hp1soc] = HSOC_pp(sA,sB);
    Hsoc = [Hpsoc,zeros(8,8);
            zeros(8,8),Hp1soc];
% B-B' Hamiltonian    
    HM = Htb+Hsoc;
% Hamiltonian for X-p - {B,B'} interaction
HMC = Hamiltonian_MH(kx,ky,kz,t4,posAB);
HMC(1:16,1:16) =  HM;
 
