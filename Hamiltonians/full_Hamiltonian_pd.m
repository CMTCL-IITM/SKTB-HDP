function H = full_Hamiltonian_pd(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB)
% Import kpath, TB parameters, and position of atoms from
% band_plot_pristine.m file. To run this code use:
% H = full_Hamiltonian_pd(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB).
% As an output it gives a 20X20 Hamiltonian (with SOC included)
% of d-p orbitals-based HDP.
%% Design the SK-TB Hamiltonian
% For s-d orbitals interaction (B-B)
HAA  = Hamiltonian_d(kx,ky,kz,t1,posAA,1/sqrt(2));
% For s-p orbitals (B'-B')
HBB  = Hamiltonian_p(kx,ky,kz,t2,posAA,1/sqrt(2));
% For {s,p}-{s,d} orbitals interaction s-p orbitals (B-B')
HAB  = Hamiltonian_spsd(kx,ky,kz,t3,posAB,1/2);
HBA = HAB';
% Hamiltonian without SOC
Htb  = [HAA, zeros(6,6), HAB, zeros(6,4);
        zeros(6,6), HAA, zeros(6,4), HAB;
        HBA, zeros(4,6), HBB, zeros(4,4);
        zeros(4,6), HBA, zeros(4,4), HBB];
% Hamiltonian of SOC
[Hdsoc,Hpsoc] = HSOC_pd(sA,sB);
Hsoc = [Hdsoc,zeros(12,8);
        zeros(8,12),Hpsoc];
% Full Hamiltonian
H = Htb+Hsoc;
