function H = full_Hamiltonian_pp(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB)
% Import kpath, TB parameters, and position of atoms from
% band_plot_pristine.m file. To run this code use:
% H = full_Hamiltonian_pp(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB).
% As an output it gives a 8X8 Hamiltonian (with SOC included)
% of d-p orbitals-based HDP.
%% Design the SK-TB Hamiltonian
% For s-p orbitals interaction (B-B)
    HAA     = Hamiltonian_p(kx,ky,kz,t1,posAA,1/sqrt(2));
% For s-p orbitals interaction (B'-B')
    HBB     = Hamiltonian_p(kx,ky,kz,t2,posAA,1/sqrt(2));
% For s-p orbitals interaction (B-B')
    HAB     = Hamiltonian_p(kx,ky,kz,t3,posAB,1/2);
% Hamiltonian without SOC
    HAA = [HAA,zeros(4);zeros(4),HAA];
    HBB = [HBB,zeros(4);zeros(4),HBB];
    HAB = [HAB,zeros(4);zeros(4),HAB];
    HBA     = transpose(conj(HAB));
    
    Htb = [HAA,HAB;
           HBA,HBB];
% Hamiltonian of SOC     
    [Hp1soc,Hp2soc] = HSOC_pp(sA,sB);
    Hsoc = [Hp1soc,zeros(8);
          zeros(8),Hp2soc];
% Full Hamiltonian
    H = Htb+Hsoc;
