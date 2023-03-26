function [Hp1soc,Hp2soc] = HSOC_pp(sA,sB)
% retun SOC Hamiltonian for p orbiral system where sA is SOC strength
% for B-p-orbitals and sB is the SOC strength of B'-p-orbitals. To run this
% code: [Hp1soc,Hp2soc] = HSOC_pp(sA,sB)
%%
% SOC Hamiltonian for B-{s,p}-orbitals in the base order {s, px, py, pz}-up,
% {s, px, py, pz}-dn.
Hp1soc  =   sA*[0,    0,    0,    0,   0,    0,   0,    0;
                0,    0,-1j*1,    0,   0,    0,   0,    1;
                0, 1j*1,    0,    0,   0,    0,   0,-1j*1;
                0,    0,    0,    0,   0,   -1,1j*1,    0;

                0,    0,    0,    0,   0,    0,   0,    0;
                0,    0,    0,   -1,   0,    0,1j*1,    0;
                0,    0,    0,-1j*1,   0,-1j*1,   0,    0;
                0,    1, 1j*1,    0,   0,    0,   0,    0];
% SOC Hamiltonian for B'-{s,p}-orbitals in the base order {s, px, py, pz}-up,
% {s, px, py, pz}-dn.
Hp2soc  =   sB*[0,    0,    0,    0,   0,    0,   0,    0;
                0,    0,-1j*1,    0,   0,    0,   0,    1;
                0, 1j*1,    0,    0,   0,    0,   0,-1j*1;
                0,    0,    0,    0,   0,   -1,1j*1,    0;

                0,    0,    0,    0,   0,    0,   0,    0;
                0,    0,    0,   -1,   0,    0,1j*1,    0;
                0,    0,    0,-1j*1,   0,-1j*1,   0,    0;
                0,    1, 1j*1,    0,   0,    0,   0,    0];
