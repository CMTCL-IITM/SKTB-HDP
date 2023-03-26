function [Hdsoc,Hpsoc] = HSOC_pd(sA,sB)
% retun SOC Hamiltonian for d and p orbiral system where sA is SOC strength
% for d-orbitals and sB is the SOC strength of p-orbitals. To run this
% code: [Hdsoc,Hpsoc] = HSOC_pd(sA,sB)
%%
% SOC Hamiltonian for {s,d}-orbitals in the base order {s, x2, z2, xy, xz,
% yz}-up, {s, x2, z2, xy, xz, yz}-dn.
Hdsoc    =  sA*[0,    0,   0,    0,      0,          0,    0,         0,    0,   0,       0,         0;
                0,    0,   0,    0,      0,          0,    0,         0,    0,   0,-sqrt(3),1j*sqrt(3);
                0,    0,   0,-2*1j,      0,          0,    0,         0,    0,   0,       1,        1j;
                0,    0,2*1j,    0,      0,          0,    0,         0,    0,   0,     -1j,         1;
                0,    0,   0,    0,      0,        -1j,    0,   sqrt(3),   -1,   1j,      0,         0;
                0,    0,   0,    0,     1j,          0,    0,-1j*sqrt(3),  -1j,  -1,      0,         0;
                0,    0,   0,    0,      0,          0,    0,         0,    0,   0,       0,         0;
                0,    0,   0,    0,sqrt(3), 1j*sqrt(3),    0,         0,    0,   0,       0,         0;
                0,    0,   0,    0,     -1,         1j,    0,         0,    0,2*1j,       0,         0;
                0,    0,   0,    0,    -1j,         -1,    0,         0,-2*1j,   0,       0,         0;
                0,-sqrt(3),1,   1j,      0,          0,    0,         0,    0,   0,       0,        1j;
                0,-1j*sqrt(3), -1j, 1,   0,          0,    0,         0,    0,   0,     -1j,        0];

% SOC Hamiltonian for {s,p}-orbitals in the base order {s, px, py, pz}-up,
% {s, px, py, pz}-dn.
Hpsoc    =  sB*[0,    0,    0,    0,   0,    0,   0,    0;
                0,    0,-1j*1,    0,   0,    0,   0,    1;
                0, 1j*1,    0,    0,   0,    0,   0,-1j*1;
                0,    0,    0,    0,   0,   -1,1j*1,    0;

                0,    0,    0,    0,   0,    0,   0,    0;
                0,    0,    0,   -1,   0,    0,1j*1,    0;
                0,    0,    0,-1j*1,   0,-1j*1,   0,    0;
                0,    1, 1j*1,    0,   0,    0,   0,    0];