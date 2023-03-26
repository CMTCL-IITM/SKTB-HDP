function H = Hamiltonian_MH(kx,ky,kz,t4,posAB)
% Import kpath, TB parameters, lattice constant, and position of atoms.
% To run this code use:
% Hpp = Hamiltonian_p(kx,ky,kz,t,pos,a).
% As an output it reture a 24X24 Hamiltonian of B-{s,p}-X-{s,p} orbitals interactions.
x = 0.25;
H3 = zeros(16,48);
% obtain fitting parameters
Ecs = t4(1);
Ecp = t4(2);
t41 = [0.0,0.0,t4(3),t4(4),0,0];
t42 = [0.0,0.0,t4(5),t4(6),0,0];
%% Design Hamiltonian for X-p orbitals
for i = 1:6

H3(1:4,4*i-3:4*i)  = Hamiltonian_p_for_H(kx,ky,kz,t41,posAB(i,:),x);
H3(5:8,4*i+21:4*i+24)  = Hamiltonian_p_for_H(kx,ky,kz,t41,posAB(i,:),x);

H3(9:12,4*i-3:4*i)  = Hamiltonian_p_for_H(kx,ky,kz,t42,-posAB(i,:),0.5-x);
H3(13:16,4*i+21:4*i+24)  = Hamiltonian_p_for_H(kx,ky,kz,t42,-posAB(i,:),0.5-x);

end

Hcl = Hamiltonian_p_for_H(kx,ky,kz,[0.00, 0.0, 0.04, 0.02, Ecp, Ecs],posAB,1/2);
H4 = kron(eye(12),Hcl);

H = [zeros(16),H3;
     transpose(conj(H3)),H4];