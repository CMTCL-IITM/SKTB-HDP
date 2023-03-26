function Hdd = Hamiltonian_d(kx,ky,kz,t,pos,a) 
% Import kpath, TB parameters, lattice constant, and position of atoms.
% To run this code use:
% Hdd = Hamiltonian_d(kx,ky,kz,t,pos,a).
% As an output it return a 6X6 Hamiltonian of {s,d}-{s,d} orbitals interactions.
%% obtain fitting parameters
Oz2=t(1);Ot2g=t(2);Vdds=t(3);Vddp=t(4);Vddd=t(5);VssA=t(6);OssA=t(7);VsdA=t(8);
% nat = no. of atoms
nat = size(pos,1);
%  generate direction cosines (dcs) from atomic positions
for i = 1: size(pos,1)
    alpha = pos(i,1)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    beta = pos(i,2)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    gamma = pos(i,3)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    dcs(i,:) = [alpha,beta,gamma];
end
%% d-d interactions
x2x2 = 0;
for j = 1:nat
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    
    x2x2 = x2x2 + (3/4*(l^2-m^2)^2*Vdds+(l^2+m^2-(l^2-m^2)^2)*Vddp+(n^2+(l^2-m^2)^2/4)*Vddd)*exp(1j*(dot(k,r)));
end
x2x2 = real(x2x2);

z2z2 = 0;
for j = 1:nat
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2z2 = z2z2 + ((n^2-(l^2+m^2)/2)*Vdds + 3*n^2*(l^2+m^2)*Vddp + 3/4*(l^2+m^2)^2*Vddd)*exp(1j*(dot(k,r)));
end
z2z2 = real(z2z2);

z2x2 = 0;
for j = 1:nat
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2x2 = z2x2 + sqrt(3)*((l^2-m^2)*(n^2-(l^2+m^2)/2)*Vdds/2+n^2*(m^2-l^2)*Vddp+(1+n^2)*(l^2-m^2)/4*Vddd)*exp(1j*(dot(k,r)));
end

x2xy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2xy = x2xy + (3/2*l*m*(l^2-m^2)*Vdds + 2*l*m*(m^2-l^2)*Vddp + l*m*(l^2-m^2)/2*Vddd)*exp(1j*(dot(k,r)));
end

x2xz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2xz = x2xz + (3/2*n*l*(l^2-m^2)*Vdds + n*l*(1-2*(l^2-m^2))*Vddp - n*l*(1-(l^2-m^2)/2)*Vddd)*exp(1j*(dot(k,r)));
end

x2yz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2yz = x2yz + (3/2*m*n*(l^2-m^2)*Vdds - m*n*(1+2*(l^2-m^2))*Vddp + m*n*(1+(l^2-m^2)/2)*Vddd)*exp(1j*(dot(k,r)));
end

z2xy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2xy = z2xy + sqrt(3)*(l*m*(n^2-(l^2+m^2)/2)*Vdds - 2*l*m*n^2*Vddp + l*m*(1+n^2)/2*Vddd)*exp(1j*(dot(k,r)));
end

z2xz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2xz = z2xz + sqrt(3)*(l*n*(n^2-(l^2+m^2)/2)*Vdds+l*n*(l^2+m^2-n^2)*Vddp-l*n*(l^2+m^2)/2*Vddd)*exp(1j*(dot(k,r)));
end

z2yz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2yz = z2yz + sqrt(3)*(m*n*(n^2-(l^2+m^2)/2)*Vdds+m*n*(l^2+m^2-n^2)*Vddp-m*n*(l^2+m^2)/2*Vddd)*exp(1j*(dot(k,r)));
end

xyxy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xyxy = xyxy + (3*l^2*m^2*Vdds+(l^2+m^2-4*l^2*m^2)*Vddp+(n^2+l^2*m^2)*Vddd)*exp(1j*(dot(k,r)));
end
xyxy = real(xyxy);

xzxz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzxz = xzxz + (3*l^2*n^2*Vdds+(l^2+n^2-4*l^2*n^2)*Vddp+(m^2+l^2*n^2)*Vddd)*exp(1j*(dot(k,r)));
end
xzxz = real(xzxz);

yzyz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    yzyz = yzyz + (3*m^2*n^2*Vdds+(m^2+n^2-4*m^2*n^2)*Vddp+(l^2+m^2*n^2)*Vddd)*exp(1j*(dot(k,r)));
end
yzyz = real(yzyz);

xyyz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xyyz = xyyz + (3*l*m^2*n*Vdds+l*n*(1-4*m^2)*Vddp+l*n*(m^2-1)*Vddd)*exp(1j*(dot(k,r)));
end

xyxz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xyxz = xyxz + (3*l^2*m*n*Vdds+m*n*(1-4*l^2)*Vddp+m*n*(l^2-1)*Vddd)*exp(1j*(dot(k,r)));
end

xzyz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzyz = xzyz + (3*l*m*n^2*Vdds+l*m*(1-4*n^2)*Vddp+l*m*(n^2-1)*Vddd)*exp(1j*(dot(k,r)));
end

%% s-d interactions
ss = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    ss = ss + VssA*exp(1j*(dot(k,r)));
end
ss = real(ss);

z2s = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2s = z2s + (n^2-(l^2+m^2)/2)*VsdA*exp(1j*(dot(k,r)));
end

x2s = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2s = x2s + (sqrt(3)/2*(l^2-m^2)*VsdA)*exp(1j*(dot(k,r)));
end

xys = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xys = xys + (sqrt(3)*l*m*VsdA)*exp(1j*(dot(k,r)));
end

xzs = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzs = xzs + (sqrt(3)*l*n*VsdA)*exp(1j*(dot(k,r)));
end

yzs = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    yzs = yzs + (sqrt(3)*n*m*VsdA)*exp(1j*(dot(k,r)));
end

%% basis order of the Hamiltonian -- {s, x2, z2, xy, xz, yz}
% Hamiltonian
Hdd =  [OssA+ss,conj(z2s),conj(x2s),conj(xys),conj(xzs),conj(yzs);
        z2s,    Oz2+z2z2,       z2x2, z2xy,  z2xz,  z2yz;
        x2s,conj(z2x2), Oz2+x2x2, x2xy,  x2xz,  x2yz;
        xys,conj(z2xy), conj(x2xy), Ot2g+xyxy,  xyxz,  xyyz;
        xzs,conj(z2xz), conj(x2xz), conj(xyxz),  Ot2g+xzxz,  xzyz;
        yzs,conj(z2yz), conj(x2yz), conj(xyyz),  conj(xzyz),  Ot2g+yzyz];
