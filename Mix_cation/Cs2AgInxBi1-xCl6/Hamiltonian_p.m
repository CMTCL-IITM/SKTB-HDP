function Hpp = Hamiltonian_p(kx,ky,kz,t,pos,a)
%
nat = size(pos,1);
% load tb parameters
Vss = t(1);Vsp=t(2);Vpps=t(3);Vppp=t(4);Opp=t(5);Oss=t(6);

%%  dcs
for i = 1: size(pos,1)
    alpha = pos(i,1)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    beta = pos(i,2)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    gamma = pos(i,3)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    dcs(i,:) = [alpha,beta,gamma];
end

ss = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    ss  = ss + Vss*exp(1j*(dot(k,r)));
end
ss = real(ss);

spx = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    spx  = spx + l*Vsp*exp(1j*(dot(k,r)));
end

spy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    spy  = spy + m*Vsp*exp(1j*(dot(k,r)));
end

spz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    spz  = spz + n*Vsp*exp(1j*(dot(k,r)));
end

xx = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    xx  = xx + ((l^2)*Vpps+(1-l^2)*Vppp)*exp(1j*(dot(k,r)));
end
xx = real(xx);

yy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    yy = yy + ((m^2)*Vpps+(1-m^2)*Vppp)*exp(1j*(dot(k,r)));
end
yy = real(yy);

xy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    xy  = xy + l*m*(Vpps-Vppp)*exp(1j*(dot(k,r)));
end

zz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];
    zz = zz + ((n^2)*Vpps+(1-n^2)*Vppp)*exp(1j*(dot(k,r)));
end
zz = real(zz);

%% basis order -- s, px, py, pz

Hpp =  [Oss+ss,    spx,spy,spz;
        conj(spx), Opp+xx, xy, 0;
        conj(spy), conj(xy), Opp+yy,0;
        conj(spz), 0, 0, Opp+zz];

