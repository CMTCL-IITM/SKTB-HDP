function Hspsd = Hamiltonian_spsd(kx,ky,kz,t,pos,a)

Vss=t(1);Vsp=t(2);Vsd=t(3);Vpds=t(4);Vpdp=t(5);
%
nat = size(pos,1);

%%  dcs
for i = 1: size(pos,1)
    alpha = pos(i,1)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    beta = pos(i,2)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    gamma = pos(i,3)/sqrt(pos(i,1)^2+pos(i,2)^2+pos(i,3)^2);
    dcs(i,:) = [alpha,beta,gamma];
end

%%  d-p interaction

z2x = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2x = z2x + (l*(n^2-(l^2+m^2)/2)*Vpds-sqrt(3)*l*n^2*Vpdp)*exp(1j*(dot(k,r)));
end


z2y = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2y = z2y + (m*(n^2-(l^2+m^2)/2)*Vpds-sqrt(3)*m*n^2*Vpdp)*exp(1j*(dot(k,r)));
end


z2z = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2z = z2z + (n*(n^2-(l^2+m^2)/2)*Vpds + sqrt(3)*n*(l^2+m^2)*Vpdp)*exp(1j*(dot(k,r)));
end


x2x = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2x = x2x + (sqrt(3)/2*l*(l^2-m^2)*Vpds+l*(1-l^2+m^2)*Vpdp)*exp(1j*(dot(k,r)));
end



x2y = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2y = x2y + (sqrt(3)/2*m*(l^2-m^2)*Vpds-m*(1+l^2-m^2)*Vpdp)*exp(1j*(dot(k,r)));
end


x2z = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2z = x2z + (sqrt(3)/2*n*(l^2-m^2)*Vpds-n*(l^2-m^2)*Vpdp)*exp(1j*(dot(k,r)));
end


xyx = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xyx = xyx + (sqrt(3)*l^2*m*Vpds + m*(1-2*l^2)*Vpdp)*exp(1j*(dot(k,r)));
end


xyy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xyy = xyy + (sqrt(3)*m^2*l*Vpds + l*(1-2*m^2)*Vpdp)*exp(1j*(dot(k,r)));
end


xyz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xyz = xyz + (sqrt(3)*l*m*n*Vpds - 2*l*m*n*Vpdp)*exp(1j*(dot(k,r)));
end


xzx = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzx = xzx + (sqrt(3)*l^2*n*Vpds + n*(1-2*l^2)*Vpdp)*exp(1j*(dot(k,r)));
end


xzy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzy = xzy + (sqrt(3)*l*m*n*Vpds - 2*l*m*n*Vpdp)*exp(1j*(dot(k,r)));
end


xzz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzz = xzz + (sqrt(3)*n^2*l*Vpds+l*(1-2*n^2)*Vpdp)*exp(1j*(dot(k,r)));
end

yzx = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    yzx = yzx + (sqrt(3)*l*m*n*Vpds - 2*l*m*n*Vpdp)*exp(1j*(dot(k,r)));
end


yzy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    yzy = yzy + (sqrt(3)*m^2*n*Vpds + n*(1-2*m^2)*Vpdp)*exp(1j*(dot(k,r)));
end


yzz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    yzz = yzz + (sqrt(3)*n^2*m*Vpds+m*(1-2*n^2)*Vpdp)*exp(1j*(dot(k,r)));
end


z2s = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    z2s = z2s + (n^2-(l^2+m^2)/2)*Vsd*exp(1j*(dot(k,r)));
end


x2s = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    x2s = x2s + (sqrt(3)/2*(l^2-m^2)*Vsd)*exp(1j*(dot(k,r)));
end


xys = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xys = xys + (sqrt(3)*l*m*Vsd)*exp(1j*(dot(k,r)));
end


xzs = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    xzs = xzs + (sqrt(3)*l*n*Vsd)*exp(1j*(dot(k,r)));
end


yzs = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    yzs = yzs + (sqrt(3)*n*m*Vsd)*exp(1j*(dot(k,r)));
end


sAsB = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    sAsB = sAsB + Vss*exp(1j*(dot(k,r)));
end


sAx = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    sAx = sAx + l*Vsp*exp(1j*(dot(k,r)));
end


sAy = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    sAy = sAy + m*Vsp*exp(1j*(dot(k,r)));
end


sAz = 0;
for j = 1:nat
    
    l = dcs(j,1);
    m = dcs(j,2);
    n = dcs(j,3);
    
    k = [kx,ky,kz];
    r = a*[l,m,n];

    sAz = sAz + n*Vsp*exp(1j*(dot(k,r)));
end

%%

Hspsd = [ sAsB,sAx,sAy,sAz;
        z2s,z2x,z2y,z2z;
        x2s,x2x,x2y,x2z;
        xys,xyx,xyy,xyz;
        xzs,xzx,xzy,xzz;
        yzs,yzx,yzy,yzz];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%