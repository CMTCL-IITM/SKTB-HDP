function H = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB,selectcompd)
% Import kpath, TB parameters, and position of atoms from
% band_plot.m file. To run this code use:
% H = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,posAA,posAB,sA,sB).
% As an output it gives a 80X80 Hamiltonian (with SOC included)
% of d-p orbitals-based HDP.
%% Desigh the SK-TB Hamiltonian
for p = 1:2
    
    EAA = zeros(12);
    EAA(1,1) = t1(p,7);EAA(7,7) = t1(p,7);
    EAA(2,2) = t1(p,1);EAA(8,8) = t1(p,1);EAA(3,3) = t1(p,1);EAA(9,9) = t1(p,1);
    EAA(4,4) = t1(p,2);EAA(5,5) = t1(p,2);EAA(6,6) = t1(p,2);EAA(10,10) = t1(p,2);
    EAA(11,11) = t1(p,2);EAA(12,12) = t1(p,2);
    
    t1(p,1)=0; t1(p,2)=0; t1(p,7)=0;
    
    EBB = zeros(8);
    EBB(1,1)=t2(p,6);EBB(5,5)=t2(p,6);
    EBB(2,2)=t2(p,5);EBB(3,3)=t2(p,5);EBB(4,4)=t2(p,5);EBB(6,6)=t2(p,5);EBB(7,7)=t2(p,5);EBB(8,8)=t2(p,5);
    
    t2(p,5)=0; t2(p,6)=0;
    
    %%
    
    HAA12  = Hamiltonian_d(kx,ky,kz,t1(p,:),posAA(1:4,:),1/sqrt(2));
    HAA13  = Hamiltonian_d(kx,ky,kz,t1(p,:),posAA(5:8,:),1/sqrt(2));
    HAA14  = Hamiltonian_d(kx,ky,kz,t1(p,:),posAA(9:12,:),1/sqrt(2));
    
    HBB12  = Hamiltonian_p(kx,ky,kz,t2(p,:),posAA(1:4,:),1/sqrt(2));
    HBB13  = Hamiltonian_p(kx,ky,kz,t2(p,:),posAA(5:8,:),1/sqrt(2));
    HBB14  = Hamiltonian_p(kx,ky,kz,t2(p,:),posAA(9:12,:),1/sqrt(2));
    
    HAB12  = Hamiltonian_spsd(kx,ky,kz,t3(p,:),posAB(1:2,:),1/2);
    HAB13  = Hamiltonian_spsd(kx,ky,kz,t3(p,:),posAB(3:4,:),1/2);
    HAB14  = Hamiltonian_spsd(kx,ky,kz,t3(p,:),posAB(5:6,:),1/2);
    
    [Hdsoc,Hpsoc] = HSOC_pd(sA,sB);
    
    %%
    
    AB12 = [HAB12,zeros(6,4);zeros(6,4),HAB12];
    AB13 = [HAB13,zeros(6,4);zeros(6,4),HAB13];
    AB14 = [HAB14,zeros(6,4);zeros(6,4),HAB14];
    
    AA12 =  [HAA12,zeros(6,6);zeros(6,6),HAA12];
    AA13 =  [HAA13,zeros(6,6);zeros(6,6),HAA13];
    AA14 =  [HAA14,zeros(6,6);zeros(6,6),HAA14];
    
    BB12 =  [HBB12,zeros(4,4);zeros(4,4),HBB12];
    BB13 =  [HBB13,zeros(4,4);zeros(4,4),HBB13];
    BB14 =  [HBB14,zeros(4,4);zeros(4,4),HBB14];
    
    %% unit cell
    %     H  = [EAA+AA12+AA13+AA14,AB12+AB13+AB14;
    %           AB12'+AB13'+AB14',EBB+BB12+BB13+BB14];
    %% supercell
    H1(:,:,p) = [Hdsoc+EAA,       AA13,            AA14,      AA12,        AB13,           AB12,  zeros(12,8),       AB14;
        AA13',      Hdsoc+EAA,            AA12,      AA13,        AB14,    zeros(12,8),         AB12,       AB13;
        AA14',          AA12',       Hdsoc+EAA,      AA14,   zeros(12,8),         AB14,         AB13,       AB12;
        AA12',           AA13',          AA14',    Hdsoc+EAA,     AB12,           AB13,         AB14,     zeros(12,8);
        AB13',           AB14',      zeros(8,12),      AB12',       Hpsoc+EBB,    BB12,         BB14,       BB13;
        AB12',            zeros(8,12),   AB14',        AB13',       BB12',    Hpsoc+EBB,        BB13,       BB14;
        zeros(8,12),     AB12',          AB13',       AB14',        BB14',        BB13',   Hpsoc+EBB,       BB12;
        AB14',           AB13',          AB12',   zeros(8,12),      BB13',        BB14',       BB12',    Hpsoc+EBB];
end
%% setup for mixed cation

%% setup for In0.75Bi0.25
if selectcompd == 3
    Hc1 = H1(:,:,1);
    Hc2 = H1(:,:,2);
    
    H = zeros(80);
    
    H(1:72,1:72) = Hc1(1:72,1:72);
    H(:,73:80) = Hc2(:,73:80);
    H(73:80,:) = Hc2(73:80,:);
    H(73:80,73:80) = H1(73:80,73:80,2);
end
%% setup for In0.5Bi0.5
if selectcompd == 2
    
    Hc1 = H1(:,:,1);
    Hc2 = H1(:,:,2);
    
    H = zeros(80);
    
    H(1:64,1:64) = Hc2(1:64,1:64);
    H(1:64,65:80) = Hc1(1:64,65:80);
    H(65:80,1:64) = Hc1(65:80,1:64);
    H(65:80,65:80) = H1(65:80,65:80,1);
end
%% setup for In0.25Bi0.75
if selectcompd == 1
    
    Hc1 = H1(:,:,1);
    Hc2 = H1(:,:,2);
    
    H = zeros(80);
    
    H(1:72,1:72) = Hc2(1:72,1:72);
    H(:,73:80) = Hc1(:,73:80);
    H(73:80,:) = Hc1(73:80,:);
    H(73:80,73:80) = H1(73:80,73:80,1);
end