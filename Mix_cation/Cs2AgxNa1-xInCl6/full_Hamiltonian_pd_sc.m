function H = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB,selectcompd)
% Import kpath, TB parameters, and position of atoms from
% band_plot.m file. To run this code use:
% H = full_Hamiltonian_pd_sc(kx,ky,kz,t1,t2,t3,tNa1,tNaIn,tNaAg,posAA,posAB,sA,sB).
% As an output it gives a 80X80 Hamiltonian (with SOC included)
% of Bp-Xp orbitals-based HDP.
%% Desigh the SK-TB Hamiltonian
%% Ag
EAA = zeros(12);
EAA(1,1) = t1(7);EAA(7,7) = t1(7);
EAA(2,2) = t1(1);EAA(8,8) = t1(1);EAA(3,3) = t1(1);EAA(9,9) = t1(1);
EAA(4,4) = t1(2);EAA(5,5) = t1(2);EAA(6,6) = t1(2);EAA(10,10) = t1(2);
EAA(11,11) = t1(2);EAA(12,12) = t1(2);

t1(1)=0; t1(2)=0; t1(7)=0;
%% In
EBB = zeros(8);
EBB(1,1)=t2(6);EBB(5,5)=t2(6);
EBB(2,2)=t2(5);EBB(3,3)=t2(5);EBB(4,4)=t2(5);EBB(6,6)=t2(5);EBB(7,7)=t2(5);EBB(8,8)=t2(5);

t2(5)=0; t2(6)=0;

%% Na
ENa = zeros(8);
ENa(1,1)=tNa1(6);ENa(5,5)=tNa1(6);
ENa(2,2)=tNa1(5);ENa(3,3)=tNa1(5);ENa(4,4)=tNa1(5);ENa(6,6)=tNa1(5);ENa(7,7)=tNa1(5);ENa(8,8)=tNa1(5);

tNa1(5)=0; tNa1(6)=0;

%% Ag

HAA12  = Hamiltonian_d(kx,ky,kz,t1,posAA(1:4,:),1/sqrt(2));
HAA13  = Hamiltonian_d(kx,ky,kz,t1,posAA(5:8,:),1/sqrt(2));
HAA14  = Hamiltonian_d(kx,ky,kz,t1,posAA(9:12,:),1/sqrt(2));

%% In

HBB12  = Hamiltonian_p(kx,ky,kz,t2,posAA(1:4,:),1/sqrt(2));
HBB13  = Hamiltonian_p(kx,ky,kz,t2,posAA(5:8,:),1/sqrt(2));
HBB14  = Hamiltonian_p(kx,ky,kz,t2,posAA(9:12,:),1/sqrt(2));

%% AB
HAB12  = Hamiltonian_spsd(kx,ky,kz,t3,posAB(1:2,:),1/2);
HAB13  = Hamiltonian_spsd(kx,ky,kz,t3,posAB(3:4,:),1/2);
HAB14  = Hamiltonian_spsd(kx,ky,kz,t3,posAB(5:6,:),1/2);

[Hdsoc,Hpsoc] = HSOC_pd(sA,sB);

%% Na

HNa12  = Hamiltonian_p(kx,ky,kz,tNa1,posAA(1:4,:),1/sqrt(2));
HNa13  = Hamiltonian_p(kx,ky,kz,tNa1,posAA(5:8,:),1/sqrt(2));
HNa14  = Hamiltonian_p(kx,ky,kz,tNa1,posAA(9:12,:),1/sqrt(2));

HNaIn12  = Hamiltonian_p(kx,ky,kz,tNaIn,posAB(1:2,:),1/2);
HNaIn13  = Hamiltonian_p(kx,ky,kz,tNaIn,posAB(3:4,:),1/2);
HNaIn14  = Hamiltonian_p(kx,ky,kz,tNaIn,posAB(5:6,:),1/2);

HNaAg12  = Hamiltonian_spsd(kx,ky,kz,tNaAg,posAA(1:4,:),1/sqrt(2));
HNaAg13  = Hamiltonian_spsd(kx,ky,kz,tNaAg,posAA(5:8,:),1/sqrt(2));
HNaAg14  = Hamiltonian_spsd(kx,ky,kz,tNaAg,posAA(9:12,:),1/sqrt(2));
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

Na12 =  [HNa12,zeros(4,4);zeros(4,4),HNa12];
Na13 =  [HNa13,zeros(4,4);zeros(4,4),HNa13];
Na14 =  [HNa14,zeros(4,4);zeros(4,4),HNa14];

NaIn12 =  [HNaIn12,zeros(4,4);zeros(4,4),HNaIn12];
NaIn13 =  [HNaIn13,zeros(4,4);zeros(4,4),HNaIn13];
NaIn14 =  [HNaIn14,zeros(4,4);zeros(4,4),HNaIn14];

NaAg12 =  [HNaAg12,zeros(6,4);zeros(6,4),HNaAg12];
NaAg13 =  [HNaAg13,zeros(6,4);zeros(6,4),HNaAg13];
NaAg14 =  [HNaAg14,zeros(6,4);zeros(6,4),HNaAg14];
%% unit cell
%     H  = [EAA+AA12+AA13+AA14,AB12+AB13+AB14;
%           AB12'+AB13'+AB14',EBB+BB12+BB13+BB14];
%% Hamiltonian for mixed cation case
if selectcompd == 1
    %% CsAg0.75Na0.25InCl
    H = [Hdsoc+EAA,      AA13,            AA14,      NaAg12,        AB13,           AB12,  zeros(12,8),       AB14;
        AA13',      Hdsoc+EAA,            AA12,      NaAg13,        AB14,    zeros(12,8),         AB12,       AB13;
        AA14',          AA12',       Hdsoc+EAA,      NaAg14,   zeros(12,8),         AB14,         AB13,       AB12;
        NaAg12',      NaAg13',          NaAg14',       ENa,       NaIn12,         NaIn13,         NaIn14,     zeros(8);
        AB13',           AB14',      zeros(8,12),      NaIn12',       Hpsoc+EBB,    BB12,         BB14,       BB13;
        AB12',            zeros(8,12),   AB14',        NaIn13',       BB12',    Hpsoc+EBB,        BB13,       BB14;
        zeros(8,12),     AB12',          AB13',        NaIn14',        BB14',       BB13',   Hpsoc+EBB,       BB12;
        AB14',           AB13',          AB12',       zeros(8),       BB13',        BB14',       BB12',    Hpsoc+EBB];
end
if selectcompd == 2
    %% CsAg0.5Na0.5InCl
    H = [Hdsoc+EAA,      AA13,            NaAg14,      NaAg12,        AB13,           AB12,  zeros(12,8),       AB14;
        AA13',      Hdsoc+EAA,            NaAg12,      NaAg13,        AB14,    zeros(12,8),         AB12,       AB13;
        NaAg14',      NaAg12',             ENa,         Na14,      zeros(8),       NaIn14,       NaIn13,     NaIn12;
        NaAg12',      NaAg13',             Na14',        ENa,        NaIn12,         NaIn13,         NaIn14,     zeros(8);
        AB13',           AB14',        zeros(8),      NaIn12',       Hpsoc+EBB,    BB12,         BB14,       BB13;
        AB12',           zeros(8,12),   NaIn14',        NaIn13',       BB12',    Hpsoc+EBB,        BB13,       BB14;
        zeros(8,12),     AB12',          NaIn13',        NaIn14',        BB14',       BB13',   Hpsoc+EBB,       BB12;
        AB14',           AB13',          NaIn12',       zeros(8),       BB13',        BB14',       BB12',    Hpsoc+EBB];
end
if selectcompd == 3
    %% CsAg0.25Na0.75InCl
    H = [ENa,            Na13,              Na14,       NaAg12',          NaIn13,         NaIn12,       zeros(8),       NaIn14;
        Na13',            ENa,              Na12,       NaAg13',          NaIn14,       zeros(8),         NaIn12,       NaIn13;
        Na14',           Na12',              ENa,       NaAg14',        zeros(8),         NaIn14,         NaIn13,       NaIn12;
        NaAg12,         NaAg13,           NaAg14,    Hdsoc+EAA,           AB12,           AB13,           AB14,     zeros(12,8);
        NaIn13',        NaIn14',         zeros(8),      AB12',       Hpsoc+EBB,           BB12,           BB14,       BB13;
        NaIn12',         zeros(8),        NaIn14',      AB13',          BB12',       Hpsoc+EBB,           BB13,       BB14;
        zeros(8),        NaIn12',         NaIn13',      AB14',          BB14',           BB13',      Hpsoc+EBB,       BB12;
        NaIn14',         NaIn13',         NaIn12',     zeros(8,12),        BB13',           BB14',          BB12',    Hpsoc+EBB];
end