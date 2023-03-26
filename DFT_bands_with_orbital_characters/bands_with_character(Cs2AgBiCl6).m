%% smooth plotting of the orbital characters over the band structure for (Cs2AgBiCl6) 
% author: Mayank Gupta and B. R. K. Nanda
% Date: 09/03/2023
% doi:
% contact: nandab@iitm.ac.in and mayank77338@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
%% load the figure
figure(1)
hold on
tic
%% Input DFT parameters
nkpt = 202;
nbnd =  294;
bstart = 0;
%% load DFT obtained orbital characters
Bis = load('Bi-s.dat');
Bip = load('Bi-p.dat');
Ags = load('Ag-s.dat');
Ageg = load('Ag-eg.dat');
Agt2g = load('Ag-t2g.dat');
Clp = load('Cl-p.dat');

a1 = [Bis(:,3);Bip(:,3);Ags(:,3);Ageg(:,3);Agt2g(:,3);Clp(:,3)];
a2 = rescale(a1);
a3 = reshape(a2,[length(Bis),6]);

Bis(:,3)    = a3(:,1);
Bip(:,3)    = a3(:,2);
Ags(:,3)    = a3(:,3);
Ageg(:,3)   = a3(:,4);
Agt2g(:,3)  = a3(:,5);
Clp(:,3)    = a3(:,6);

%% plot band without character

eigval = reshape(Bis(:,2),[nkpt,nbnd]);
plot(Bis(1:nkpt),eigval,'color',[0.5 0.5 0.5],'LineWidth',0.5)

%% plot orbital character
orb = {Bis;Bip;Ags;Ageg;Agt2g;Clp};
color = {'g','c',[0.49 0.18 0.55],'r','b','y'};

for M = [2,3,4,5,1,6]
    
    k=bstart*nkpt+1;
    abc = orb{M};
    
    for j=1:nbnd
        
        yy1=abc(k:k+nkpt-1,1);
        yy2=abc(k:k+nkpt-1,2);
        yy3=abc(k:k+nkpt-1,3)*2;
        
        yy2(end) = NaN;
        
        patch(yy1,yy2,'r','EdgeColor',color{M},...
            'FaceVertexAlphaData',yy3,'AlphaDataMapping','none',...
            'EdgeAlpha','interp','LineWidth',5)
        k=k+nkpt;
    end
end

%%  Figure smoothing and cleaning
ax=gca;

ax.Box = 'off';
ax.LineWidth = 0.005;
ax.FontSize = 0.001;
ax.TickDir = 'in';
ax.TickLength = [0.001 0.001];
ax.YLim = [-10 10];
ax.XLim = [0 1.07596];
% ax.YTick = [-12 -8 -4 0 4 8];
ax.XColor = 'none';
ax.YColor = 'none';

toc
