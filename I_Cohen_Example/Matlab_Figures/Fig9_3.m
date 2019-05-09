close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=cos(phi_rad).^3;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);


% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=0.125+3/8*cos(phi_rad)+3/8*cos(phi_rad).^2+1/8*cos(phi_rad).^3;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=-3/32-15/32*cos(phi_rad)+15/32*cos(phi_rad).^2+35/32*cos(phi_rad).^3;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
[jj,ii]=meshgrid(0:3,0:3);
H1=(-1).^(ii+jj)./(1+ii+jj);
H2=1./(1+ii+jj);
[a_max,~]=eigs(H1\H2,1);
c=a_max/sum(a_max);
B=c(1)+c(2)*cos(phi_rad)+c(3)*cos(phi_rad).^2+c(4)*cos(phi_rad).^3;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);
