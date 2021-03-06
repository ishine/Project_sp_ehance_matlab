close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=cos(phi_rad);
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=0.5+0.5*cos(phi_rad);
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=0.25+0.75*cos(phi_rad);
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=0.5*(3^0.5-1)+0.5*(3-3^0.5)*cos(phi_rad);
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

