close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=cos(phi_rad).^2;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);


% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=0.25+0.5*cos(phi_rad)+0.25*cos(phi_rad).^2;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=-1/6+1/3*cos(phi_rad)+5/6*cos(phi_rad).^2;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
B=0.5/(3+7^0.5)+7^0.5/(3+7^0.5)*cos(phi_rad)+2.5/(3+7^0.5)*cos(phi_rad).^2;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

