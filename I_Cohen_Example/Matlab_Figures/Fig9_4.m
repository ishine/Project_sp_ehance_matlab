close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity

delta=0.01;     % distance between microphones
f=0.5e3;
theta_d=0/180*pi;
M=2;     % number of microphones
tau0=delta/c;

alpha11=0; % dipole
h=1/(1-exp(1i*2*pi*f*tau0*(1-alpha11)))*[1; -exp(-1i*2*pi*f*tau0*alpha11)];

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*delta/c*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

alpha11=-1; % cardioid
h=1/(1-exp(1i*2*pi*f*tau0*(1-alpha11)))*[1; -exp(-1i*2*pi*f*tau0*alpha11)];

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*delta/c*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

alpha11=-1/3; % hypercardioid
h=1/(1-exp(1i*2*pi*f*tau0*(1-alpha11)))*[1; -exp(-1i*2*pi*f*tau0*alpha11)];

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*delta/c*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

alpha11=(1-3^0.5)/(3-3^0.5); % supercardioid
h=1/(1-exp(1i*2*pi*f*tau0*(1-alpha11)))*[1; -exp(-1i*2*pi*f*tau0*alpha11)];

% Beampattern
phi=(-180:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*delta/c*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
set(gcf, 'Color', [1, 1, 1]);
db_polar_m(phi_rad, Ba, -50, 0, linewd);

