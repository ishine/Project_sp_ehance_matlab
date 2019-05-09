close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity

delta=0.005;     % distance between microphones
f=1e3;
theta_d=0/180*pi;

M=6;     % number of microphones
[m_mat,n_mat]=meshgrid(1:M,1:M);
GammaC=besseli(0,2i*pi*f*(m_mat-n_mat)*delta/c);
b1=0.5*[3^0.5-1; 3-3^0.5];  % supercardioid
f_bar=2*pi*delta/c*f*(0:M-1)';
b0bar=besselj(0,-f_bar);
b1bar=1i*besselj(1,-f_bar); 
Gamma_dpc=[b0bar b1bar];
d=exp(-1i*f_bar);

epsilon=0;
GammaC_eps=GammaC+epsilon*eye(M);
hLS=GammaC_eps\Gamma_dpc*b1;
h=hLS-(1-d'*hLS)/(d'/GammaC_eps*d)*(GammaC_eps\d);

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

epsilon=1e-5;
GammaC_eps=GammaC+epsilon*eye(M);
hLS=GammaC_eps\Gamma_dpc*b1;
h=hLS-(1-d'*hLS)/(d'/GammaC_eps*d)*(GammaC_eps\d);

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

epsilon=1e-3;
GammaC_eps=GammaC+epsilon*eye(M);
hLS=GammaC_eps\Gamma_dpc*b1;
h=hLS-(1-d'*hLS)/(d'/GammaC_eps*d)*(GammaC_eps\d);

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

epsilon=0.1;
GammaC_eps=GammaC+epsilon*eye(M);
hLS=GammaC_eps\Gamma_dpc*b1;
h=hLS-(1-d'*hLS)/(d'/GammaC_eps*d)*(GammaC_eps\d);

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

