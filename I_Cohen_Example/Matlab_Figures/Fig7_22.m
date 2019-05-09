close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity
theta_d=0/180*pi;
theta1=45/180*pi;
theta2=90/180*pi;
delta=1e-2;
alpha=1e-3;
M=8;     % number of microphones
f_vec=[1 2 4 8]*1e3;  % frequencies
i_c=[1 0 0]';

idx_f=1;
f=f_vec(idx_f);
[m_mat,n_mat]=meshgrid(1:M,1:M);
Gamma0=sinc(2*f*delta/c*(m_mat-n_mat));
Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
C=[exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta1)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta2))];
h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;

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

idx_f=2;
f=f_vec(idx_f);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
[m_mat,n_mat]=meshgrid(1:M,1:M);
Gamma0=sinc(2*f*delta/c*(m_mat-n_mat));
Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
C=[exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta1)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta2))];
h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;

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

idx_f=3;
f=f_vec(idx_f);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
[m_mat,n_mat]=meshgrid(1:M,1:M);
Gamma0=sinc(2*f*delta/c*(m_mat-n_mat));
Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
C=[exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta1)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta2))];
h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;

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

idx_f=4;
f=f_vec(idx_f);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
[m_mat,n_mat]=meshgrid(1:M,1:M);
Gamma0=sinc(2*f*delta/c*(m_mat-n_mat));
Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
C=[exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta1)) exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta2))];
h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;

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
