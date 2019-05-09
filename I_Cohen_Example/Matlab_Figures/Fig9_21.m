close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity
theta_d=0/180*pi;

delta=0.005;     % distance between microphones
M=4;     % number of microphones
f=0.5e3;

alpha31=-1;
S=diag(0:M-1);
D=[exp(-1i*2*pi*f*delta/c*(0:M-1)')'; exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31)'; (S*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'; (S^2*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'];
h=D'/(D*D')*[1 0 0 0]';

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

M=4;     % number of microphones
f=7e3;

alpha31=-1;
S=diag(0:M-1);
D=[exp(-1i*2*pi*f*delta/c*(0:M-1)')'; exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31)'; (S*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'; (S^2*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'];
h=D'/(D*D')*[1 0 0 0]';

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

M=8;     % number of microphones
f=0.5e3;

alpha31=-1;
S=diag(0:M-1);
D=[exp(-1i*2*pi*f*delta/c*(0:M-1)')'; exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31)'; (S*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'; (S^2*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'];
h=D'/(D*D')*[1 0 0 0]';

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

M=8;     % number of microphones
f=7e3;

alpha31=-1;
S=diag(0:M-1);
D=[exp(-1i*2*pi*f*delta/c*(0:M-1)')'; exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31)'; (S*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'; (S^2*exp(-1i*2*pi*f*delta/c*(0:M-1)'*alpha31))'];
h=D'/(D*D')*[1 0 0 0]';

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

