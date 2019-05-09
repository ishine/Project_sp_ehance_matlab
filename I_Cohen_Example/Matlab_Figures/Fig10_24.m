close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity

delta=0.005;     % distance between microphones
f=1e3;
theta_d=0/180*pi;

b1=0.5*[3^0.5-1; 3-3^0.5];  % supercardioid
f_bar=2*pi*delta/c*f*(0:1)';
b0bar=besselj(0,f_bar); 
b1bar=2i*besselj(1,f_bar); 
B1bar=[b0bar.'; b1bar.'];
hNR=B1bar\b1;

M=6;     % number of microphones

aleph=0;
H=toeplitz([conj(hNR(1)); zeros(M-2,1)],[hNR' zeros(1,M-2)])';
f_bar=2*pi*delta/c*f*(0:M-1)';
d=exp(-1i*f_bar);
dt=H'*d;
[m_mat,n_mat]=meshgrid(1:M,1:M);
GammaC=besseli(0,2i*pi*f*(m_mat-n_mat)*delta/c);
R=H'*GammaC*H;
b0bar=besselj(0,-f_bar);
b1bar=1i*besselj(1,-f_bar); 
Gamma_dpc=[b0bar b1bar];
Ra=aleph*R+(1-aleph)*H'*H;
g_Ua=aleph*(Ra\H')*Gamma_dpc*b1;
g_Ta=g_Ua+(1-dt'*g_Ua)/(dt'/Ra*dt)*(Ra\dt);
h=H*g_Ta;

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

aleph=0.5;
H=toeplitz([conj(hNR(1)); zeros(M-2,1)],[hNR' zeros(1,M-2)])';
f_bar=2*pi*delta/c*f*(0:M-1)';
d=exp(-1i*f_bar);
dt=H'*d;
[m_mat,n_mat]=meshgrid(1:M,1:M);
GammaC=besseli(0,2i*pi*f*(m_mat-n_mat)*delta/c);
R=H'*GammaC*H;
b0bar=besselj(0,-f_bar);
b1bar=1i*besselj(1,-f_bar); 
Gamma_dpc=[b0bar b1bar];
Ra=aleph*R+(1-aleph)*H'*H;
g_Ua=aleph*(Ra\H')*Gamma_dpc*b1;
g_Ta=g_Ua+(1-dt'*g_Ua)/(dt'/Ra*dt)*(Ra\dt);
h=H*g_Ta;

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

aleph=0.9;
H=toeplitz([conj(hNR(1)); zeros(M-2,1)],[hNR' zeros(1,M-2)])';
f_bar=2*pi*delta/c*f*(0:M-1)';
d=exp(-1i*f_bar);
dt=H'*d;
[m_mat,n_mat]=meshgrid(1:M,1:M);
GammaC=besseli(0,2i*pi*f*(m_mat-n_mat)*delta/c);
R=H'*GammaC*H;
b0bar=besselj(0,-f_bar);
b1bar=1i*besselj(1,-f_bar); 
Gamma_dpc=[b0bar b1bar];
Ra=aleph*R+(1-aleph)*H'*H;
g_Ua=aleph*(Ra\H')*Gamma_dpc*b1;
g_Ta=g_Ua+(1-dt'*g_Ua)/(dt'/Ra*dt)*(Ra\dt);
h=H*g_Ta;

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

aleph=0.99;
H=toeplitz([conj(hNR(1)); zeros(M-2,1)],[hNR' zeros(1,M-2)])';
f_bar=2*pi*delta/c*f*(0:M-1)';
d=exp(-1i*f_bar);
dt=H'*d;
[m_mat,n_mat]=meshgrid(1:M,1:M);
GammaC=besseli(0,2i*pi*f*(m_mat-n_mat)*delta/c);
R=H'*GammaC*H;
b0bar=besselj(0,-f_bar);
b1bar=1i*besselj(1,-f_bar); 
Gamma_dpc=[b0bar b1bar];
Ra=aleph*R+(1-aleph)*H'*H;
g_Ua=aleph*(Ra\H')*Gamma_dpc*b1;
g_Ta=g_Ua+(1-dt'*g_Ua)/(dt'/Ra*dt)*(Ra\dt);
h=H*g_Ta;

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

