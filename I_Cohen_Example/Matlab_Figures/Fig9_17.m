close all; clc; clear all;
linewd = 0.8;

% Input:
c=340; % voice velocity
theta_d=0/180*pi;
M=4;     % number of microphones

delta=0.01;     % distance between microphones
f=0.5e3;
[m_mat,n_mat]=meshgrid(1:M,1:M);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
psi1=0; psi2=pi/2;
Gamma_0_pi2=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
psi1=pi/2; psi2=pi;
Gamma_pi2_pi=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
[T,Lambda]=jeig(Gamma_0_pi2,Gamma_pi2_pi);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);
h=t1/(d'*t1);

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


delta=0.01;     % distance between microphones
f=7e3;
[m_mat,n_mat]=meshgrid(1:M,1:M);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
psi1=0; psi2=pi/2;
Gamma_0_pi2=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
psi1=pi/2; psi2=pi;
Gamma_pi2_pi=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
[T,Lambda]=jeig(Gamma_0_pi2,Gamma_pi2_pi);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);
h=t1/(d'*t1);

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

delta=0.04;     % distance between microphones
f=0.5e3;
[m_mat,n_mat]=meshgrid(1:M,1:M);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
psi1=0; psi2=pi/2;
Gamma_0_pi2=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
psi1=pi/2; psi2=pi;
Gamma_pi2_pi=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
[T,Lambda]=jeig(Gamma_0_pi2,Gamma_pi2_pi);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);
h=t1/(d'*t1);

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

delta=0.04;     % distance between microphones
f=7e3;
[m_mat,n_mat]=meshgrid(1:M,1:M);
d=exp(-1i*2*pi*f*delta/c*(0:M-1)'*cos(theta_d));
psi1=0; psi2=pi/2;
Gamma_0_pi2=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
psi1=pi/2; psi2=pi;
Gamma_pi2_pi=(exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*f*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*f*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
[T,Lambda]=jeig(Gamma_0_pi2,Gamma_pi2_pi);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);
h=t1/(d'*t1);

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

