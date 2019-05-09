close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

M_values=[10 20 30 40]; % Number of sensors

% Input:
c=340; % voice velocity

A0=0.5;
T=500;
f0=0.1;
f=0.1;

theta_d=70/180*pi;
theta_u1=50/180*pi;
theta_u2=30/180*pi;

alpha=0.01; % sigma_w^2/sigma_u^2

M=M_values(1);

d=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_d));
du1=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u1));
du2=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u2));
phiV0=T*du1*du1'+T*du2*du2'+T*alpha*eye(M);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_dB=0;
iSNR=10^(iSNR_dB/10);

phiV1=T*A0^2/2./iSNR;
sigma2=T*A0^2/2./iSNR/phiV0(1,1);
phiV=sigma2*phiV0;

I_M=eye(M);
i_i=I_M(:,1);

phiY=phiV+phiX1*d*d';
i_c=[1 0 0]';
C=[d du1 du2];
h=inv(phiY)*C*inv(C'*inv(phiY)*C)*i_c;


% Beampattern
phi=(0:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba_dB=20*log10(Ba);

figure
plot(phi_rad/pi*180, Ba_dB); grid
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',[0 30 50 70 90:30:180]);
set(gca,'YLim',[-100 20]);
box on; grid on;

M=M_values(2);

d=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_d));
du1=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u1));
du2=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u2));
phiV0=T*du1*du1'+T*du2*du2'+T*alpha*eye(M);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_dB=0;
iSNR=10^(iSNR_dB/10);

phiV1=T*A0^2/2./iSNR;
sigma2=T*A0^2/2./iSNR/phiV0(1,1);
phiV=sigma2*phiV0;

I_M=eye(M);
i_i=I_M(:,1);

phiY=phiV+phiX1*d*d';
C=[d du1 du2];
h=inv(phiY)*C*inv(C'*inv(phiY)*C)*i_c;


% Beampattern
phi=(0:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba_dB=20*log10(Ba);

figure
plot(phi_rad/pi*180, Ba_dB); grid
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',[0 30 50 70 90:30:180]);
set(gca,'YLim',[-100 20]);
box on; grid on;

M=M_values(3);

d=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_d));
du1=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u1));
du2=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u2));
phiV0=T*du1*du1'+T*du2*du2'+T*alpha*eye(M);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_dB=0;
iSNR=10^(iSNR_dB/10);

phiV1=T*A0^2/2./iSNR;
sigma2=T*A0^2/2./iSNR/phiV0(1,1);
phiV=sigma2*phiV0;

I_M=eye(M);
i_i=I_M(:,1);

phiY=phiV+phiX1*d*d';
C=[d du1 du2];
h=inv(phiY)*C*inv(C'*inv(phiY)*C)*i_c;

% Beampattern
phi=(0:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba_dB=20*log10(Ba);

figure
plot(phi_rad/pi*180, Ba_dB); grid
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',[0 30 50 70 90:30:180]);
set(gca,'YLim',[-100 20]);
box on; grid on;

M=M_values(4);

d=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_d));
du1=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u1));
du2=exp(-1i*2*pi*f.*(0:(M-1))'*cos(theta_u2));
phiV0=T*du1*du1'+T*du2*du2'+T*alpha*eye(M);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_dB=0;
iSNR=10^(iSNR_dB/10);

phiV1=T*A0^2/2./iSNR;
sigma2=T*A0^2/2./iSNR/phiV0(1,1);
phiV=sigma2*phiV0;

I_M=eye(M);
i_i=I_M(:,1);

phiY=phiV+phiX1*d*d';
C=[d du1 du2];
h=inv(phiY)*C*inv(C'*inv(phiY)*C)*i_c;

% Beampattern
phi=(0:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*m_mat.*cos(phi_mat));
B=D_phi'*h;
Ba=abs(B);
Ba_dB=20*log10(Ba);

figure
plot(phi_rad/pi*180, Ba_dB); grid
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',[0 30 50 70 90:30:180]);
set(gca,'YLim',[-100 20]);
box on; grid on;

