close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A=0.5;
f0=0.1;
L=30;

broadband_iSNR_dB=0;
broadband_iSNR=10.^(broadband_iSNR_dB/10);
sigma_v2=A^2/2/broadband_iSNR;

T=500;
t=(0:499)';
x=A*cos(2*pi*f0*t+pi/3);
y=x+sigma_v2^0.5*randn(size(t));

Y=fft(y,512);
f=[0:256 -255:-1]'/512;
df=f(2)-f(1);

phiX=A^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A^2/4*(sin(T*pi*(f-f0))./sin(pi*(f-f0))).^2;
phiV=T*A^2/2./broadband_iSNR;

narrowband_iSNR=(1./phiV) * phiX;
H_W=narrowband_iSNR./(1+narrowband_iSNR);

Z=H_W.*Y;
z=ifft(Z);

f1=fftshift(f);
Ya=fftshift(abs(Y));

figure
plot(f1(2:end),Ya(2:end),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',-0.5:0.1:0.5); 
box on; grid on;

Za=fftshift(abs(Z));

figure
plot(f1(2:end),Za(2:end),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',-0.5:0.1:0.5); 
box on; grid on;

locplot('fig3_4c');
plot(t,y,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-2 2]);
box on; grid on;

figure
plot(t,z(1:T),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-1 1]);
box on; grid on;


