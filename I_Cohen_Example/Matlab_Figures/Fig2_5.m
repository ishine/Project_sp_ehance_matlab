close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=12;

A=0.5;
f0=0.1;
L=30;

iSNR_dB=0;
iSNR=10^(iSNR_dB/10);
sigma_v2=A^2/2/iSNR;

t=(0:500)';
x=A*cos(2*pi*f0*t+pi/3);
y=x+sigma_v2^0.5*randn(size(t));

Rv=sigma_v2*eye(L);
[t1,t2]=meshgrid(1:L,1:L);
Rx=0.5*A^2*cos(2*pi*f0*(t1-t2));
Ry=Rx+Rv;
i_i=[1; zeros(L-1,1)];
hW=Ry\Rx*i_i;
z=conv(y,hW,'same');

figure
plot(t,y,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-2 2]);
box on; grid on;

figure
plot(t,z,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-1 1]);
box on; grid on;

