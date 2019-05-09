close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

iSNR_dB=(-10:20)';
iSNR=10.^(iSNR_dB/10);

H_W=iSNR./(1+iSNR);
teta=acos(H_W.^0.5);
xi_n=10*log10(1./(H_W.^2));
nu_d=10*log10((1-H_W).^2);

figure
plot(iSNR_dB,H_W,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_dB,teta,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_dB,xi_n,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[0 30]);
box on; grid on;

figure
plot(iSNR_dB,nu_d,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;
