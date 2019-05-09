close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

iSNR_dB=(-10:20)';
iSNR=10.^(iSNR_dB/10);

mu=[0.5 1 2 5]';
[A_mu,A_iSNR]=meshgrid(mu,iSNR);

H_T=A_iSNR./(A_mu+A_iSNR);
teta=acos(H_T.^0.5);
xi_n=10*log10(1./(H_T.^2));
nu_d=10*log10((1-H_T).^2);

figure
plot(iSNR_dB(1:2:end),H_T(1:2:end,1),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_dB(1:2:end),H_T(1:2:end,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),H_T(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),H_T(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_dB(1:2:end),teta(1:2:end,1),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_dB(1:2:end),teta(1:2:end,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),teta(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),teta(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_dB(1:2:end),xi_n(1:2:end,1),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_dB(1:2:end),xi_n(1:2:end,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),xi_n(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),xi_n(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[0 30]);
box on; grid on;

figure
plot(iSNR_dB(1:2:end),nu_d(1:2:end,1),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_dB(1:2:end),nu_d(1:2:end,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),nu_d(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),nu_d(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;
