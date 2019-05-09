close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

iSNR_dB=(-10:20)';
iSNR=10.^(iSNR_dB/10);

mu=[0.5 1 2 5]';
[A_mu,A_iSNR]=meshgrid(mu,iSNR);

H_W=iSNR./(1+iSNR);
teta=acos(H_W.^0.5);

H=H_W*ones(1,3);
H(:,1)=1-sin(teta); % mag
H(:,3)=H_W.^0.5; % pow

xi_n=10*log10(1./(H.^2));
nu_d=10*log10((1-H).^2);

figure
plot(iSNR_dB(1:2:end),H(1:2:end,1),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_dB(1:2:end),H(1:2:end,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),H(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
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
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_dB(1:2:end),nu_d(1:2:end,1),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_dB(1:2:end),nu_d(1:2:end,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_dB(1:2:end),nu_d(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;
