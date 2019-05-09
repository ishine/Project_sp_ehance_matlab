close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

broadband_iSNR_dB=(-5:15)';
broadband_iSNR=10.^(broadband_iSNR_dB/10);

f=linspace(0,0.5,1e4);  % frequency
f=f(1:end-1);
df=f(2)-f(1);

phiX=zeros(size(f));
phiX(f<=0.25)=1;

narrowband_iSNR=(broadband_iSNR * (phiX./(1-2*f)));
H_W=narrowband_iSNR./(1+narrowband_iSNR);

int_phiX=df*ones(size(broadband_iSNR))*sum(phiX);
int_H2_phiX=df*sum(H_W.^2 .* (ones(size(broadband_iSNR))*phiX),2);
int_phiV=(1./broadband_iSNR) .* df .* (ones(size(broadband_iSNR))*sum(1-2*f));
int_H2_phiV=df* sum(H_W.^2 .* ((1./broadband_iSNR)*(1-2*f)),2);
int_Hm12_phiX=df*sum(abs(H_W-1).^2 .* (ones(size(broadband_iSNR))*phiX),2);

oSNR=int_H2_phiX./int_H2_phiV;
oSNR_dB=10*log10(oSNR);
gain=oSNR_dB-broadband_iSNR_dB;
MSE=10*log10(int_H2_phiV+int_Hm12_phiX);
xi_n=10*log10(int_phiV./int_H2_phiV);
xi_d=10*log10(int_phiX./int_H2_phiX);

figure
plot(broadband_iSNR_dB,gain,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(broadband_iSNR_dB,MSE,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(broadband_iSNR_dB,xi_n,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(broadband_iSNR_dB,xi_d,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;
