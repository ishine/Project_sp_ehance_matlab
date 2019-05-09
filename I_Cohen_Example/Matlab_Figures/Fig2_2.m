close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=12;

A=0.5;
f0=0.1;
sigma_v2=0.3;

L_values=1:60;  % length of the filter
Gain_dB_values=zeros(size(L_values));
MSE_values=zeros(size(L_values));
for idx=1:length(L_values)
    L=L_values(idx);
    
    iSNR=A^2/2/sigma_v2;
    iSNR_dB=10*log10(iSNR);
    Rv=sigma_v2*eye(L);
    
    [t1,t2]=meshgrid(1:L,1:L);
    Rx=0.5*A^2*cos(2*pi*f0*(t1-t2));
    Ry=Rx+Rv;
    i_i=[1; zeros(L-1,1)];
    hW=Ry\Rx*i_i;
    
    oSNR=(hW'*Rx*hW)/(hW'*Rv*hW);
    oSNR_dB=10*log10(oSNR);
    Gain_dB_values(idx)=oSNR_dB-iSNR_dB;
    
    MSE_values(idx)=10*log10(A^2/2-2*hW'*Rx*i_i+hW'*Ry*hW);
end

figure
plot(L_values,Gain_dB_values,'.','linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;

figure
plot(L_values,MSE_values,'.','linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;
