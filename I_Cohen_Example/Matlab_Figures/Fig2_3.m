close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=12;

A=0.5;
f0=0.1;
L=30;

iSNR_db_values=0:20;  
oSNR_dB_values=zeros(size(iSNR_db_values));
MSE_values=zeros(size(iSNR_db_values));
for idx=1:length(iSNR_db_values)
    iSNR_dB=iSNR_db_values(idx);
    iSNR=10^(iSNR_dB/10);
    
    sigma_v2=A^2/2/iSNR;
    Rv=sigma_v2*eye(L);
    
    [t1,t2]=meshgrid(1:L,1:L);
    Rx=0.5*A^2*cos(2*pi*f0*(t1-t2));
    Ry=Rx+Rv;
    i_i=[1; zeros(L-1,1)];
    hW=Ry\Rx*i_i;
    
    oSNR=(hW'*Rx*hW)/(hW'*Rv*hW);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx)=oSNR_dB;
    
    MSE_values(idx)=10*log10(A^2/2-2*hW'*Rx*i_i+hW'*Ry*hW);
end

figure
plot(iSNR_db_values,oSNR_dB_values,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values,MSE_values,'linewidth',linewd, 'MarkerSize',MarkerSize);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

