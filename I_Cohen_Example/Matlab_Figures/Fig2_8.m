close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A=0.5;
f_vec=0.1:0.03:0.45;
L=20; 
K=4;

hidx_values=1:4;  % index for different filters
iSNR_db_values=(-10:10)';  
oSNR_dB_values=zeros(length(iSNR_db_values),length(hidx_values));
MSE_values=zeros(length(iSNR_db_values),length(hidx_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(hidx_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(hidx_values));

i_i=[1; zeros(L-1,1)];
[t1,t2]=meshgrid(1:L,1:L);
Rx=0.5*A^2*cos(2*pi*f_vec(1)*(t1-t2));
for k=2:K
    Rx=Rx+0.5*A^2*cos(2*pi*f_vec(k)*(t1-t2));
end

for idx_snr=1:length(iSNR_db_values)
    iSNR_dB=iSNR_db_values(idx_snr);
    iSNR=10^(iSNR_dB/10);
    
    sigma_v2=K*A^2/2/iSNR;
    Rv=sigma_v2*eye(L);
    Ry=Rx+Rv;
    
    % h_max
    idx1=1;
    [Tx,~]=eigs(Rv\Rx,1);
    zeta=(Tx'*Rx*i_i)/(Tx'*Ry*Tx);
    h=zeta*Tx;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
    
    % h_Wiener
    idx1=2;
    h=Ry\Rx*i_i;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
    
    % h_MVDR
    idx1=3;
    [Tx,Xi_x]=eigs(Rx,2*K);
    h=Tx*Tx'*i_i;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
    
    % h_T
    idx1=4;
    mu=0.5;
    h=(Rx+mu*Rv)\Rx*i_i;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
end

figure
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,1)-iSNR_db_values(1:2:end),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,2)-iSNR_db_values(1:2:end),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,3)-iSNR_db_values(1:2:end),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,4)-iSNR_db_values(1:2:end),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;
