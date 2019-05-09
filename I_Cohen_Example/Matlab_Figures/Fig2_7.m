close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A=0.5;
f_vec=[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4];
L=30;

K_values=[1 2 4 8];
iSNR_db_values=(0:20)';  
oSNR_dB_values=zeros(length(iSNR_db_values),length(K_values));
MSE_values=zeros(length(iSNR_db_values),length(K_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(K_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(K_values));
for idx1=1:length(K_values)
    K=K_values(idx1);
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
        
    end
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
set(gca,'YLim',[-1 1]);
box on; grid on;
