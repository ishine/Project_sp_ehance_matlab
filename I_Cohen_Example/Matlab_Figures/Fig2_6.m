close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

alpha=0.8;
L=30;

xin_dB_values=[12 13 14 15];

iSNR_db_values=(0:20)';  
oSNR_dB_values=zeros(length(iSNR_db_values),length(xin_dB_values));
MSE_values=zeros(length(iSNR_db_values),length(xin_dB_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(xin_dB_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(xin_dB_values));
for idx_xi=1:length(xin_dB_values)
    xin_dB=xin_dB_values(idx_xi);
    xin=10^(xin_dB/10);
    fun=@(x) myfunc2(x,xin);
    mu_sig2=fzero(fun,0.1);
    
    i_i=[1; zeros(L-1,1)];
    [t1,t2]=meshgrid(1:L,1:L);
    Rx=alpha.^(abs(t1-t2));
    for idx_snr=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idx_snr);
        iSNR=10^(iSNR_dB/10);
        
        sigma_v2=1/iSNR;
        Rv=sigma_v2*eye(L);
        Ry=Rx+Rv;
        
        mu=mu_sig2/sigma_v2;
        h=(Rx+mu*Rv)\Rx*i_i;
        
        oSNR=(h'*Rx*h)/(h'*Rv*h);
        oSNR_dB=10*log10(oSNR);
        oSNR_dB_values(idx_snr,idx_xi)=oSNR_dB;
        
        MSE_values(idx_snr,idx_xi)=10*log10(1-2*h'*Rx*i_i+h'*Ry*h);
        
        xi_n=sigma_v2/(h'*Rv*h);
        xi_n_dB_values(idx_snr,idx_xi)=10*log10(xi_n);
        
        xi_d=1/(h'*Rx*h);
        xi_d_dB_values(idx_snr,idx_xi)=10*log10(xi_d);
    end
end

figure
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,1)-iSNR_db_values(1:2:end),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,2)-iSNR_db_values(1:2:end),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,3)-iSNR_db_values(1:2:end),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),oSNR_dB_values(1:2:end,4)-iSNR_db_values(1:2:end),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[4.8 5.8]);
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),MSE_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
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
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[11 16]);
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[6 10]);
box on; grid on;

