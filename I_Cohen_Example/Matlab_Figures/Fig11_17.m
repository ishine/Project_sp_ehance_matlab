close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

P=20;
Lh=30;
beta=0.1; % sigma_w^2/sigma_u^2
alpha=0.8;
theta_d=0/180*pi;
theta_n=90/180*pi;
c=340; % voice velocity
fs=8e3;
delta=0.03;

L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);

M_values=[4 6 10 15]; % Number of sensors
[t1,t2]=meshgrid(1:L,1:L);
Rx=alpha.^(abs(t1-t2));

iSNR_db_values=-5:2:15;  
Gain_dB_values=zeros(length(iSNR_db_values),length(M_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(M_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(M_values));
nu_dB_values=zeros(length(iSNR_db_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    
    Rvw=beta*eye(Lh*M);
    Rvi=kron(ones(M),eye(Lh));
    Rv0=Rvi+Rvw;
    
    G_theta_d=zeros(M*Lh,L);
    for m=1:M
        Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
        G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
    end
    diag_GG=diag(G_theta_d'*G_theta_d);
    idx_max=find(diag_GG==max(diag_GG),1);
    i_ell=zeros(L,1);
    i_ell(idx_max)=1;
    

    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        sigma_i2=Rx(1,1)/Rv0(1,1)/iSNR;
        Rv_bar=sigma_i2*Rv0;
        
        h=pinv(G_theta_d*Rx*G_theta_d'+Rv_bar)*G_theta_d*Rx*i_ell;
        
        oSNR=(h'*G_theta_d*Rx*G_theta_d'*h)/(h'*Rv_bar*h);
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxM)=oSNR_dB-iSNR_dB;
        
        xi_n=Rv_bar(1,1)/(h'*Rv_bar*h);
        xi_n_dB_values(idxS,idxM)=10*log10(xi_n);
        
        xi_d=Rx(1,1)/(h'*G_theta_d*Rx*G_theta_d'*h);
        xi_d_dB_values(idxS,idxM)=10*log10(xi_d);
        
        nu=(G_theta_d'*h-i_ell)'*Rx*(G_theta_d'*h-i_ell)/Rx(1,1);
        nu_dB_values(idxS,idxM)=10*log10(nu);
    end
end

figure
plot(iSNR_db_values,Gain_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values,Gain_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,Gain_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,Gain_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values,xi_n_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values,xi_n_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,xi_n_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,xi_n_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values,xi_d_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values,xi_d_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,xi_d_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,xi_d_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values,nu_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values,nu_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,nu_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values,nu_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

