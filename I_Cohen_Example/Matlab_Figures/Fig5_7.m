close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
A0=0.5;
T=500;
f0=0.1;
theta=70/180*pi;
alpha=0.01; 

M=5; % Number of sensors
mu_values=[0.5 1 2 5];

f=linspace(-0.5,0.5,1e4)';  % frequency
f=f(2:end-1);
df=f(2)-f(1);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0))./sin(pi*(f-f0))).^2;

iSNR_db_values=-5:15;  
Gain_dB_values=zeros(length(iSNR_db_values),length(mu_values));
MSE_values=zeros(length(iSNR_db_values),length(mu_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(mu_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(mu_values));
for idxM=1:length(mu_values)
    mu=mu_values(idxM);
    
    [f_mat,m_mat]=meshgrid(f,1:M);
    d=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta));
    du=exp(-1i*2*pi*f_mat.*(m_mat-1));
    phiV0=zeros(M,M,length(f));
    for idxf=1:length(f)
        phiV0(:,:,idxf)=T*du(:,idxf)*du(:,idxf)'+T*alpha*eye(M);
    end
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        phiV1=T*A0^2/2./iSNR*ones(size(f));
        sigma2=A0^2/2./iSNR/(1+alpha);
        phiV=sigma2*phiV0;
        
        h=zeros(M,length(f));
        for idxf=1:length(f)
            h(:,idxf)=phiX1(idxf)*inv(phiX1(idxf)*d(:,idxf)*d(:,idxf)'+mu* phiV(:,:,idxf))*d(:,idxf);
        end
        
        int_phiX1=df*sum(phiX1);
        int_phiV1=df*sum(phiV1);
        int_H2_phiX=df*sum(phiX1 .* abs(sum(conj(h).*d,1)').^2 );
        int_H2_phiV=0;
        for idxf=1:length(f)
            int_H2_phiV=int_H2_phiV+df* h(:,idxf)'*phiV(:,:,idxf)*h(:,idxf);
        end
        int_Hm12_phiX=df*sum( phiX1 .* abs(sum(conj(h).*d,1)'-1).^2);
        
        oSNR=int_H2_phiX./int_H2_phiV;
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxM)=oSNR_dB-iSNR_dB;
        
        MSE_values(idxS,idxM)=10*log10(int_Hm12_phiX./int_phiX1);
        
        xi_n_dB_values(idxS,idxM)=10*log10(int_phiV1./int_H2_phiV);
        xi_d_dB_values(idxS,idxM)=10*log10(int_phiX1./int_H2_phiX);
    end
end

figure
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
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
box on; grid on;

