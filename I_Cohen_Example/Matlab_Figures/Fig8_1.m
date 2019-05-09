close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
A0=0.5;
T=500;
f0=0.1;

theta_d=70/180*pi;
theta_u1=50/180*pi;
theta_u2=30/180*pi;

alpha=0.01; % sigma_w^2/sigma_u^2

M_values=[10 20 30 40]; % Number of sensors

f=linspace(-0.5,0.5,1e3)';  % frequency

f=f(2:end-1);
df=f(2)-f(1);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0))./sin(pi*(f-f0))).^2;

iSNR_db_values=-5:15;  
Gain_dB_values=zeros(length(iSNR_db_values),length(M_values));
nu_values=zeros(length(iSNR_db_values),length(M_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(M_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    [f_mat,m_mat]=meshgrid(f,1:M);
    d=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta_d));
    du1=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta_u1));
    du2=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta_u2));
    phiV0=zeros(M,M,length(f));
    for idxf=1:length(f)
        phiV0(:,:,idxf)=T*du1(:,idxf)*du1(:,idxf)'+T*du2(:,idxf)*du2(:,idxf)'+T*alpha*eye(M);
    end
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        phiV1=T*A0^2/2./iSNR*ones(size(f));
        sigma2=T*A0^2/2./iSNR/phiV0(1,1,1);
        phiV=sigma2*phiV0;
        
        h=zeros(M,length(f));
        for idxf=1:length(f)
            phiY=phiV(:,:,idxf)+phiX1(idxf)*d(:,idxf)*d(:,idxf)';
            h(:,idxf)=(I_M-inv(phiY)*phiV(:,:,idxf))*i_i;
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
        
        xi_n_dB_values(idxS,idxM)=10*log10(int_phiV1./int_H2_phiV);
        xi_d_dB_values(idxS,idxM)=10*log10(int_phiX1./int_H2_phiX);
        nu_values(idxS,idxM)=10*log10(int_Hm12_phiX./int_phiX1);
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

figure
plot(iSNR_db_values(1:2:end),nu_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),nu_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),nu_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),nu_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

