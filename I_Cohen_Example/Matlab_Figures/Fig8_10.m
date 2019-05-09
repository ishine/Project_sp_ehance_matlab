close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
A0=0.5;
T=512;
f0=0.1;

theta_d=70/180*pi;

alpha=0.001; % sigma_w^2/sigma_u^2
M_values=[2 5 10 20]; % Number of sensors
f=f0;

phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_db_values=-20:20;
Hw_dB_values=zeros(length(iSNR_db_values),length(M_values)+1);

iSNR=10.^(iSNR_db_values/10);
Hw_dB_values(:,end)=20*log10(iSNR./(iSNR+1)); % theoretical value
for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    idx_ij=find(m_mat>n_mat);
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        phiV1=phiX1/iSNR;
        sigma2=phiX1/iSNR/phiV0(1,1);
        phiV=sigma2*phiV0;
        
        phiY=phiX1*d*d'+phiV;
        GammaY=phiY/phiY(1,1);
        
        Hw_dB_values(idxS,idxM)=20*log10((d'*(GammaY-Gamma0)*d)/(M^2-d'*Gamma0*d));
    end
end

figure
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,5),'-b','linewidth',2*linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-50 0]);
box on; grid on;

alpha=0.005; % sigma_w^2/sigma_u^2
M_values=[2 5 10 20]; % Number of sensors
f=f0;

phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_db_values=-20:20;
Hw_dB_values=zeros(length(iSNR_db_values),length(M_values)+1);

iSNR=10.^(iSNR_db_values/10);
Hw_dB_values(:,end)=20*log10(iSNR./(iSNR+1)); % theoretical value
for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    idx_ij=find(m_mat>n_mat);
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        phiV1=phiX1/iSNR;
        sigma2=phiX1/iSNR/phiV0(1,1);
        phiV=sigma2*phiV0;
        
        phiY=phiX1*d*d'+phiV;
        GammaY=phiY/phiY(1,1);
        
        Hw_dB_values(idxS,idxM)=20*log10((d'*(GammaY-Gamma0)*d)/(M^2-d'*Gamma0*d));
    end
end

figure
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,5),'-b','linewidth',2*linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-50 0]);
box on; grid on;

alpha=0.01; % sigma_w^2/sigma_u^2
M_values=[2 5 10 20]; % Number of sensors
f=f0;

phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_db_values=-20:20;
Hw_dB_values=zeros(length(iSNR_db_values),length(M_values)+1);

iSNR=10.^(iSNR_db_values/10);
Hw_dB_values(:,end)=20*log10(iSNR./(iSNR+1)); % theoretical value
for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    idx_ij=find(m_mat>n_mat);
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        phiV1=phiX1/iSNR;
        sigma2=phiX1/iSNR/phiV0(1,1);
        phiV=sigma2*phiV0;
        
        phiY=phiX1*d*d'+phiV;
        GammaY=phiY/phiY(1,1);
        
        Hw_dB_values(idxS,idxM)=20*log10((d'*(GammaY-Gamma0)*d)/(M^2-d'*Gamma0*d));
    end
end

figure
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,5),'-b','linewidth',2*linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-50 0]);
box on; grid on;


alpha=0.02; % sigma_w^2/sigma_u^2
M_values=[2 5 10 20]; % Number of sensors
f=f0;

phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_db_values=-20:20;
Hw_dB_values=zeros(length(iSNR_db_values),length(M_values)+1);

iSNR=10.^(iSNR_db_values/10);
Hw_dB_values(:,end)=20*log10(iSNR./(iSNR+1)); % theoretical value
for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    idx_ij=find(m_mat>n_mat);
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        phiV1=phiX1/iSNR;
        sigma2=phiX1/iSNR/phiV0(1,1);
        phiV=sigma2*phiV0;
        
        phiY=phiX1*d*d'+phiV;
        GammaY=phiY/phiY(1,1);
        
        Hw_dB_values(idxS,idxM)=20*log10((d'*(GammaY-Gamma0)*d)/(M^2-d'*Gamma0*d));
    end
end


figure
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Hw_dB_values(1:2:end,5),'-b','linewidth',2*linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-50 0]);
box on; grid on;

