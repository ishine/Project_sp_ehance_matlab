close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A0=0.5;
T0=500;
f0=0.1;
theta=120/180*pi;
alpha=0.01; 

M=6;
N_values=[1 2 3 4]; % Number of sensors
i_i=[1;zeros(M-1,1)];

f=f0; 
phiX1=A0^2/4*T0^2;

[f_mat,m_mat]=meshgrid(f,1:M);
d=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta));

phiV0=zeros(M,M,length(f));
[mm,nn]=meshgrid(1:M,1:M);
for idxf=1:length(f)
    phiV0(:,:,idxf)=T0*sinc(2*f(idxf)*(mm-nn))+T0*alpha*eye(M);
end

iSNR_db_values=-5:15;  

Gain_dB_values=zeros(length(iSNR_db_values),length(N_values));
MSE_values=zeros(length(iSNR_db_values),length(N_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(N_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(N_values));
for idxN=1:length(N_values)
    N=N_values(idxN);
    
    phiV=phiV0;
    idxf=1;
    phiX=phiX1(idxf)*(d(:,idxf)*d(:,idxf)');
    phiX=phiX.*(1-eye(M))+real(diag(diag(phiX)));
    phiV(:,:,idxf)=phiV(:,:,idxf).*(1-eye(M))+real(diag(diag(phiV(:,:,idxf))));
    [T,Lambda]=jeig(phiX,phiV(:,:,idxf));
    [diagL,idx]=sort(diag(Lambda),'descend');
    Lambda=diag(diagL);
    T=T(:,idx);
    
    T1N=T(:,1:N);
    P1N=T1N*inv(T1N'*T1N)*T1N';

    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        sigma2=phiX1./phiV0(1,1)./iSNR;
        phiV=sigma2*phiV0;
        
        h=zeros(M,length(f));
        for idxf=1:length(f)
            h(:,idxf)=P1N*(d(:,idxf)*d(:,idxf)')*P1N*(phiX+phiV(:,:,idxf))*i_i/(d(:,idxf)'*P1N*(phiX+phiV(:,:,idxf))*P1N*d(:,idxf));
        end
        int_phiX1=phiX1;
        int_phiV1=phiV(1,1,:);
        int_H2_phiX=phiX1 .* abs(sum(conj(h).*d,1)').^2 ;
        int_H2_phiV=0;
        for idxf=1:length(f)
            int_H2_phiV=int_H2_phiV+ h(:,idxf)'*phiV(:,:,idxf)*h(:,idxf);
        end
        int_Hm12_phiX= phiX1 .* abs(sum(conj(h).*d,1)'-1).^2;
        
        oSNR=int_H2_phiX./int_H2_phiV;
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxN)=oSNR_dB-iSNR_dB;
        
        MSE_values(idxS,idxN)=10*log10(int_H2_phiV+int_Hm12_phiX);
        
        xi_n_dB_values(idxS,idxN)=10*log10(int_phiV1./int_H2_phiV);
        xi_d_dB_values(idxS,idxN)=10*log10(int_phiX1./int_H2_phiX);
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

