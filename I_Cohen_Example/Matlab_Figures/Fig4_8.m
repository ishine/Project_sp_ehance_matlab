close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

L=10;
M=5; % Number of sensors
mu_values=[0.5 1 2 5];

[t1,t2]=meshgrid(1:L,1:L);
alpha=0.8;
Rx1=alpha.^(abs(t1-t2));
Rx=kron(ones(M),Rx1);

Rvw=0.01*eye(L*M);
Rvi=eye(L*M);
for m1=0:M-1
    for m2=0:M-1
        if L>abs(m2-m1)
            Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
        end
    end
end
Rv0=Rvi+Rvw;

iSNR_db_values=-5:15;
Gain_dB_values=zeros(length(iSNR_db_values),length(mu_values));
MSE_values=zeros(length(iSNR_db_values),length(mu_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(mu_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(mu_values));
for idx_mu=1:length(mu_values)
    mu=mu_values(idx_mu);
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        sigma_i2=trace(Rx1)/trace(Rv0(1:L,1:L))/iSNR;
        Rv=sigma_i2*Rv0;
        
        [T,Lambda]=jeig(Rx,Rv);
        [diagL,idx1]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx1);
        
        I_i=[eye(L) zeros(L,(M-1)*L)];
        hW=0;
        for q=1:(M*L)
            hW=hW+I_i*Rx*(T(:,q)*T(:,q)')/(mu+diagL(q));
        end
        oSNR=trace(hW*Rx*hW')/trace(hW*Rv*hW');
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idx_mu)=oSNR_dB-iSNR_dB;
        
        MSE_values(idxS,idx_mu)=10*log10( trace(Rx1-2*hW*Rx*I_i'+hW*(Rx+Rv)*hW'));
        
        xi_n=trace(Rv(1:L,1:L))/trace(hW*Rv*hW');
        xi_n_dB_values(idxS,idx_mu)=10*log10(xi_n);
        
        xi_d=trace(Rx1)/trace(hW*Rx*hW');
        xi_d_dB_values(idxS,idx_mu)=10*log10(xi_d);
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

















