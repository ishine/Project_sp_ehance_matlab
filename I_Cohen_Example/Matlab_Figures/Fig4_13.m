close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

K=3;
A0=0.5./(1:K)';
f0=0.1;
f_vec=0.05:0.1:0.45;
L=10;
M=10; % Number of sensors
mu_values=[0.5 1 2 5];
theta_values=[60 60 60 60 60]/180*pi;  % directions of desired signals

iSNR_db_values=-5:15;
Gain_dB_values=zeros(length(iSNR_db_values),length(mu_values));
MSE_values=zeros(length(iSNR_db_values),length(mu_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(mu_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(mu_values));

[t2,t1]=meshgrid(1:L,1:L);
for idxM=1:length(mu_values)
    mu=mu_values(idxM);
    
    Rx=zeros(M*L);
    for m1=0:M-1
        for m2=0:M-1
            for k1=1:K
                Rx((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=Rx((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))+0.5*A0(k1)*A0(k1)*cos(2*pi*f_vec(k1)*(t1-t2+m1*cos(theta_values(k1))-m2*cos(theta_values(k1))));
            end
        end
    end
    
    Rvw=0.01*eye(L*M);
    Rvi=zeros(L*M);
    for m1=0:M-1
        for m2=0:M-1
            Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=0.1*0.5.^abs(t1-t2+m1-m2);
        end
    end
    Rv0=Rvi+Rvw;
    
    T0=zeros(L,L,M);
    Lambda0=zeros(L,L,M);
    for m=0:M-1
        [T,Lambda]=jeig(Rx((m*L+1):(m*L+L),(m*L+1):(m*L+L)),Rv0((m*L+1):(m*L+L),(m*L+1):(m*L+L)));
        [diagL,idx1]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx1);
        T0(:,:,m+1)=T;
        Lambda0(:,:,m+1)=Lambda;
    end
    
    i_i=[1; zeros(M-1,1)];
    Rxl0=zeros(M,M,L);  % spatial correlation matrix
    for m1=0:M-1
        for m2=0:M-1
            for l=1:L
                Rxl0(m1+1,m2+1,l)=T0(:,l,m1+1)'*Rx((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))*T0(:,l,m2+1);
            end
        end
    end
    
    Rvl=zeros(M,M,L);  % spatial correlation matrix
    for m1=0:M-1
        for m2=0:M-1
            for l=1:L
                Rvl(m1+1,m2+1,l)=T0(:,l,m1+1)'*Rv0((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))*T0(:,l,m2+1);
            end
        end
    end
   
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        sigma_i2=trace(Lambda0(:,:,1))/L/iSNR;
        
        Rxl=1/sigma_i2*Rxl0;  % spatial correlation matrix
        sum_hRxh=0;
        sum_hRvh=0;
        sum_hRxi=0;
        for l=1:L
            [S,Omega]=jeig(Rxl(:,:,l),Rvl(:,:,l));
            [diagL,idx1]=sort(diag(Omega),'descend');
            Omega=diag(diagL);
            S=S(:,idx1);
            
            h=zeros(M,1);
            for m=1:M
                h=h+Omega(m,m)/(mu+Omega(m,m))*S(:,m)*S(:,m)'*Rvl(:,:,l)*i_i;
            end

            sum_hRxh=sum_hRxh+h'*Rxl(:,:,l)*h;
            sum_hRvh=sum_hRvh+h'*Rvl(:,:,l)*h;
            sum_hRxi=sum_hRxi+h'*Rxl(:,:,l)*i_i;
        end
        
        oSNR=sum_hRxh/sum_hRvh;
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxM)=oSNR_dB-iSNR_dB;
        
        MSE=sum(Rxl(1,1,:))+sum_hRxh+sum_hRvh-2*sum_hRxi;
        MSE_values(idxS,idxM)=10*log10(MSE);
        
        xi_n=L/sum_hRvh;
        xi_n_dB_values(idxS,idxM)=10*log10(xi_n);
        
        xi_d=sum(Rxl(1,1,:))/sum_hRxh;
        xi_d_dB_values(idxS,idxM)=10*log10(xi_d);
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


















