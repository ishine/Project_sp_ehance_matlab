close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

K=5;
A0=0.5./(1:K)';
f_vec=0.05:0.1:0.45;
L=10;
M_values=[2 5 10 20]; % Number of sensors

[t2,t1]=meshgrid(1:L,1:L);
Rx1=0.5*A0(1)^2*cos(2*pi*f_vec(1)*(t1-t2));
for k=2:length(f_vec)
    Rx1=Rx1+0.5*A0(k)^2*cos(2*pi*f_vec(k)*(t1-t2));
end

Rvw1=0.01*eye(L);
Rvi1=0.1*0.5.^abs(t1-t2);
Rv01=Rvi1+Rvw1;

[T10,Lambda10]=jeig(Rx1,Rv01);
[diagL10,idx1]=sort(diag(Lambda10),'descend');
Lambda10=diag(diagL10);
T10=T10(:,idx1);

iSNR_db_values=-5:15;
Gain_dB_values=zeros(length(iSNR_db_values),length(M_values));
MSE_values=zeros(length(iSNR_db_values),length(M_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(M_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    i_i=[1; zeros(M-1,1)];
    Rxl0=ones(M,M,L);  % spatial correlation matrix
    for l=1:L
        Rxl0(:,:,l)=diagL10(l);
    end
    
    Rvw=kron(eye(M),Rvw1);
    Rvi=eye(L*M);
    for m1=0:M-1
        for m2=0:M-1
            Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=0.1*0.5.^abs(t1-t2+m1-m2);
        end
    end
    Rv0=Rvi+Rvw;
    
    Rvl=zeros(M,M,L);  % spatial correlation matrix
    for m1=0:M-1
        for m2=0:M-1
            Rvl(m1+1,m2+1,:)=diag(T10'*Rv0((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))*T10);
        end
    end
    
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        sigma_i2=trace(Rx1/Rv01)/L/iSNR;

        Rxl=1/sigma_i2*Rxl0;  % spatial correlation matrix
        sum_hRxh=0;
        sum_hRvh=0;
        sum_hRxi=0;
        for l=1:L
            [S,Omega]=jeig(Rxl(:,:,l),Rvl(:,:,l));
            [diagL,idx1]=sort(diag(Omega),'descend');
            Omega=diag(diagL);
            S=S(:,idx1);
            
            h=(Rxl(:,:,l)+Rvl(:,:,l))\Rxl(:,:,l)*i_i;
            
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


















