close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=5;

A0=0.5;
f0=0.1;
sigma_i2=0.5;
sigma_w2=sigma_i2/100;

L_values=1:60;  % length of the filter
M_values=[1 2 5 10]; % Number of sensors
Gain_dB_values=zeros(length(L_values),length(M_values));
MSE_values=zeros(length(L_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    for idx=1:length(L_values)
        L=L_values(idx);
        
        Rvw=sigma_w2*eye(L*M);
        Rvi=eye(L*M);
        for m1=0:M-1
            for m2=0:M-1
                if L>abs(m2-m1)
                    Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
                end
            end
        end
        Rvi=sigma_i2*Rvi;
        Rv=Rvi+Rvw;
        
        [t1,t2]=meshgrid(1:L,1:L);
        Rx1=0.5*A0^2*cos(2*pi*f0*(t1-t2));
        Rx=kron(ones(M),Rx1);
        Ry=Rx+Rv;
        
        [T,Lambda]=jeig(Rx,Rv);
        [diagL,idx1]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx1);
        
        I_i=[eye(L) zeros(L,(M-1)*L)];
        A=I_i*Rx*T/(Lambda+eye(L*M));
        
        iSNR=trace(Rx1)/trace(Rv(1:L,1:L));
        iSNR_dB=10*log10(iSNR);
        oSNR=trace(A*Lambda*A')/trace(A*A');
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idx,idxM)=oSNR_dB-iSNR_dB;
        
        MSE_values(idx,idxM)=10*log10(1/L*trace(Rx1-2*A*T'*Rx*I_i'+A*(Lambda+eye(L*M))*A'));
    end
end

figure
plot(L_values,Gain_dB_values(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(L_values,Gain_dB_values(:,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,Gain_dB_values(:,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,Gain_dB_values(:,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;

figure
plot(L_values,MSE_values(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(L_values,MSE_values(:,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,MSE_values(:,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,MSE_values(:,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;
