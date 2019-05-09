close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A0=0.5;
f0=0.1;
alpha=0.01; % sigma_w^2/sigma_u^2
theta_values=[0 30 60 90]/180*pi;
sigmav2=0.1;
M=12;
T0=500;
f=f0;

phiX1=A0^2/4*T0^2;
phiV0=T0*alpha*eye(M);
[mm,nn]=meshgrid(1:M,1:M);
phiV0=phiV0+T0*sinc(2*f*(mm-nn));
phiV=sigmav2*phiV0;
phiV=phiV.*(1-eye(M))+real(diag(diag(phiV)));

N_values=1:M;
gain=zeros(length(N_values),length(theta_values));
Wgain=zeros(length(N_values),length(theta_values));
for idx2=1:length(theta_values)
    theta=theta_values(idx2);
    d=exp(-1i*2*pi*f*(0:M-1)'*cos(theta));
    phiX=phiX1*(d*d');
    phiX=phiX.*(1-eye(M))+real(diag(diag(phiX)));
    
    [T,Lambda]=jeig(phiX,phiV);
    [diagL,idx]=sort(abs(diag(Lambda)),'descend');
    Lambda=diag(diagL);
    T=T(:,idx);

    for idx1=1:length(N_values)
        N=N_values(idx1);
        
        T1N=T(:,1:N);
        P1N=T1N*pinv(T1N'*T1N)*T1N';
        gain(idx1,idx2)=10*log10(abs(d'*P1N'*d)^2/(d'*P1N'*(phiV/phiV(1,1))*P1N*d));
        Wgain(idx1,idx2)=10*log10(abs(d'*P1N'*d)^2/(d'*P1N'*P1N*d));
    end
end

figure
plot(N_values(1:1:end),Wgain(1:1:end,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(N_values(1:1:end),Wgain(1:1:end,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(N_values(1:1:end),Wgain(1:1:end,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(N_values(1:1:end),Wgain(1:1:end,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(N_values)]);
box on; grid on;

figure
plot(N_values(1:1:end),gain(1:1:end,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(N_values(1:1:end),gain(1:1:end,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(N_values(1:1:end),gain(1:1:end,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(N_values(1:1:end),gain(1:1:end,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 10]);
box on; grid on;
