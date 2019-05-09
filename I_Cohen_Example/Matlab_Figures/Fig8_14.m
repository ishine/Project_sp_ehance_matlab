close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
A0=0.5;
T=512;
f0=0.1;

theta_d=70/180*pi;
M_values=[3 6 9 12]; % Number of sensors
f=f0;

phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0+eps))./sin(pi*(f-f0+eps))).^2;

iSNR_dB=-10;
iSNR=10^(iSNR_dB/10);

theta_values=(40:5:100)'/180*pi;

alpha=1e-20; % sigma_w^2/sigma_u^2
s=1253;

Ks=30;

J_values=zeros(length(theta_values),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    C=cholcov(phiV0);
    rng(s);
    V0= C'*randn(M,Ks);
    
    phiV1=phiX1/iSNR;
    sigma2=phiX1/iSNR/phiV0(1,1);
    V=sigma2^0.5*V0;
    phiV=V*V'/Ks;
    
    phiY=phiX1*d*d'+phiV;
    GammaY=phiY/phiY(1,1);
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    [T,Lambda]=jeig(GammaY,Gamma0);
    [diagL,idx]=sort(diag(Lambda),'descend');
    Lambda=diag(diagL);
    T=T(:,idx);
    
    for idxS=1:length(theta_values)
        theta=theta_values(idxS);
        d_theta=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta));
        for m=2:M
            J_values(idxS,idxM)=J_values(idxS,idxM)+abs(T(:,m)'*d_theta)^2;
        end
    end
end

figure
plot(theta_values/pi*180,J_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(theta_values/pi*180,J_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XTick',40:10:100);
box on; grid on;


Ks=50;
J_values=zeros(length(theta_values),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    C=cholcov(phiV0);
    rng(s);
    V0= C'*randn(M,Ks);
    
    phiV1=phiX1/iSNR;
    sigma2=phiX1/iSNR/phiV0(1,1);
    V=sigma2^0.5*V0;
    phiV=V*V'/Ks;
    
    phiY=phiX1*d*d'+phiV;
    GammaY=phiY/phiY(1,1);
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    [T,Lambda]=jeig(GammaY,Gamma0);
    [diagL,idx]=sort(diag(Lambda),'descend');
    Lambda=diag(diagL);
    T=T(:,idx);
    
    for idxS=1:length(theta_values)
        theta=theta_values(idxS);
        d_theta=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta));
        for m=2:M
            J_values(idxS,idxM)=J_values(idxS,idxM)+abs(T(:,m)'*d_theta)^2;
        end
    end
end

figure
plot(theta_values/pi*180,J_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(theta_values/pi*180,J_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XTick',40:10:100);
box on; grid on;

Ks=300;
J_values=zeros(length(theta_values),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    C=cholcov(phiV0);
    rng(s);
    V0= C'*randn(M,Ks);
    
    phiV1=phiX1/iSNR;
    sigma2=phiX1/iSNR/phiV0(1,1);
    V=sigma2^0.5*V0;
    phiV=V*V'/Ks;
    
    phiY=phiX1*d*d'+phiV;
    GammaY=phiY/phiY(1,1);
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    [T,Lambda]=jeig(GammaY,Gamma0);
    [diagL,idx]=sort(diag(Lambda),'descend');
    Lambda=diag(diagL);
    T=T(:,idx);
    
    for idxS=1:length(theta_values)
        theta=theta_values(idxS);
        d_theta=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta));
        for m=2:M
            J_values(idxS,idxM)=J_values(idxS,idxM)+abs(T(:,m)'*d_theta)^2;
        end
    end
end

figure
plot(theta_values/pi*180,J_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(theta_values/pi*180,J_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XTick',40:10:100);
box on; grid on;
print -depsc2 fig8_14c.eps

Ks=1000;
J_values=zeros(length(theta_values),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    
    I_M=eye(M);
    i_i=I_M(:,1);
    d=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta_d));
    
    phiV0=(1-alpha)*sinc(2*f*(m_mat-n_mat))+alpha*eye(M); % diffuse noise
    C=cholcov(phiV0);
    rng(s);
    V0= C'*randn(M,Ks);
    
    phiV1=phiX1/iSNR;
    sigma2=phiX1/iSNR/phiV0(1,1);
    V=sigma2^0.5*V0;
    phiV=V*V'/Ks;
    
    phiY=phiX1*d*d'+phiV;
    GammaY=phiY/phiY(1,1);
    
    Gamma0=sinc(2*f*(m_mat-n_mat));
    [T,Lambda]=jeig(GammaY,Gamma0);
    [diagL,idx]=sort(diag(Lambda),'descend');
    Lambda=diag(diagL);
    T=T(:,idx);
    
    for idxS=1:length(theta_values)
        theta=theta_values(idxS);
        d_theta=exp(-1i*2*pi*f.*(0:M-1)'*cos(theta));
        for m=2:M
            J_values(idxS,idxM)=J_values(idxS,idxM)+abs(T(:,m)'*d_theta)^2;
        end
    end
end

figure
plot(theta_values/pi*180,J_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(theta_values/pi*180,J_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(theta_values/pi*180,J_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XTick',40:10:100);
box on; grid on;

