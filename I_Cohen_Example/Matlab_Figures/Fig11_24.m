close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

% Input:
P=20;
Lh=30;
beta=0.1; % sigma_w^2/sigma_u^2
alpha=0.8;
theta_d=0/180*pi;
theta_n=90/180*pi;
c=340; % voice velocity
fs=8e3;
iSNR_dB=0;

delta_vec=[0.03 0.03 0.03 0.03];     % distance between microphones
M_vec=[4 6 10 15];     % number of microphones
L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);

idx_fig=1;
delta=delta_vec(idx_fig);
M=M_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

[t1,t2]=meshgrid(1:L,1:L);
Rx=alpha.^(abs(t1-t2));
Rvw=beta*eye(Lh*M);
Rvi=kron(ones(M),eye(Lh));
Rv0=Rvi+Rvw;
iSNR=10^(iSNR_dB/10);
sigma_i2=Rx(1,1)/Rv0(1,1)/iSNR;
Rv_bar=sigma_i2*Rv0;

[T,Lambda]=jeig(G_theta_d*Rx*G_theta_d',Rv_bar);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);

zeta=(t1'*G_theta_d*Rx*i_ell)/(t1'*(G_theta_d*Rx*G_theta_d'+Rv_bar)*t1);
h=zeta*t1; 

% Beampattern
phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
Ba2=zeros(Lphi,1);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    Ba2(idx_phi)=h'*G_phi*G_phi'*h;
end
Ba=Ba2.^0.5;
Ba2=Ba2/max(Ba2);
Ba_dB=10*log10(Ba2);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180); 
set(gca,'YLim',[-50 0]);
box on; grid on;

idx_fig=2;
delta=delta_vec(idx_fig);
M=M_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

[t1,t2]=meshgrid(1:L,1:L);
Rx=alpha.^(abs(t1-t2));
Rvw=beta*eye(Lh*M);
Rvi=kron(ones(M),eye(Lh));
Rv0=Rvi+Rvw;
iSNR=10^(iSNR_dB/10);
sigma_i2=Rx(1,1)/Rv0(1,1)/iSNR;
Rv_bar=sigma_i2*Rv0;

[T,Lambda]=jeig(G_theta_d*Rx*G_theta_d',Rv_bar);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);

zeta=(t1'*G_theta_d*Rx*i_ell)/(t1'*(G_theta_d*Rx*G_theta_d'+Rv_bar)*t1);
h=zeta*t1; 

% Beampattern
phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
Ba2=zeros(Lphi,1);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    Ba2(idx_phi)=h'*G_phi*G_phi'*h;
end
Ba=Ba2.^0.5;
Ba2=Ba2/max(Ba2);
Ba_dB=10*log10(Ba2);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180); 
set(gca,'YLim',[-50 0]);
box on; grid on;

idx_fig=3;
delta=delta_vec(idx_fig);
M=M_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

[t1,t2]=meshgrid(1:L,1:L);
Rx=alpha.^(abs(t1-t2));
Rvw=beta*eye(Lh*M);
Rvi=kron(ones(M),eye(Lh));
Rv0=Rvi+Rvw;
iSNR=10^(iSNR_dB/10);
sigma_i2=Rx(1,1)/Rv0(1,1)/iSNR;
Rv_bar=sigma_i2*Rv0;

[T,Lambda]=jeig(G_theta_d*Rx*G_theta_d',Rv_bar);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);

zeta=(t1'*G_theta_d*Rx*i_ell)/(t1'*(G_theta_d*Rx*G_theta_d'+Rv_bar)*t1);
h=zeta*t1; 

% Beampattern
phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
Ba2=zeros(Lphi,1);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    Ba2(idx_phi)=h'*G_phi*G_phi'*h;
end
Ba=Ba2.^0.5;
Ba2=Ba2/max(Ba2);
Ba_dB=10*log10(Ba2);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180); 
set(gca,'YLim',[-50 0]);
box on; grid on;

idx_fig=4;
delta=delta_vec(idx_fig);
M=M_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

[t1,t2]=meshgrid(1:L,1:L);
Rx=alpha.^(abs(t1-t2));
Rvw=beta*eye(Lh*M);
Rvi=kron(ones(M),eye(Lh));
Rv0=Rvi+Rvw;
iSNR=10^(iSNR_dB/10);
sigma_i2=Rx(1,1)/Rv0(1,1)/iSNR;
Rv_bar=sigma_i2*Rv0;

[T,Lambda]=jeig(G_theta_d*Rx*G_theta_d',Rv_bar);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);
t1=T(:,1);

zeta=(t1'*G_theta_d*Rx*i_ell)/(t1'*(G_theta_d*Rx*G_theta_d'+Rv_bar)*t1);
h=zeta*t1; 

% Beampattern
phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
Ba2=zeros(Lphi,1);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    Ba2(idx_phi)=h'*G_phi*G_phi'*h;
end
Ba=Ba2.^0.5;
Ba2=Ba2/max(Ba2);
Ba_dB=10*log10(Ba2);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180); 
set(gca,'YLim',[-50 0]);
box on; grid on;

