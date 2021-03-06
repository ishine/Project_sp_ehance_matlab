close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

% Input:
c=340; % voice velocity
fs=8e3;
theta_d=0/180*pi;
delta=0.01;    % distance between microphones
M=10;    % number of microphones
eps_vec=[1e-5 1e-3 1e-1 1];

P=25;
Lh=30;
L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);

idx_fig=1;
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
GGT_phi=zeros(Lphi,M*Lh,M*Lh);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    GGT_phi(idx_phi,:,:)=G_phi*G_phi';
end
dphi=phi_vec(2)-phi_vec(1);
GammaT_0_pi=zeros(M*Lh,M*Lh);
for idx1=1:M*Lh
    for idx2=1:M*Lh
        GammaT_0_pi(idx1,idx2)=sum(GGT_phi(:,idx1,idx2).*sin(phi_vec)*dphi);
    end
end
GammaT_0_pi_eps=GammaT_0_pi+eps_vec(idx_fig)*eye(M*Lh);
h=pinv(GammaT_0_pi_eps)*G_theta_d*pinv(G_theta_d'*pinv(GammaT_0_pi_eps)*G_theta_d)*i_ell;

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
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180);
set(gca,'YLim',[-30 0]);
box on; grid on;

idx_fig=2;
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
GGT_phi=zeros(Lphi,M*Lh,M*Lh);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    GGT_phi(idx_phi,:,:)=G_phi*G_phi';
end
dphi=phi_vec(2)-phi_vec(1);
GammaT_0_pi=zeros(M*Lh,M*Lh);
for idx1=1:M*Lh
    for idx2=1:M*Lh
        GammaT_0_pi(idx1,idx2)=sum(GGT_phi(:,idx1,idx2).*sin(phi_vec)*dphi);
    end
end
GammaT_0_pi_eps=GammaT_0_pi+eps_vec(idx_fig)*eye(M*Lh);
h=pinv(GammaT_0_pi_eps)*G_theta_d*pinv(G_theta_d'*pinv(GammaT_0_pi_eps)*G_theta_d)*i_ell;

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
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180);
set(gca,'YLim',[-30 0]);
box on; grid on;

idx_fig=3;
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
GGT_phi=zeros(Lphi,M*Lh,M*Lh);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    GGT_phi(idx_phi,:,:)=G_phi*G_phi';
end
dphi=phi_vec(2)-phi_vec(1);
GammaT_0_pi=zeros(M*Lh,M*Lh);
for idx1=1:M*Lh
    for idx2=1:M*Lh
        GammaT_0_pi(idx1,idx2)=sum(GGT_phi(:,idx1,idx2).*sin(phi_vec)*dphi);
    end
end
GammaT_0_pi_eps=GammaT_0_pi+eps_vec(idx_fig)*eye(M*Lh);
h=pinv(GammaT_0_pi_eps)*G_theta_d*pinv(G_theta_d'*pinv(GammaT_0_pi_eps)*G_theta_d)*i_ell;

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
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180);
set(gca,'YLim',[-30 0]);
box on; grid on;

idx_fig=4;
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;

phi_vec=(0:1:180)'/180*pi;
Lphi=length(phi_vec);
GGT_phi=zeros(Lphi,M*Lh,M*Lh);
G_phi=zeros(M*Lh,L);
for idx_phi=1:Lphi
    phi=phi_vec(idx_phi);
    for m=1:M
        Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
        G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
    end
    GGT_phi(idx_phi,:,:)=G_phi*G_phi';
end
dphi=phi_vec(2)-phi_vec(1);
GammaT_0_pi=zeros(M*Lh,M*Lh);
for idx1=1:M*Lh
    for idx2=1:M*Lh
        GammaT_0_pi(idx1,idx2)=sum(GGT_phi(:,idx1,idx2).*sin(phi_vec)*dphi);
    end
end
GammaT_0_pi_eps=GammaT_0_pi+eps_vec(idx_fig)*eye(M*Lh);
h=pinv(GammaT_0_pi_eps)*G_theta_d*pinv(G_theta_d'*pinv(GammaT_0_pi_eps)*G_theta_d)*i_ell;

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
Ba=Ba/max(Ba);
Ba_dB=20*log10(Ba);

figure
plot(phi_vec/pi*180,Ba_dB,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[0 180]);
set(gca,'XTick',0:30:180);
set(gca,'YLim',[-30 0]);
box on; grid on;

