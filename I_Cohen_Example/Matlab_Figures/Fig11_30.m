close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

% Input:
c=340; % voice velocity
fs=8e3;
theta_d=0/180*pi;
theta_1=180/180*pi;   % cardioid
M_vec=[4 6 8 10];     % number of microphones
delta_vec=[0.01 0.01 0.01 0.01]; % distance between microphones
P=10;
Lh_vec=[25 25 25 25];

idx_fig=1;
Lh=Lh_vec(idx_fig);
L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);
M=M_vec(idx_fig);
delta=delta_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;
G_theta_1=zeros(M*Lh,L);
for m=1:M
    Gm_theta_1=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_1));
    G_theta_1(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_1;
end
C=[G_theta_d G_theta_1];
i_1=[i_ell; zeros(L,1)];
h=C*pinv(C'*C+1e-4*eye(2*L))*i_1;

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
    Ba2(idx_phi)=h'*(G_phi*G_phi')*h;
end
Ba=Ba2.^0.5;
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
Lh=Lh_vec(idx_fig);
L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);
M=M_vec(idx_fig);
delta=delta_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;
G_theta_1=zeros(M*Lh,L);
for m=1:M
    Gm_theta_1=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_1));
    G_theta_1(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_1;
end
C=[G_theta_d G_theta_1];
i_1=[i_ell; zeros(L,1)];
h=C*pinv(C'*C+1e-4*eye(2*L))*i_1;

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
    Ba2(idx_phi)=h'*(G_phi*G_phi')*h;
end
Ba=Ba2.^0.5;
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
Lh=Lh_vec(idx_fig);
L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);
M=M_vec(idx_fig);
delta=delta_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;
G_theta_1=zeros(M*Lh,L);
for m=1:M
    Gm_theta_1=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_1));
    G_theta_1(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_1;
end
C=[G_theta_d G_theta_1];
i_1=[i_ell; zeros(L,1)];
h=C*pinv(C'*C+1e-4*eye(2*L))*i_1;

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
    Ba2(idx_phi)=h'*(G_phi*G_phi')*h;
end
Ba=Ba2.^0.5;
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
Lh=Lh_vec(idx_fig);
L=2*P+Lh;
[jj,ii]=meshgrid(1:L,1:Lh);
M=M_vec(idx_fig);
delta=delta_vec(idx_fig);
G_theta_d=zeros(M*Lh,L);
for m=1:M
    Gm_theta_d=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_d));
    G_theta_d(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_d;
end
diag_GG=diag(G_theta_d'*G_theta_d);
idx_max=find(diag_GG==max(diag_GG),1);
i_ell=zeros(L,1);
i_ell(idx_max)=1;
G_theta_1=zeros(M*Lh,L);
for m=1:M
    Gm_theta_1=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(theta_1));
    G_theta_1(((m-1)*Lh+1):(m*Lh),:)=Gm_theta_1;
end
C=[G_theta_d G_theta_1];
i_1=[i_ell; zeros(L,1)];
h=C*pinv(C'*C+1e-4*eye(2*L))*i_1;

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
    Ba2(idx_phi)=h'*(G_phi*G_phi')*h;
end
Ba=Ba2.^0.5;
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

