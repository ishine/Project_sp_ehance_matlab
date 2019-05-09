close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
M=6;     % number of microphones
[m_mat,n_mat]=meshgrid(1:M,1:M);

delta=0.005;
eps_values=[0 1e-5 1e-3 1e-1];     % number of microphones
f=linspace(100,8e3,15)';  % frequency

D1_dB_values=zeros(length(f),length(eps_values));
D2_dB_values=zeros(length(f),length(eps_values));
D3_dB_values=zeros(length(f),length(eps_values));
D4_dB_values=zeros(length(f),length(eps_values));

for idxM=1:length(eps_values)
    epsilon=eps_values(idxM);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        GammaC=besseli(0,2i*pi*fk*(m_mat-n_mat)*delta/c);
        b1=[0; 1];  % dipole
        f_bar=2*pi*delta/c*fk*(0:M-1)';
        b0bar=besselj(0,f_bar);
        b1bar=2i*besselj(1,f_bar);
        B1bar=[b0bar.'; b1bar.'];
        GammaC_eps=GammaC+epsilon*eye(M);
        h=(GammaC_eps\B1bar')/(B1bar/GammaC_eps*B1bar')*b1;
        D1_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        b1=[0.5; 0.5];  % cardioid
        f_bar=2*pi*delta/c*fk*(0:M-1)';
        b0bar=besselj(0,f_bar);
        b1bar=2i*besselj(1,f_bar);
        B1bar=[b0bar.'; b1bar.'];
        GammaC_eps=GammaC+epsilon*eye(M);
        h=(GammaC_eps\B1bar')/(B1bar/GammaC_eps*B1bar')*b1;
        D2_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        b1=[0.25; 0.75];  % hypercardioid
        f_bar=2*pi*delta/c*fk*(0:M-1)';
        b0bar=besselj(0,f_bar);
        b1bar=2i*besselj(1,f_bar);
        B1bar=[b0bar.'; b1bar.'];
        GammaC_eps=GammaC+epsilon*eye(M);
        h=(GammaC_eps\B1bar')/(B1bar/GammaC_eps*B1bar')*b1;
        D3_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        b1=0.5*[3^0.5-1; 3-3^0.5];  % supercardioid
        f_bar=2*pi*delta/c*fk*(0:M-1)';
        b0bar=besselj(0,f_bar);
        b1bar=2i*besselj(1,f_bar);
        B1bar=[b0bar.'; b1bar.'];
        GammaC_eps=GammaC+epsilon*eye(M);
        h=(GammaC_eps\B1bar')/(B1bar/GammaC_eps*B1bar')*b1;
        D4_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
    end
end

figure
plot(f/1e3,D1_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,D1_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D1_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D1_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[0 7]);
box on; grid on;

figure
plot(f/1e3,D2_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,D2_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D2_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D2_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[0 7]);
box on; grid on;

figure
plot(f/1e3,D3_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,D3_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D3_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D3_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[0 7]);
box on; grid on;

figure
plot(f/1e3,D4_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,D4_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D4_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D4_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[0 7]);
box on; grid on;

