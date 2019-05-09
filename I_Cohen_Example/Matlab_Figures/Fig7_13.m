close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
delta=1e-2;
alpha_values=[0.001 0.01 0.1 1];

M_values=[2 4 6 8]; % Number of sensors

f=linspace(0,8e3,15)';  % frequency

D1_dB_values=zeros(length(f),length(M_values));
D2_dB_values=zeros(length(f),length(M_values));
D3_dB_values=zeros(length(f),length(M_values));
D4_dB_values=zeros(length(f),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        alpha=alpha_values(1);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D1_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        alpha=alpha_values(2);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D2_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        alpha=alpha_values(3);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D3_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        alpha=alpha_values(4);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
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
set(gca,'YLim',[0 15]);
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
set(gca,'YLim',[0 15]);
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
set(gca,'YLim',[0 15]);
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
set(gca,'YLim',[0 15]);
box on; grid on;
