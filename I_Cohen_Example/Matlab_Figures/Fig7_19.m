close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
theta1=45/180*pi;
theta2=90/180*pi;
delta=1e-2;
alpha_values=[1e-5 1e-3 1];

M_values=[6 8 10 12]; % Number of sensors

f=linspace(0.001,8e3,15)';  % frequency

W1_dB_values=zeros(length(f),length(M_values));
W2_dB_values=zeros(length(f),length(M_values));
W3_dB_values=zeros(length(f),length(M_values));
W4_dB_values=zeros(length(f),length(M_values));

i_c=[1 0 0]';
for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        C=[exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d)) exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta1)) exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta2))];
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        alpha=alpha_values(1);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;
        W1_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        alpha=alpha_values(2);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;
        W2_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        alpha=alpha_values(3);
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\C)/(C'/Gamma_alpha*C)*i_c;
        W3_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        % Delay-Sum
        h=d/M;
        W4_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
    end
end

figure
plot(f/1e3,W1_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W2_dB_values(:,1),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W3_dB_values(:,1),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W4_dB_values(:,1),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-40 20]);
box on; grid on;

figure
plot(f/1e3,W1_dB_values(:,2),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W2_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W3_dB_values(:,2),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W4_dB_values(:,2),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-40 20]);
box on; grid on;

figure
plot(f/1e3,W1_dB_values(:,3),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W2_dB_values(:,3),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W3_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W4_dB_values(:,3),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-40 20]);
box on; grid on;

figure
plot(f/1e3,W1_dB_values(:,4),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W2_dB_values(:,4),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W3_dB_values(:,4),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W4_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'YLim',[-40 20]);
box on; grid on;

