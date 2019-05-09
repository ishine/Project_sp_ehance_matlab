close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
delta=1e-2;
f_vec=[1 2 4 8]*1e3;  % frequencies

M_values=[2 4 6 8]; % Number of sensors

alpha_values=logspace(-2,0,15)';  % frequency

D1_dB_values=zeros(length(alpha_values),length(M_values));
D2_dB_values=zeros(length(alpha_values),length(M_values));
D3_dB_values=zeros(length(alpha_values),length(M_values));
D4_dB_values=zeros(length(alpha_values),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    for idx_a=1:length(alpha_values)
        alpha=alpha_values(idx_a);
        fk=f_vec(1);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D1_dB_values(idx_a,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        fk=f_vec(2);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D2_dB_values(idx_a,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        fk=f_vec(3);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D3_dB_values(idx_a,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        fk=f_vec(4);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        Gamma_alpha=(1-alpha)*Gamma0+alpha*eye(M);
        h=(Gamma_alpha\d)/(d'/Gamma_alpha*d);
        D4_dB_values(idx_a,idxM)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
    end
end

figure
semilogx(alpha_values,D1_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
semilogx(alpha_values,D1_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D1_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D1_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[0 15]);
box on; grid on;

figure
semilogx(alpha_values,D2_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
semilogx(alpha_values,D2_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D2_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D2_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[0 15]);
box on; grid on;

figure
semilogx(alpha_values,D3_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
semilogx(alpha_values,D3_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D3_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D3_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[0 15]);
box on; grid on;

figure
semilogx(alpha_values,D4_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
semilogx(alpha_values,D4_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D4_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
semilogx(alpha_values,D4_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[0 15]);
box on; grid on;


