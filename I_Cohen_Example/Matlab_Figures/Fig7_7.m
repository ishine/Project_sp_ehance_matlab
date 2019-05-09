close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta1_d=90/180*pi;
delta1=3e-2;

theta2_d=0/180*pi;
delta2=1e-2;

theta3_d=0/180*pi;
delta3=3e-2;

M_values=[2 4 6 8]; % Number of sensors

f=linspace(0,8e3,15)';  % frequency

W1_dB_values=zeros(length(f),length(M_values));
W2_dB_values=zeros(length(f),length(M_values));
W3_dB_values=zeros(length(f),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta1/c*(0:M-1)'*cos(theta1_d));
        Gamma0=sinc(2*fk*delta1/c*(m_mat-n_mat));
        h=(Gamma0\d)/(d'/Gamma0*d);
        W1_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        d=exp(-1i*2*pi*fk*delta2/c*(0:M-1)'*cos(theta2_d));
        Gamma0=sinc(2*fk*delta2/c*(m_mat-n_mat));
        h=(Gamma0\d)/(d'/Gamma0*d);
        W2_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        d=exp(-1i*2*pi*fk*delta3/c*(0:M-1)'*cos(theta3_d));
        Gamma0=sinc(2*fk*delta3/c*(m_mat-n_mat));
        h=(Gamma0\d)/(d'/Gamma0*d);
        W3_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
    end
end

figure
plot(f/1e3,W1_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W1_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W1_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W1_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-100 20]);
box on; grid on;

figure
plot(f/1e3,W2_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W2_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W2_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W2_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-100 20]);
box on; grid on;

figure
plot(f/1e3,W3_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,W3_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W3_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,W3_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'YLim',[-100 20]);
box on; grid on;

