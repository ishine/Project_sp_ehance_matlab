close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
delta=0.005;     % distance between microphones

M_values=[4 6 8 10];     % distance between microphones
f=linspace(200,8e3,15)';  % frequency

W1_dB_values=zeros(length(f),length(M_values));
W2_dB_values=zeros(length(f),length(M_values));

for idxM=1:length(M_values)
    M=M_values(idxM);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        alpha31=0;
        alpha32=-0.5;
        alpha33=-1;
        D=[exp(-1i*2*pi*fk*delta/c*(0:M-1)')'; exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha31)'; exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha32)'; exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha33)'];
        h=D'/(D*D')*[1 0 0 0]';
        W1_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
        alpha31=-1;
        S=diag(0:M-1);
        D=[exp(-1i*2*pi*fk*delta/c*(0:M-1)')'; exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha31)'; (S*exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha31))'; (S^2*exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha31))'];
        h=D'/(D*D')*[1 0 0 0]';
        W2_dB_values(idx_f,idxM)=10*log10(abs(h'*d)^2/(h'*h));
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
set(gca,'YLim',[-50 10]);
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
set(gca,'YLim',[-50 10]);
box on; grid on;

