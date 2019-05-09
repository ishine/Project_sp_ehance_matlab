close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
M=3;     % number of microphones
[m_mat,n_mat]=meshgrid(1:M,1:M);

delta_values=[0.01 0.02 0.03 0.04];     % distance between microphones
f=linspace(0,8e3,15)';  % frequency

D1_dB_values=zeros(length(f),length(delta_values));
D2_dB_values=zeros(length(f),length(delta_values));

for idxD=1:length(delta_values)
    delta=delta_values(idxD);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        alpha21=-1; % dipole
        h=1/(1-exp(1i*2*pi*fk*delta/c*(1-alpha21)))*[1; -2*exp(-1i*2*pi*fk*delta/c*alpha21); exp(-1i*4*pi*fk*delta/c*alpha21)];
        D1_dB_values(idx_f,idxD)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        alpha21=-1;
        alpha22=0;
        V=[exp(-1i*2*pi*fk*delta/c*(0:M-1)')'; exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha21)'; exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha22)'];
        h=V\[1; 0; 0];
        D2_dB_values(idx_f,idxD)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
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
set(gca,'YLim',[-10 10]);
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
set(gca,'YLim',[-10 10]);
box on; grid on;
