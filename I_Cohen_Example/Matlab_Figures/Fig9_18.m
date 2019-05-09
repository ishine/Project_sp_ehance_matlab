close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
M=4;     % number of microphones
[m_mat,n_mat]=meshgrid(1:M,1:M);

delta_values=[0.01 0.02 0.03 0.04];     % distance between microphones
f=linspace(350,8e3,15)';  % frequency

D1_dB_values=zeros(length(f),length(delta_values));
D2_dB_values=zeros(length(f),length(delta_values));
D3_dB_values=zeros(length(f),length(delta_values));
D4_dB_values=zeros(length(f),length(delta_values));

for idxD=1:length(delta_values)
    delta=delta_values(idxD);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        alpha31=-1;
        h=1/(1-exp(1i*2*pi*fk*delta/c*(1-alpha31)))^3*[1; -3*exp(-1i*2*pi*fk*delta/c*alpha31);  3*exp(-1i*4*pi*fk*delta/c*alpha31);  -exp(-1i*6*pi*fk*delta/c*alpha31)];
        D1_dB_values(idx_f,idxD)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        alpha31=0;
        alpha32=-0.5;
        alpha33=-1;
        V=[exp(-1i*2*pi*fk*delta/c*(0:M-1)') exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha31) exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha32) exp(-1i*2*pi*fk*delta/c*(0:M-1)'*alpha33)]';
        h=V\[1; 0; 0; 0];
        D2_dB_values(idx_f,idxD)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        h=Gamma0\d/(d'/Gamma0*d);
        D3_dB_values(idx_f,idxD)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        psi1=0; psi2=pi/2;
        Gamma_0_pi2=(exp(1i*2*pi*fk*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*fk*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*fk*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
        psi1=pi/2; psi2=pi;
        Gamma_pi2_pi=(exp(1i*2*pi*fk*delta/c*(m_mat-n_mat+eps)*cos(psi1))-exp(1i*2*pi*fk*delta/c*(m_mat-n_mat+eps)*cos(psi2)))./(1i*2*pi*fk*delta/c*(m_mat-n_mat+eps))/(cos(psi1)-cos(psi2));
        [T,Lambda]=jeig(Gamma_0_pi2,Gamma_pi2_pi);
        [diagL,idx]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx);
        t1=T(:,1);
        h=t1/(d'*t1);
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        D4_dB_values(idx_f,idxD)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
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
set(gca,'YLim',[-10 15]);
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
set(gca,'YLim',[-10 15]);
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
set(gca,'YLim',[-10 15]);
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
set(gca,'YLim',[-10 15]);
box on; grid on;


