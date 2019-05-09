close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;

delta=0.005;
M=6;     % number of microphones
[m_mat,n_mat]=meshgrid(1:M,1:M);
a_values=[0 0.5 0.9 0.99];
f=linspace(100,8e3,15)';  % frequency

D1_dB_values=zeros(length(f),length(a_values));
W2_dB_values=zeros(length(f),length(a_values));

for idxa=1:length(a_values)
    aleph=a_values(idxa);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        b1=0.5*[3^0.5-1; 3-3^0.5];  % supercardioid
        f_bar=2*pi*delta/c*fk*(0:1)';
        b0bar=besselj(0,f_bar);
        b1bar=2i*besselj(1,f_bar);
        B1bar=[b0bar.'; b1bar.'];
        hNR=B1bar\b1;
        H=toeplitz([conj(hNR(1)); zeros(M-2,1)],[hNR' zeros(1,M-2)])';
        f_bar=2*pi*delta/c*fk*(0:M-1)';
        dt=H'*d;
        [m_mat,n_mat]=meshgrid(1:M,1:M);
        GammaC=besseli(0,2i*pi*fk*(m_mat-n_mat)*delta/c);
        R=H'*GammaC*H;
        b0bar=besselj(0,-f_bar);
        b1bar=1i*besselj(1,-f_bar);
        Gamma_dpc=[b0bar b1bar];
        Ra=aleph*R+(1-aleph)*H'*H;
        g_Ua=aleph*(Ra\H')*Gamma_dpc*b1;
        g_Ta=g_Ua+(1-dt'*g_Ua)/(dt'/Ra*dt)*(Ra\dt);
        h=H*g_Ta;
        D1_dB_values(idx_f,idxa)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
        W2_dB_values(idx_f,idxa)=10*log10(abs(h'*d)^2/(h'*h));
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
set(gca,'YLim',[-30 10]);
box on; grid on;

