close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

c=340; % voice velocity
fs=8e3;
theta_d=0/180*pi;
theta_1=90/180*pi;   % dipole
P=10;

Lh_vec=(1:3:60)';
M_values=[4 6 8 10]; % Number of sensors
delta_values=[0.01 0.01 0.01 0.01]; % distance between microphones

W_dB_values=zeros(length(Lh_vec),length(M_values));
D_dB_values=zeros(length(Lh_vec),length(M_values));
for idxM=1:length(M_values)
    for idxD=1:length(Lh_vec)
        
        Lh=Lh_vec(idxD);
        M=M_values(idxM);
        delta=delta_values(idxM);
        L=2*P+Lh;
        [jj,ii]=meshgrid(1:L,1:Lh);
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
        
        phi_vec=(0:1:180)'/180*pi;
        Lphi=length(phi_vec);
        GGT_phi=zeros(Lphi,M*Lh,M*Lh);
        G_phi=zeros(M*Lh,L);
        for idx_phi=1:Lphi
            phi=phi_vec(idx_phi);
            for m=1:M
                Gm_phi=sinc(-P-Lh+1-ii+jj-fs*(m-1)*delta/c*cos(phi));
                G_phi(((m-1)*Lh+1):(m*Lh),:)=Gm_phi;
            end
            GGT_phi(idx_phi,:,:)=G_phi*G_phi';
        end
        dphi=phi_vec(2)-phi_vec(1);
        GammaT_0_pi=zeros(M*Lh,M*Lh);
        for idx1=1:M*Lh
            for idx2=1:M*Lh
                GammaT_0_pi(idx1,idx2)=sum(GGT_phi(:,idx1,idx2).*sin(phi_vec)*dphi);
            end
        end
        
        W=(h'*(G_theta_d*G_theta_d')*h)/(h'*h);
        W_dB_values(idxD,idxM)=10*log10(W);
        D=(h'*(G_theta_d*G_theta_d')*h)/(h'*GammaT_0_pi*h);
        D_dB_values(idxD,idxM)=10*log10(D);
    end
end

figure
plot(Lh_vec,W_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(Lh_vec,W_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(Lh_vec,W_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(Lh_vec,W_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[Lh_vec(1) Lh_vec(end)]);
box on; grid on;

figure
plot(Lh_vec,D_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(Lh_vec,D_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(Lh_vec,D_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(Lh_vec,D_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[Lh_vec(1) Lh_vec(end)]);
box on; grid on;

