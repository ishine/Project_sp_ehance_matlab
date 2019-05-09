close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

theta_x=70/180*pi;
theta_u=20/180*pi;
alpha=0.1;

fname='f1.wav';
[x1,Fs]=audioread(fname);  % read clean data
Nx=length(x1);
phi_x1=mean(x1.^2);

u1=randn(Nx,1);
phi_u1=mean(u1.^2);

v0=u1+(alpha*phi_u1)^0.5*randn(Nx,1);

nfft=512;
dM=0.25*nfft;
dN=1;
wintype='hamming';
X1=stft(x1,nfft,dM,dN,wintype);
[K,R]=size(X1);
U1=stft(u1,nfft,dM,dN,wintype);
phiU10_k=mean(abs(U1).^2,2);
phiU10=phiU10_k*ones(1,R);
phiU10=phiU10(:);

w1=randn(Nx,1);
W1=stft(w1,nfft,dM,dN,wintype);
sigma2_stft=mean(abs(W1(:)).^2);

broadband_iSNR_dB=(-5:15)';
broadband_iSNR=10.^(broadband_iSNR_dB/10);
sigma_v2_values=phi_x1/phi_u1/(1+alpha)./broadband_iSNR;
N=length(broadband_iSNR_dB);

win=hamming(3)*hamming(10)';
win=win./sum(win(:));

beta=1;
M_values=[1 2 5 10]; % Number of sensors
gain=zeros(N,length(M_values));
MSE=zeros(N,length(M_values));
xi_n=zeros(N,length(M_values));
xi_d=zeros(N,length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    
    [r_mat,k_mat]=meshgrid(1:R,(0:K-1)/nfft);
    rk=k_mat(:);
    [rk_mat,m_mat]=meshgrid(rk,1:M);
    d=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_x));
    du=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_u));
    phiV0=zeros(M,M,length(rk));
    for idx_rk=1:length(rk)
        phiV0(:,:,idx_rk)=phiU10(idx_rk)*du(:,idx_rk)*du(:,idx_rk)'+alpha*phi_u1*sigma2_stft*eye(M);
    end
    
    for idxS=1:N
        
        sigma_v2=sigma_v2_values(idxS);
        y1=x1+sigma_v2^0.5*v0;
        Y1=stft(y1,nfft,dM,dN,wintype);
        Y1a2_smooth=conv2(abs(Y1).^2,win,'same');
        Y1a2_smooth=Y1a2_smooth(:);

        phiV1=sigma_v2*phiV0(1,1,:);
        phiV1=phiV1(:);
        phiX1=max(Y1a2_smooth-beta*phiV1,0);
        phiV=sigma_v2*phiV0;
        
        h=zeros(M,length(rk));
        for idx_rk=1:length(rk)
            phiY=phiV(:,:,idx_rk)+phiX1(idx_rk)*d(:,idx_rk)*d(:,idx_rk)';
            h(:,idx_rk)=phiX1(idx_rk)*(phiY\d(:,idx_rk));
        end
        
        int_phiX1=sum(phiX1);
        int_phiV1=sum(phiV1);
        int_H2_phiX=sum(phiX1 .* abs(sum(conj(h).*d,1)').^2 );
        int_H2_phiV=0;
        for idx_rk=1:length(rk)
            int_H2_phiV=int_H2_phiV+ h(:,idx_rk)'*phiV(:,:,idx_rk)*h(:,idx_rk);
        end
        int_Hm12_phiX=sum( phiX1 .* abs(sum(conj(h).*d,1)'-1).^2);
        
        oSNR=int_H2_phiX./int_H2_phiV;
        oSNR_dB=10*log10(oSNR);
        gain(idxS,idxM)=oSNR_dB-broadband_iSNR_dB(idxS);
        MSE(idxS,idxM)=10*log10(int_H2_phiV+int_Hm12_phiX);
        xi_n(idxS,idxM)=10*log10(int_phiV1./int_H2_phiV);
        xi_d(idxS,idxM)=10*log10(int_phiX1./int_H2_phiX);
    end
end

figure
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(broadband_iSNR_dB(1:2:end),MSE(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(broadband_iSNR_dB(1:2:end),MSE(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),MSE(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),MSE(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(broadband_iSNR_dB(1:2:end),xi_n(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(broadband_iSNR_dB(1:2:end),xi_n(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_n(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_n(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

