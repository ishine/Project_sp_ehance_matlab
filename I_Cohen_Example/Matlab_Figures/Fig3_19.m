close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

fname='f6.wav';
[x,Fs]=audioread(fname);  % read clean data
phi_x=mean(x.^2);

v0=randn(size(x));

nfft=512;
dM=0.25*nfft;
dN=1;
wintype='hamming';
X=stft(x,nfft,dM,dN,wintype);
V0=stft(v0,nfft,dM,dN,wintype);
phiV0=mean(abs(V0(:)).^2);

broadband_iSNR_dB=(-5:15)';
broadband_iSNR=10.^(broadband_iSNR_dB/10);
sigma_v2_values=phi_x./broadband_iSNR;
N=length(broadband_iSNR_dB);

win=hamming(3)*hamming(10)';
win=win./sum(win(:));

beta_values=[1 2 3 4]';
Nb=length(beta_values);
gain=zeros(N,Nb);
MSE=zeros(N,Nb);
xi_n=zeros(N,Nb);
xi_d=zeros(N,Nb);
for k=1:Nb
    beta=beta_values(k);
    for idx=1:N
        sigma_v2=sigma_v2_values(idx);
        y=x+sigma_v2^0.5*v0;
        Y=stft(y,nfft,dM,dN,wintype);
        Ya2_smooth=conv2(abs(Y).^2,win,'same');

        phiV=sigma_v2*phiV0;
        %phiX=max(abs(Y).^2-phiV,0);
        phiX=max(Ya2_smooth-beta*phiV,0);
        narrowband_iSNR=phiX/phiV;
        H_W=narrowband_iSNR./(1+narrowband_iSNR);
        
        int_phiX=mean(phiX(:));
        int_phiV=phiV;
        int_H2_phiX=mean(H_W(:).^2 .* phiX(:));
        int_H2_phiV=mean(H_W(:).^2 .* phiV);
        int_Hm12_phiX=mean(abs(H_W(:)-1).^2 .* phiX(:));
        
        oSNR=int_H2_phiX./int_H2_phiV;
        oSNR_dB=10*log10(oSNR);
        gain(idx,k)=oSNR_dB-broadband_iSNR_dB(idx);
        MSE(idx,k)=10*log10(int_H2_phiV+int_Hm12_phiX);
        xi_n(idx,k)=10*log10(int_phiV./int_H2_phiV);
        xi_d(idx,k)=10*log10(int_phiX./int_H2_phiX);
    end
end

figure
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),gain(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
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
plot(broadband_iSNR_dB(1:2:end),MSE(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),MSE(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
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
plot(broadband_iSNR_dB(1:2:end),xi_n(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_n(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
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
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,2),'--r*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,3),':gs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(broadband_iSNR_dB(1:2:end),xi_d(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;
