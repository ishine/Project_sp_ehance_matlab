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

broadband_iSNR_dB=-5;
broadband_iSNR=10.^(broadband_iSNR_dB/10);
sigma_v2=phi_x1/phi_u1/(1+alpha)./broadband_iSNR;

win=hamming(3)*hamming(10)';
win=win./sum(win(:));

beta=1;
M=1; % Number of sensors

w=(alpha*phi_u1)^0.5*randn(Nx,M);
W=zeros(M,K,R);
for m=1:M
    W(m,:,:)=stft(w(:,m),nfft,dM,dN,wintype);
end
W=reshape(W,M,K*R);

[r_mat,k_mat]=meshgrid(1:R,(0:K-1)/nfft);
rk=k_mat(:);
[rk_mat,m_mat]=meshgrid(rk,1:M);
d=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_x));
du=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_u));
phiV0=zeros(M,M,length(rk));
for idx_rk=1:length(rk)
    phiV0(:,:,idx_rk)=phiU10(idx_rk)*du(:,idx_rk)*du(:,idx_rk)'+alpha*phi_u1*sigma2_stft*eye(M);
end
y1=x1+sigma_v2^0.5*v0;
Y1=stft(y1,nfft,dM,dN,wintype);
Y1a2_smooth=conv2(abs(Y1).^2,win,'same');
Y1a2_smooth=Y1a2_smooth(:);

phiV1=sigma_v2*phiV0(1,1,:);
phiV1=phiV1(:);
phiX1=max(Y1a2_smooth-beta*phiV1,0);

phiV=sigma_v2*phiV0;

h=zeros(M,length(rk));
Z=zeros(length(rk),1);
for idx_rk=1:length(rk)
    phiY=phiV(:,:,idx_rk)+phiX1(idx_rk)*d(:,idx_rk)*d(:,idx_rk)';
    h(:,idx_rk)=phiX1(idx_rk)*inv(phiY)*d(:,idx_rk);
    Z(idx_rk)=h(:,idx_rk)'*(X1(idx_rk)*d(:,idx_rk)+sigma_v2*U1(idx_rk)*du(:,idx_rk)+sigma_v2*W(:,idx_rk) );
end
Z=reshape(Z,[K R]);
z=istft(Z,nfft,dM,dN,wintype);

int_H2_phiV=0;
for idx_rk=1:length(rk)
    int_H2_phiV=int_H2_phiV+ h(:,idx_rk)'*phiV(:,:,idx_rk)*h(:,idx_rk);
end

fontsize=20;
pos1=[0.08 0.355 0.9 0.58];
pos2=[0.08 0.15 0.9 0.2];
ts=0; te=3.35;
fs=0; fe=Fs/2;
RdB=65;


[S,t,f]=spec(y1,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(2,1,1)
imagesc(t,f,S)
gca1=gca;
colormap(1-gray)
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(2,1,2)
len=length(x1);
tt=(0:len-1)/Fs;
plot(tt,10*y1);
gca2=gca;
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)

[S,t,f]=spec(z,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(2,1,1)
imagesc(t,f,S)
gca1=gca;
colormap(1-gray)
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(2,1,2)
len=length(z);
tt=(0:len-1)/Fs;
plot(tt,10*z);
gca2=gca;
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)

M=2; % Number of sensors

w=(alpha*phi_u1)^0.5*randn(Nx,M);
W=zeros(M,K,R);
for m=1:M
    W(m,:,:)=stft(w(:,m),nfft,dM,dN,wintype);
end
W=reshape(W,M,K*R);

[r_mat,k_mat]=meshgrid(1:R,(0:K-1)/nfft);
rk=k_mat(:);
[rk_mat,m_mat]=meshgrid(rk,1:M);
d=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_x));
du=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_u));
phiV0=zeros(M,M,length(rk));
for idx_rk=1:length(rk)
    phiV0(:,:,idx_rk)=phiU10(idx_rk)*du(:,idx_rk)*du(:,idx_rk)'+alpha*phi_u1*sigma2_stft*eye(M);
end
y1=x1+sigma_v2^0.5*v0;
Y1=stft(y1,nfft,dM,dN,wintype);
Y1a2_smooth=conv2(abs(Y1).^2,win,'same');
Y1a2_smooth=Y1a2_smooth(:);

phiV1=sigma_v2*phiV0(1,1,:);
phiV1=phiV1(:);
phiX1=max(Y1a2_smooth-beta*phiV1,0);
phiV=sigma_v2*phiV0;

h=zeros(M,length(rk));
Z=zeros(length(rk),1);
for idx_rk=1:length(rk)
    phiY=phiV(:,:,idx_rk)+phiX1(idx_rk)*d(:,idx_rk)*d(:,idx_rk)';
    h(:,idx_rk)=phiX1(idx_rk)*inv(phiY)*d(:,idx_rk);
    Z(idx_rk)=h(:,idx_rk)'*(X1(idx_rk)*d(:,idx_rk)+sigma_v2*U1(idx_rk)*du(:,idx_rk)+sigma_v2*W(:,idx_rk) );
end
Z=reshape(Z,[K R]);
z=istft(Z,nfft,dM,dN,wintype);

int_H2_phiV=0;
for idx_rk=1:length(rk)
    int_H2_phiV=int_H2_phiV+ h(:,idx_rk)'*phiV(:,:,idx_rk)*h(:,idx_rk);
end

[S,t,f]=spec(z,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(2,1,1)
imagesc(t,f,S)
gca1=gca;
colormap(1-gray)
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(2,1,2)
len=length(z);
tt=(0:len-1)/Fs;
plot(tt,10*z);
gca2=gca;
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)

M=5; % Number of sensors

w=(alpha*phi_u1)^0.5*randn(Nx,M);
W=zeros(M,K,R);
for m=1:M
    W(m,:,:)=stft(w(:,m),nfft,dM,dN,wintype);
end
W=reshape(W,M,K*R);

[r_mat,k_mat]=meshgrid(1:R,(0:K-1)/nfft);
rk=k_mat(:);
[rk_mat,m_mat]=meshgrid(rk,1:M);
d=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_x));
du=exp(-1i*2*pi*rk_mat.*(m_mat-1)*cos(theta_u));
phiV0=zeros(M,M,length(rk));
for idx_rk=1:length(rk)
    phiV0(:,:,idx_rk)=phiU10(idx_rk)*du(:,idx_rk)*du(:,idx_rk)'+alpha*phi_u1*sigma2_stft*eye(M);
end
y1=x1+sigma_v2^0.5*v0;
Y1=stft(y1,nfft,dM,dN,wintype);
Y1a2_smooth=conv2(abs(Y1).^2,win,'same');
Y1a2_smooth=Y1a2_smooth(:);

phiV1=sigma_v2*phiV0(1,1,:);
phiV1=phiV1(:);
phiX1=max(Y1a2_smooth-beta*phiV1,0);
phiV=sigma_v2*phiV0;

h=zeros(M,length(rk));
Z=zeros(length(rk),1);
for idx_rk=1:length(rk)
    phiY=phiV(:,:,idx_rk)+phiX1(idx_rk)*d(:,idx_rk)*d(:,idx_rk)';
    h(:,idx_rk)=phiX1(idx_rk)*inv(phiY)*d(:,idx_rk);
    Z(idx_rk)=h(:,idx_rk)'*(X1(idx_rk)*d(:,idx_rk)+sigma_v2*U1(idx_rk)*du(:,idx_rk)+sigma_v2*W(:,idx_rk) );
end
Z=reshape(Z,[K R]);
z=istft(Z,nfft,dM,dN,wintype);

int_H2_phiV=0;
for idx_rk=1:length(rk)
    int_H2_phiV=int_H2_phiV+ h(:,idx_rk)'*phiV(:,:,idx_rk)*h(:,idx_rk);
end

[S,t,f]=spec(z,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(2,1,1)
imagesc(t,f,S)
gca1=gca;
colormap(1-gray)
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(2,1,2)
len=length(z);
tt=(0:len-1)/Fs;
plot(tt,10*z);
gca2=gca;
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)


