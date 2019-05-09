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

fontsize=20;
pos1=[0.08 0.358 0.8 0.58];
pos2=[0.08 0.15 0.8 0.2];
pos3=[0.888 0.15 0.05 0.785];
ts=0; te=2.28;
fs=0; fe=Fs/2;
RdB=50;

broadband_iSNR_dB=10;
broadband_iSNR=10.^(broadband_iSNR_dB/10);
sigma_v2=phi_x./broadband_iSNR;

y=x+sigma_v2^0.5*v0;
Y=stft(y,nfft,dM,dN,wintype);

[S,t,f]=spec(y,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(3,1,1)
imagesc(t,f,S)
gca1=gca;
ylabel('tag1','fontsize',fontsize);
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(3,1,2)
h2=imagesc(t,f,S);
cbar_axes=colorbar;
set(h2,'visible','off')
set(gca,'visible','off')
subplot(3,1,3)
len=length(y);
tt=(0:len-1)/Fs;
plot(tt,8*y);
gca2=gca;
xlabel('tag2','fontsize',fontsize);
ylabel('tag3','fontsize',fontsize);
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)
set(cbar_axes,'position',pos3,'fontsize',fontsize)

phiV=sigma_v2*phiV0;
win=hamming(3)*hamming(10)';
win=win./sum(win(:));
Ya2_smooth=conv2(abs(Y).^2,win,'same');
phiX=max(Ya2_smooth-1*phiV,0);
narrowband_iSNR=phiX/phiV;
H_W=narrowband_iSNR./(1+narrowband_iSNR);
Z=H_W.*Y;
z=istft(Z,nfft,dM,dN,wintype);

[S,t,f]=spec(z,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(3,1,1)
imagesc(t,f,S)
gca1=gca;
ylabel('tag1','fontsize',fontsize);
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(3,1,2)
h2=imagesc(t,f,S);
cbar_axes=colorbar;
set(h2,'visible','off')
set(gca,'visible','off')
subplot(3,1,3)
len=length(y);
tt=(0:len-1)/Fs;
plot(tt,8*y);
gca2=gca;
xlabel('tag2','fontsize',fontsize);
ylabel('tag3','fontsize',fontsize);
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)
set(cbar_axes,'position',pos3,'fontsize',fontsize)

phiX=max(Ya2_smooth-2*phiV,0);
narrowband_iSNR=phiX/phiV;
H_W=narrowband_iSNR./(1+narrowband_iSNR);
Z=H_W.*Y;
z=istft(Z,nfft,dM,dN,wintype);

[S,t,f]=spec(z,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(3,1,1)
imagesc(t,f,S)
gca1=gca;
ylabel('tag1','fontsize',fontsize);
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(3,1,2)
h2=imagesc(t,f,S);
cbar_axes=colorbar;
set(h2,'visible','off')
set(gca,'visible','off')
subplot(3,1,3)
len=length(y);
tt=(0:len-1)/Fs;
plot(tt,8*y);
gca2=gca;
xlabel('tag2','fontsize',fontsize);
ylabel('tag3','fontsize',fontsize);
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)
set(cbar_axes,'position',pos3,'fontsize',fontsize)

phiX=max(Ya2_smooth-3*phiV,0);
narrowband_iSNR=phiX/phiV;
H_W=narrowband_iSNR./(1+narrowband_iSNR);
Z=H_W.*Y;
z=istft(Z,nfft,dM,dN,wintype);

[S,t,f]=spec(z,Fs,[ts te],[fs fe],nfft,nfft-dM);
t=t-t(1);
f=f/1e3;
% plot spectrum
MaxS=max(S(:));
S=S-MaxS;
S(S<-RdB)=-RdB;

figure
subplot(3,1,1)
imagesc(t,f,S)
gca1=gca;
ylabel('tag1','fontsize',fontsize);
set(gca1,'ydir','normal','XTickLabel',[])
xtick=get(gca1,'XTick');
subplot(3,1,2)
h2=imagesc(t,f,S);
cbar_axes=colorbar;
set(h2,'visible','off')
set(gca,'visible','off')
subplot(3,1,3)
len=length(y);
tt=(0:len-1)/Fs;
plot(tt,8*y);
gca2=gca;
xlabel('tag2','fontsize',fontsize);
ylabel('tag3','fontsize',fontsize);
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)
set(cbar_axes,'position',pos3,'fontsize',fontsize)

