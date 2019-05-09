close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

fname='f1.wav';
[x,Fs]=audioread(fname);  % read clean data

nfft=512;
dM=0.25*nfft;

fontsize=20;
pos1=[0.08 0.355 0.9 0.58];
pos2=[0.08 0.15 0.9 0.2];
ts=0; te=3.35;
fs=0; fe=Fs/2;
RdB=65;

[S,t,f]=spec(x,Fs,[ts te],[fs fe],nfft,nfft-dM);
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
len=length(x);
tt=(0:len-1)/Fs;
plot(tt,10*x);
gca2=gca;
set(gca2,'YTickLabel',[],'XTick',xtick,'xlim',[t(1) t(end)])
ylim_time=get(gca2,'ylim');
set(gca1,'position',pos1,'fontsize',fontsize)
set(gca2,'position',pos2,'fontsize',fontsize)
set(gca2,'YTick',[-1 0 1],'ylim',ylim_time)

