function [Sy,t,f]=spec(y,Fs,tI,fI,M,Mo)	% compute spectrum in dB
win=hamming(M);	% window function
% time interval
ls=round((tI(1)*Fs-M/2)/(M-Mo));
ls=max(ls,0);
le=round((tI(2)*Fs-M/2)/(M-Mo));
t=((ls:le)*(M-Mo)+M/2)/Fs;
% frequency interval
ks=round(fI(1)*M/Fs);
fe=min(fI(2),Fs/2);
ke=round(fe*M/Fs);
f=(ks:ke)*Fs/M;
Sy=zeros(ke-ks+1,le-ls+1);
start=ls*(M-Mo);
for l=1:le-ls+1
    % 1. Short Time Fourier Analysis
    Y=fft(win.*y(start+(1:M)));
    Sy(:,l)=20*log10(abs(Y(ks+1:ke+1)));
    start=start+M-Mo;
end
