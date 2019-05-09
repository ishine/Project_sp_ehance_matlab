close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A0=0.5;
T=500;
f0=0.1;
theta=70/180*pi;
alpha=0.01; 

iSNR_dB=-5;
iSNR=10^(iSNR_dB/10);
sigmau2=A0^2/2./iSNR/(1+alpha);
sigmaw2=alpha*sigmau2;

Nfft=2^9;
Nfft2=Nfft/2;
f=[0:Nfft2 -Nfft2+1:-1]'/Nfft;
df=f(2)-f(1);

t=(0:T-1)';
x1=A0*cos(2*pi*f0*t+pi/3);
u1=sigmau2^0.5*randn(size(t));
X1=fft(x1,Nfft);
U1=fft(u1,Nfft);
phiX1=A0^2/4*(sin(T*pi*(f+f0))./sin(pi*(f+f0))).^2+A0^2/4*(sin(T*pi*(f-f0))./sin(pi*(f-f0))).^2;

M=1; % Number of sensors

[f_mat,m_mat]=meshgrid(f,1:M);
d=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta));
du=exp(-1i*2*pi*f_mat.*(m_mat-1));
X=(ones(M,1)*X1.').*d;
U=(ones(M,1)*U1.').*du;
w=sigmaw2^0.5*randn(M,T);
W=fft(w,Nfft,2);
Y=X+U+W;
y1=ifft(Y(1,:));

phiV0=zeros(M,M,length(f));
for idxf=1:length(f)
    phiV0(:,:,idxf)=T*du(:,idxf)*du(:,idxf)'+T*alpha*eye(M);
end
phiV=sigmau2*phiV0;

h=zeros(M,length(f));
for idxf=1:length(f)
    phiY=phiV(:,:,idxf)+phiX1(idxf)*d(:,idxf)*d(:,idxf)';
    h(:,idxf)=phiX1(idxf)*inv(phiY)*d(:,idxf);
end

Z=sum(conj(h).*Y,1);

int_H2_phiV=0;
for idxf=1:length(f)
    int_H2_phiV=int_H2_phiV+df* h(:,idxf)'*phiV(:,:,idxf)*h(:,idxf);
end

f1=fftshift(f);
Ya=fftshift(abs(Y(1,:)));

figure
plot(f1(2:end),Ya(2:end),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',-0.5:0.1:0.5); 
set(gca,'YLim',[0 140]);
box on; grid on;

Ea=fftshift(abs(Z.'-X1));

figure
plot(f1(2:end),Ea(2:end),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',-0.5:0.1:0.5); 
set(gca,'YLim',[0 20]);
box on; grid on;

M=2; % Number of sensors
[f_mat,m_mat]=meshgrid(f,1:M);
d=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta));
du=exp(-1i*2*pi*f_mat.*(m_mat-1));
X=(ones(M,1)*X1.').*d;
U=(ones(M,1)*U1.').*du;
w=sigmaw2^0.5*randn(M,T);
W=fft(w,Nfft,2);
Y=X+U+W;

phiV0=zeros(M,M,length(f));
for idxf=1:length(f)
    phiV0(:,:,idxf)=T*du(:,idxf)*du(:,idxf)'+T*alpha*eye(M);
end
phiV=sigmau2*phiV0;

h=zeros(M,length(f));
for idxf=1:length(f)
    phiY=phiV(:,:,idxf)+phiX1(idxf)*d(:,idxf)*d(:,idxf)';
    h(:,idxf)=phiX1(idxf)*inv(phiY)*d(:,idxf);
end

Z=sum(conj(h).*Y);

int_H2_phiV=0;
for idxf=1:length(f)
    int_H2_phiV=int_H2_phiV+df* h(:,idxf)'*phiV(:,:,idxf)*h(:,idxf);
end

figure
Ea=fftshift(abs(Z.'-X1));
plot(f1(2:end),Ea(2:end),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',-0.5:0.1:0.5); 
set(gca,'YLim',[0 20]);
box on; grid on;

M=5; % Number of sensors
[f_mat,m_mat]=meshgrid(f,1:M);
d=exp(-1i*2*pi*f_mat.*(m_mat-1)*cos(theta));
du=exp(-1i*2*pi*f_mat.*(m_mat-1));
X=(ones(M,1)*X1.').*d;
U=(ones(M,1)*U1.').*du;
w=sigmaw2^0.5*randn(M,T);
W=fft(w,Nfft,2);
Y=X+U+W;

phiV0=zeros(M,M,length(f));
for idxf=1:length(f)
    phiV0(:,:,idxf)=T*du(:,idxf)*du(:,idxf)'+T*alpha*eye(M);
end
phiV1=T*A0^2/2./iSNR*ones(size(f));
phiV=sigmau2*phiV0;

h=zeros(M,length(f));
for idxf=1:length(f)
    phiY=phiV(:,:,idxf)+phiX1(idxf)*d(:,idxf)*d(:,idxf)';
    h(:,idxf)=phiX1(idxf)*inv(phiY)*d(:,idxf);
end

Z=sum(conj(h).*Y);
z=real(ifft(Z));

int_H2_phiX=df*sum(phiX1 .* abs(sum(conj(h).*d,1)').^2 );
int_H2_phiV=0;
for idxf=1:length(f)
    int_H2_phiV=int_H2_phiV+df* h(:,idxf)'*phiV(:,:,idxf)*h(:,idxf);
end

Za=fftshift(abs(Z));

figure
Ea=fftshift(abs(Z.'-X1));
plot(f1(2:end),Ea(2:end),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',-0.5:0.1:0.5); 
set(gca,'YLim',[0 20]);
box on; grid on;


