close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A0=0.5;
f0=0.1;
L=30;

t=(0:550)';
x=A0*cos(2*pi*f0*t+pi/3);

M=1; % Number of sensors
[t1,t2]=meshgrid(1:L,1:L);
Rx1=0.5*A0^2*cos(2*pi*f0*(t1-t2));

Rvw=0.01*eye(L*M);
Rvi=eye(L*M);
for m1=0:M-1
    for m2=0:M-1
        if L>abs(m2-m1)
            Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
        end
    end
end
Rv0=Rvi+Rvw;
Rx=kron(ones(M),Rx1);

iSNR_dB=-5;
iSNR=10^(iSNR_dB/10);

sigma_i2=A0^2/2/(1.01*iSNR);
Rv=sigma_i2*Rv0;
Ry=Rx+Rv;

I_i=[eye(L) zeros(L,(M-1)*L)];
hW=I_i*Rx/Ry;

vw=(sigma_i2/100)^0.5*randn(length(t),M);
vi0=sigma_i2^0.5*randn(length(t),1);
vi=hankel(vi0,[vi0(end) zeros(1,M-1)]);
y=x*ones(1,M)+vw+vi;

Y=zeros(M*L,length(t));
for k=1:length(t)-L
    Y(:,k)=reshape(y(k:(k+L-1),:),L*M,1);
end
Z=hW*Y;
z=Z(1,:);

xx=x*ones(1,M);
X=zeros(M*L,length(t));
for k=1:length(t)-L
    X(:,k)=reshape(xx(k:(k+L-1),:),L*M,1);
end
V=zeros(M*L,length(t));
for k=1:length(t)-L
    V(:,k)=reshape(vw(k:(k+L-1),:)+vi(k:(k+L-1),:),L*M,1);
end

figure
plot(t,y(:,1),'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
set(gca,'XLim',[1 500]);
set(gca,'YLim',[-2 2]);
box on; grid on;

figure
plot(t,z,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 500]);
set(gca,'YLim',[-1 1]);
box on; grid on;

M=2; % Number of sensors
[t1,t2]=meshgrid(1:L,1:L);
Rx1=0.5*A0^2*cos(2*pi*f0*(t1-t2));

Rvw=0.01*eye(L*M);
Rvi=eye(L*M);
for m1=0:M-1
    for m2=0:M-1
        if L>abs(m2-m1)
            Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
        end
    end
end
Rv0=Rvi+Rvw;
Rx=kron(ones(M),Rx1);

iSNR_dB=-5;
iSNR=10^(iSNR_dB/10);

sigma_i2=A0^2/2/(1.01*iSNR);
Rv=sigma_i2*Rv0;
Ry=Rx+Rv;

I_i=[eye(L) zeros(L,(M-1)*L)];
hW=I_i*Rx/Ry;

vw=(sigma_i2/100)^0.5*randn(length(t),M);
vi0=sigma_i2^0.5*randn(length(t),1);
vi=hankel(vi0,[vi0(end) zeros(1,M-1)]);
y=x*ones(1,M)+vw+vi;

Y=zeros(M*L,length(t));
for k=1:length(t)-L
    Y(:,k)=reshape(y(k:(k+L-1),:),L*M,1);
end
Z=hW*Y;
z=Z(1,:);

xx=x*ones(1,M);
X=zeros(M*L,length(t));
for k=1:length(t)-L
    X(:,k)=reshape(xx(k:(k+L-1),:),L*M,1);
end
V=zeros(M*L,length(t));
for k=1:length(t)-L
    V(:,k)=reshape(vw(k:(k+L-1),:)+vi(k:(k+L-1),:),L*M,1);
end

figure
plot(t,z,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 500]);
set(gca,'YLim',[-1 1]);
box on; grid on;

M=5; % Number of sensors
[t1,t2]=meshgrid(1:L,1:L);
Rx1=0.5*A0^2*cos(2*pi*f0*(t1-t2));

Rvw=0.01*eye(L*M);
Rvi=eye(L*M);
for m1=0:M-1
    for m2=0:M-1
        if L>abs(m2-m1)
            Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
        end
    end
end
Rv0=Rvi+Rvw;
Rx=kron(ones(M),Rx1);

iSNR_dB=-5;
iSNR=10^(iSNR_dB/10);

sigma_i2=A0^2/2/(1.01*iSNR);
Rv=sigma_i2*Rv0;
Ry=Rx+Rv;

I_i=[eye(L) zeros(L,(M-1)*L)];
hW=I_i*Rx/Ry;

vw=(sigma_i2/100)^0.5*randn(length(t),M);
vi0=sigma_i2^0.5*randn(length(t),1);
vi=hankel(vi0,[vi0(end) zeros(1,M-1)]);
y=x*ones(1,M)+vw+vi;

Y=zeros(M*L,length(t));
for k=1:length(t)-L
    Y(:,k)=reshape(y(k:(k+L-1),:),L*M,1);
end
Z=hW*Y;
z=Z(1,:);

xx=x*ones(1,M);
X=zeros(M*L,length(t));
for k=1:length(t)-L
    X(:,k)=reshape(xx(k:(k+L-1),:),L*M,1);
end
V=zeros(M*L,length(t));
for k=1:length(t)-L
    V(:,k)=reshape(vw(k:(k+L-1),:)+vi(k:(k+L-1),:),L*M,1);
end

figure
plot(t,z,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 500]);
set(gca,'YLim',[-1 1]);
box on; grid on;


