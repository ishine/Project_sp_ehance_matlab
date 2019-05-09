close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=12;

L=10;
K=5;
A=0.5./(1:K)';
f_vec=0.05:0.1:0.45;

[t1,t2]=meshgrid(1:L,1:L);
Rx=0.5*A(1)^2*cos(2*pi*f_vec(1)*(t1-t2));
for k=2:K
    Rx=Rx+0.5*A(k)^2*cos(2*pi*f_vec(k)*(t1-t2));
end
Rv=0.1*0.5.^abs(t1-t2);
[T,Lambda]=jeig(Rx,Rv);
[diagL,idx]=sort(diag(Lambda),'descend');
Lambda=diag(diagL);
T=T(:,idx);

iSNR=trace(Lambda)/L;
iSNR_dB=10*log10(iSNR);

pidx_values=1:L;  % index for different P
oSNR_values=zeros(length(pidx_values),2);
MSE_values=zeros(length(pidx_values),2);

for P=1:length(pidx_values)
    h=diagL./(1+diagL);
    h(P+1:L)=0;
    
    oSNR=(h'*Lambda*h)/(h'*h);
    oSNR_values(P,1)=10*log10(oSNR)-iSNR_dB;
    
    MSE=sum((1-h).^2.*diagL+h.^2);
    MSE_values(P,1)=10*log10(MSE);
    
    h=ones(L,1);
    h(P+1:L)=0;
    
    oSNR=(h'*Lambda*h)/(h'*h);
    oSNR_values(P,2)=10*log10(oSNR)-iSNR_dB;
    
    MSE=sum((1-h).^2.*diagL+h.^2);
    MSE_values(P,2)=10*log10(MSE);
end

figure
plot(pidx_values,oSNR_values(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(pidx_values,oSNR_values(:,2),'r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(pidx_values)]);
set(gca,'XTick',pidx_values); 
box on; grid on;

figure
plot(pidx_values,MSE_values(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(pidx_values,MSE_values(:,2),'r*','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(pidx_values)]);
set(gca,'XTick',pidx_values); 
set(gca,'YLim',[7 12]);
box on; grid on;

