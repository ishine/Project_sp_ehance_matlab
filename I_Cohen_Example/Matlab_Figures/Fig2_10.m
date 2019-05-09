close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A=0.5;
f_vec=0.1:0.03:0.45;
L=20; 
K=4;

hidx_values=1:4;  % index for different filters
iSNR_db_values=(-10:10)';  
oSNR_dB_values=zeros(length(iSNR_db_values),length(hidx_values));
MSE_values=zeros(length(iSNR_db_values),length(hidx_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(hidx_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(hidx_values));

i_i=[1; zeros(L-1,1)];
[t1,t2]=meshgrid(1:L,1:L);
Rx=0.5*A^2*cos(2*pi*f_vec(1)*(t1-t2));
for k=2:K
    Rx=Rx+0.5*A^2*cos(2*pi*f_vec(k)*(t1-t2));
end

xin_dB=10;
xin=10^(xin_dB/10);
fun=@(x) myfunc3(x,xin);
mu_sig2=fzero(fun,0.1);

mu_values=zeros(size(iSNR_db_values));

for idx_snr=1:length(iSNR_db_values)
    iSNR_dB=iSNR_db_values(idx_snr);
    iSNR=10^(iSNR_dB/10);
    
    sigma_v2=K*A^2/2/iSNR;
    Rv=sigma_v2*eye(L);
    Ry=Rx+Rv;
    
    % h_max
    idx1=1;
    [Tx,Xi_x]=eigs(Rv\Rx,1);
    zeta=(Tx'*Rx*i_i)/(Tx'*Ry*Tx);
    h=zeta*Tx;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
    
    % h_Wiener
    idx1=2;
    h=Ry\Rx*i_i;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
    
    % h_MVDR
    idx1=3;
    [Tx,Xi_x]=eigs(Rx,2*K);
    h=Tx*Tx'*i_i;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
    
    % h_T
    idx1=4;
    mu=mu_sig2/sigma_v2;
    mu_values(idx_snr)=mu;
    h=(Rx+mu*Rv)\Rx*i_i;
    
    oSNR=(h'*Rx*h)/(h'*Rv*h);
    oSNR_dB=10*log10(oSNR);
    oSNR_dB_values(idx_snr,idx1)=oSNR_dB;
    
    MSE_values(idx_snr,idx1)=10*log10(Rx(1,1)-2*h'*Rx*i_i+h'*Ry*h);
    
    xi_n=sigma_v2/(h'*Rv*h);
    xi_n_dB_values(idx_snr,idx1)=10*log10(xi_n);
    
    xi_d=Rx(1,1)/(h'*Rx*h);
    xi_d_dB_values(idx_snr,idx1)=10*log10(xi_d);
end

figure
semilogy(iSNR_db_values,mu_values,'linewidth',linewd, 'MarkerSize',MarkerSize);
idx0=find(iSNR_db_values==0);
a=(log10(mu_values(end))-log10(mu_values(idx0)))/iSNR_db_values(end);
iSNR_db0=-log10(mu_values(idx0))/a;
iSNR_db0=round(iSNR_db0*100)/100;
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XTick',sort([-10 0:5:10 iSNR_db0])); 
box on; grid on;
