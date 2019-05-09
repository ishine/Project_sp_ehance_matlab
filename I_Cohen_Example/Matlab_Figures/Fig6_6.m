close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A0=0.5;
f0=0.1;
theta0=90*pi/180;   % direction of the desired source
theta1=0*pi/180;   % direction of the interference
theta2=45*pi/180;   % direction of the interference
alpha=0.1; % ratio between the variance on the white noise and the interference

M_values=[10 20 50 100]; % Number of sensors
Rx1=A0^2;

iSNR_db_values=-5:15;  
Gain_dB_values=zeros(length(iSNR_db_values),length(M_values));
MSE_dB_values=zeros(length(iSNR_db_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    
    Rvw=alpha*eye(M);
    dv1=exp(-1i*2*pi*f0*(0:M-1)'*cos(theta1));
    Rv10=dv1*dv1';
    dv2=exp(-1i*2*pi*f0*(0:M-1)'*cos(theta2));
    Rv20=dv2*dv2';
    
    Rv0=Rv10+Rv20+Rvw;
    d=exp(-1i*2*pi*f0*(0:M-1)'*cos(theta0));
    Rx=Rx1*(d*d');
    i_i=[1; zeros(M-1,1)];
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        sigma_i2=Rx1/Rv0(1,1)/iSNR;
        Rv=sigma_i2*Rv0;
        Ry=Rx+Rv;
        
        [Qx,Lamba_x]=eig(Rx);
        idx=find(diag(Lamba_x)>1e-10);
        Qxp=Qx(:,idx);
        Rv1=sigma_i2*Rv10;
        [Qv1,Lamba_v1]=eig(Rv1);
        idx=find(diag(Lamba_v1)>1e-10);
        Qv1p=Qv1(:,idx);
        
        Cxv1=[Qxp Qv1p];
        i_c=[Qxp'*i_i; zeros(size(Qv1p,2),1)];
        h=Rv\Cxv1/(Cxv1'/Rv*Cxv1)*i_c;
        
        oSNR=(h'*Rx*h)/(h'*Rv*h);
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxM)=oSNR_dB-iSNR_dB;
        
        MSE_dB_values(idxS,idxM)=10*log10( Rx1+h'*Ry*h-h'*Rx*i_i-i_i'*Rx*h );
    end
end

figure
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),Gain_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),MSE_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),MSE_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),MSE_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),MSE_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

