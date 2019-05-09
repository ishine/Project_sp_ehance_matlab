close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

A0=0.5;
f0=0.1;
theta0=90*pi/180;   % direction of the desired source
theta1=0*pi/180;   % direction of the interference
alpha=0.1; % ratio between the variance on the white noise and the interference
beta=0.9;
mu=0.5;

M_values=[10 20 50 100]; % Number of sensors
Rx1=A0^2;

iSNR_db_values=-5:15;  
Gain_dB_values=zeros(length(iSNR_db_values),length(M_values));
MSE_dB_values=zeros(length(iSNR_db_values),length(M_values));
nu_dB_values=zeros(length(iSNR_db_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    
    Rvw=alpha*eye(M);
    [t1,t2]=meshgrid(1:M,1:M);
    Rvi=beta.^(abs(t1-t2)*cos(theta1));
    Rv0=Rvi+Rvw;
    d=exp(-1i*2*pi*f0*(0:M-1)'*cos(theta0));
    Rx=Rx1*(d*d');
    i_i=[1; zeros(M-1,1)];
    for idxS=1:length(iSNR_db_values)
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);
        
        sigma_i2=Rx1/(1+alpha)/iSNR;
        Rv=sigma_i2*Rv0;
        Ry=Rx+Rv;
        
        % h=Ry\Rx*i_i;  % Wiener filter
        h=(Rx+mu*Rv)\Rx*i_i;
        
        %iSNR=trace(Rx1)/trace(Rv(1:L,1:L));
        %iSNR_dB=10*log10(iSNR);
        oSNR=(h'*Rx*h)/(h'*Rv*h);
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxM)=oSNR_dB-iSNR_dB;
        
        MSE_dB_values(idxS,idxM)=10*log10( Rx1+h'*Ry*h-h'*Rx*i_i-i_i'*Rx*h );
        
        nu_dB_values(idxS,idxM)=10*log10( (h-i_i)'*Rx*(h-i_i)/Rx1 );
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

figure
plot(iSNR_db_values(1:2:end),nu_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),nu_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),nu_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),nu_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

