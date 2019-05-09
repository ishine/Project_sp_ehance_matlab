close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;
c=340;

theta_d=0/180*pi;
delta=1e-2;
M=8;

N_values=[1 2 4 8]; 

f=(2000:500:8000)';
D_dB_values=zeros(length(f),length(N_values));
for idxN=1:length(N_values)
    N=N_values(idxN);
    [m_mat,n_mat]=meshgrid(1:M,1:M);
    for idx_f=1:length(f)
        fk=f(idx_f);
        d=exp(-1i*2*pi*fk*delta/c*(0:M-1)'*cos(theta_d));
        Gamma0=sinc(2*fk*delta/c*(m_mat-n_mat));
        
        [T,Lambda]=jeig(d*d',Gamma0);
        [diagL,idx]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx);
        
        T1N=T(:,1:N);
        P1N=T1N/(T1N'*T1N)*T1N';
        
        h=(P1N*d)/(d'*P1N*d);
        D_dB_values(idx_f,idxN)=10*log10(abs(h'*d)^2/(h'*Gamma0*h));
    end
end

figure
plot(f/1e3,D_dB_values(:,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(f/1e3,D_dB_values(:,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D_dB_values(:,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(f/1e3,D_dB_values(:,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;


