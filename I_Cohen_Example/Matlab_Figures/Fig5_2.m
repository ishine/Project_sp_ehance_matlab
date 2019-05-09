close all; clc; clear all;
linewd = 0.8;
hcfontsize = 20;
MarkerSize=9;

f=0.01;
alpha=0.01;

M_values=1:20;
theta_values=[30 50 70 90]/180*pi;
gain=zeros(length(M_values),length(theta_values));
for idx1=1:length(M_values)
    M=M_values(idx1);
    for idx2=1:length(theta_values)
        theta=theta_values(idx2);
        d=exp(-1i*2*pi*f*(0:M-1)'*cos(theta));
        du=exp(-1i*2*pi*f*(0:M-1)');
        gain(idx1,idx2)=10*log10(d'*inv(1/(1+alpha^2)*(du*du')+alpha^2/(1+alpha^2)*eye(M))*d);
    end
end

figure
plot(M_values(1:1:end),gain(1:1:end,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(M_values(1:1:end),gain(1:1:end,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(M_values)]);
set(gca,'YLim',[0 60]);
box on; grid on;

f=0.05;
alpha=0.01;

M_values=1:20;
theta_values=[30 50 70 90]/180*pi;
gain=zeros(length(M_values),length(theta_values));
for idx1=1:length(M_values)
    M=M_values(idx1);
    for idx2=1:length(theta_values)
        theta=theta_values(idx2);
        d=exp(-1i*2*pi*f*(0:M-1)'*cos(theta));
        du=exp(-1i*2*pi*f*(0:M-1)');
        gain(idx1,idx2)=10*log10(d'*inv(1/(1+alpha^2)*(du*du')+alpha^2/(1+alpha^2)*eye(M))*d);
    end
end

figure
plot(M_values(1:1:end),gain(1:1:end,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(M_values(1:1:end),gain(1:1:end,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(M_values)]);
box on; grid on;

f=0.1;
alpha=0.01;

M_values=1:20;
theta_values=[30 50 70 90]/180*pi;
gain=zeros(length(M_values),length(theta_values));
for idx1=1:length(M_values)
    M=M_values(idx1);
    for idx2=1:length(theta_values)
        theta=theta_values(idx2);
        d=exp(-1i*2*pi*f*(0:M-1)'*cos(theta));
        du=exp(-1i*2*pi*f*(0:M-1)');
        gain(idx1,idx2)=10*log10(d'*inv(1/(1+alpha^2)*(du*du')+alpha^2/(1+alpha^2)*eye(M))*d);
    end
end

figure
plot(M_values(1:1:end),gain(1:1:end,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(M_values(1:1:end),gain(1:1:end,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(M_values)]);
box on; grid on;

f=0.2;
alpha=0.01;

M_values=1:20;
theta_values=[30 50 70 90]/180*pi;
gain=zeros(length(M_values),length(theta_values));
for idx1=1:length(M_values)
    M=M_values(idx1);
    for idx2=1:length(theta_values)
        theta=theta_values(idx2);
        d=exp(-1i*2*pi*f*(0:M-1)'*cos(theta));
        du=exp(-1i*2*pi*f*(0:M-1)');
        gain(idx1,idx2)=10*log10(d'*inv(1/(1+alpha^2)*(du*du')+alpha^2/(1+alpha^2)*eye(M))*d);
    end
end

figure
plot(M_values(1:1:end),gain(1:1:end,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(M_values(1:1:end),gain(1:1:end,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(M_values(1:1:end),gain(1:1:end,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(M_values)]);
box on; grid on;


