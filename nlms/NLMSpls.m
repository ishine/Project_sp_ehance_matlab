
%  目前该版本为 lms 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

% nlms flow
%1）给定w(0)w(0)；

%2）计算输出值：y(k)=w(k)Tx(k)y(k)=w(k)Tx(k);

%3）计算估计误差：e(k)=d(k)?y(k)e(k)=d(k)?y(k);

%4）权重更新：w(k+1)=w(k)+ (μ / (α+|x(k)|)) * 2x(k)e(k)

close all; clc; clear all;

fin0 = '..\voice\T16L';
fin1 = '..\voice\T16R';
fout = [fin0 'NLSM_OUT_t'];

omega = 0.75;
 
Nmic = 2; % mic number 
Nflt = 1024; % fir length
Num_c = Nmic*Nflt;
 W = zeros(Nflt, 1);   
 R_st = zeros(Nflt, 1);   
 
[X1,fs1]   = audioread([fin0 '.wav']);  % main mic
[X2,fs2] = audioread([fin1 '.wav']); % ref mic


Y_Up = X1;%(X1 + X2).*0.5;
Y_Down = X2;%X1 - X2;
fs = fs1;
lenS1 =length(Y_Up);
lenS2 =length(Y_Down);
lenS = min(lenS1,lenS2); 
 
u = 0.01;
omega = 0;
mu = 0.073;
k = 1;
K = 256;
alpha =1;% 1e-2;

lambda = 0.05;
 
En =zeros(lenS,1); 
y =zeros(lenS,1); 
  
w = zeros(Nmic - 1, K);
 
 En(k:k+Nflt-1) =   Y_Up(k:k+Nflt-1) ;
  
 yB = 0;
 Pest =0;
 for k=Nflt:lenS - Nflt+1    
    
       x = Y_Down(k-Nflt+1 :k);         
       yB =    ( x' * x);           
       y(k) = W'*x;   
       e = Y_Up(k) - y(k) ;  
       mu = (alpha / (yB+1e-10));   
       W = W + mu*e*x;   
       En(k) = e;         
 end
 audiowrite([fout '_y.wav'],y,fs);
 audiowrite([fout '_e.wav'],En,fs);
 %audiowrite('f:/Work/2018/Beamforming/matlab/GSCLMS/voice/NLMSEn_t.wav',En,fs);
 
fprintf('nlms end\n');
