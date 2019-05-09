
%  目前该版本为 lms 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

% nlms flow
%1）给定w(0)w(0)；

%2）计算输出值：y(k)=w(k)Tx(k)y(k)=w(k)Tx(k);

%3）计算估计误差：e(k)=d(k)?y(k)e(k)=d(k)?y(k);

%4）权重更新：w(k+1)=w(k)+ (μ / (α+|x(k)|)) * 2x(k)e(k)

close all; clc; clear all;
fin0 = 'D:\F_Drive\Work\2018\Beamforming\matlab\voice\T15L';
fin1 = 'D:\F_Drive\Work\2018\Beamforming\matlab\voice\T15R';
fout = [fin0 'rls'];
fout1 = [fin0 'B'];

omega = 0.75;
 
Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
 A_st = zeros(Nflt, 1);   
 R_st = zeros(Nflt, 1);   
 
[Y_Up,fs1]= audioread([fin0 '.wav']);  % main mic
[Y_Down,fs2]= audioread([fin1 '.wav']); % ref mic
fs = fs1;
lenS1 =length(Y_Up);
lenS2 =length(Y_Down);
 lenS = min(lenS1,lenS2); 
 
u = 0.5;
omega = 0;
mu = 0.073;
k = 1;
K = 256;
alpha = 0.1;

lambda = 0.05;
 
En =zeros(lenS,1); 
Y_LMS =zeros(lenS,1); 
Y_RLS = zeros(lenS,1);


w = zeros(Nmic - 1, K);
 
 En(k:k+Nflt-1) =   Y_Up(k:k+Nflt-1) ;
  
 yB = 0;
 Pest =0;
 
 % rls

%--------------------------------------------------------------------------
% Filtering
%--------------------------------------------------------------------------

% Filter Parameters
p       = 128;                % filter order
lambda  = 0.99;              % forgetting factor
laminv  = 1/lambda;
delta   = 1.0;              % initialization parameter

% Filter Initialization
w       = zeros(p,1);       % filter coefficients
P       = delta*eye(p);     % inverse correlation matrix
e       = Y_Down*0;              % error signal

for m = p:lenS
    
    % Acquire chunk of data
    x = Y_Down(m:-1:m-p+1);
    
   %  x_2 = Y_Down(m-Nflt+1 :m);  
    
    g = (laminv * P * x) / (1 + laminv * x'*P*x);
            
    d = Y_Up(m);
    
    e = d * w' * x;
    
    w = w   + g * e;
     
    P = laminv * P - laminv * g * x' * P;
    
    
    Y_RLS(m) = e;
end
    
    %rls 
 
 
   

 
 audiowrite([fout '.wav'],Y_RLS,fs); 
 
%audiowrite('f:/Work/2018/Beamforming/matlab/LMS/GSC_down.wav',Y_Down,fs);
%audiowrite('f:/Work/2018/Beamforming/matlab/LMS/GSC_B.wav',Y_Block,fs);

fprintf('rls end\n');
