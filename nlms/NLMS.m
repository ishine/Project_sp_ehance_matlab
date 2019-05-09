
%  目前该版本为 lms 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

% nlms flow
%1）给定w(0)w(0)；

%2）计算输出值：y(k)=w(k)Tx(k)y(k)=w(k)Tx(k);

%3）计算估计误差：e(k)=d(k)?y(k)e(k)=d(k)?y(k);

%4）权重更新：w(k+1)=w(k)+ (μ / (α+|x(k)|)) * 2x(k)e(k)

close all; clc; clear all;

fin0 = '..\..\voice\t164l';
fin1 = '..\..\voice\t164r';
fout = [fin0 '_NLSM'];

omega = 0.1;
 
Nmic = 2; % mic number 
Nflt = 1024; % fir length
Num_c = Nmic*Nflt;
 A_st = zeros(Nflt, 1);   
 R_st = zeros(Nflt, 1);   
 
[X1,fs1]   = audioread([fin0 '.wav']);  % main mic
 [X2,fs2]   = audioread([fin1 '.wav']); % ref mic


Y_Up   = (X1 ) ;
lenS1 =length(Y_Up);
Y_Down =  X2;
fs = fs1;

lenS2 =length(Y_Down);
lenS = min(lenS1,lenS2); 
 
u = 0.01;
omega = 0;
mu = 0.073;
k = 1;
K = 256;
alpha = 1e-2;

lambda = 0.05;
 
En =zeros(lenS,1); 
Y_LMS =zeros(lenS,1); 
  
w = zeros(Nmic - 1, K);
 
 En(k:k+Nflt-1) =   Y_Up(k:k+Nflt-1) ;
  
 yB = 0;
 Pest =0;
 for k=Nflt:lenS - Nflt+1        
      Y_Frame_Block = Y_Down(k-Nflt+1 :k);      
      yB = ( Y_Frame_Block' * Y_Frame_Block);           
      Y_LMS(k) = A_st'*Y_Frame_Block;   
      err = Y_Up(k) - Y_LMS(k) ;     
      mu = (alpha /yB);   
      A_st = A_st + mu*err*Y_Frame_Block;  
      En(k) = err;         
 end
 audiowrite([fout '_y.wav'],Y_LMS,fs);
 audiowrite([fout '_out.wav'],En,fs); 
 fprintf('nlms end\n');
