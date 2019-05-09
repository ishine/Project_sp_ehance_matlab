
%  目前该版本为 lms 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

% nlms flow
%1）给定w(0)w(0)；

%2）计算输出值：y(k)=w(k)Tx(k)y(k)=w(k)Tx(k);

%3）计算估计误差：e(k)=d(k)?y(k)e(k)=d(k)?y(k);

%4）权重更新：w(k+1)=w(k)+ (μ / (α+|x(k)|)) * 2x(k)e(k)

close all; clc; clear all;

omega = 0.75;
 
Nmic = 2; % mic number 
Nflt = 512; % fir length
Num_c = Nmic*Nflt;
 A_st = zeros(Nflt, 1);   
 
 
[Y_Up,fs1]= audioread('../voice/T7L1_GSC_UP.wav');  % main mic
[Y_Down,fs2]= audioread('../voice/T7L1_GSC_B.wav'); % ref mic
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
 
En =zeros(lenS,1); 
Y_LMS =zeros(lenS,1); 
w = zeros(Nmic - 1, K);
 
 En(k:k+Nflt-1) =   Y_Up(k:k+Nflt-1) ;
  
 yB = 0;
 Pest =0;
 for k=Nflt:lenS - Nflt+1    
    
     Y_Frame_Block = Y_Down(k-Nflt+1 :k);  
     % yB =   omega *  yB  +  (1- omega) *( Y_Frame_Block' * Y_Frame_Block);       
       yB =    ( Y_Frame_Block' * Y_Frame_Block);       
      
      Y_LMS(k) = A_st'*Y_Frame_Block;   
      err = Y_Up(k) - Y_LMS(k) ;
     % --- lms --- %         
     %  update mu
     %  mu = u / (alpha + |X(K)*X(K)|)   
     
     mu = u / (alpha + yB);
     A_st = A_st + mu*err*Y_Frame_Block;  
    
    %   Pest = omega * Pest + (1-omega) * yB;
    %   A_st = A_st + mu*err  * Y_Frame_Block / Pest;   
    
      
     En(k) = err;    
     
 end
   
audiowrite('../voice/NLSM_OUT_t.wav',Y_LMS,fs);
 audiowrite('../voice/NLMSEn_t.wav',En,fs);
 
fprintf('nlms end\n');
