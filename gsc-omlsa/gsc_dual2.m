
%  目前该版本为 GSC 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');


global SSP  ;
SSP = 0;
%    -3  -2   -1  0  1  2  3 

Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.02;
Speed = 340;
 
finl = './voice/t10';
%finr = './voice/t14r';


%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_90_gsc2'];
 
 
[x,fs1]= audioread([finl '.wav']);  % main mic
%[x2,fs2]= audioread([finr '.wav']); % ref mic

x1 = x(:,1);
x2 = x(:,2);

fs = fs1;


fprintf('1st do gsc filter \n');
%x3 = x1 + x2;

%x4 = x1-x2;

%define the parameter


tau_dis =  cos(Angle*pi/180.0)*Dmic;
tau = tau_dis /Speed;  %delay 
 
taulen = fs* tau;

tauinx =[-3,-2,-1, 0,1,2,3];
taudelay = tauinx ./ fs * Speed;

tau_angle = acos(taudelay./Dmic)  *180.0/pi;

lenS1 = length(x1);
lenS2 = length(x2);
lenS = min(lenS1,lenS2);
%delay one buffer
if(taulen<0)  % delay x1  
   delayPoint = round(-taulen);     
   xTwo = x2(1:lenS);  
   xOne =  zeros(lenS,1);
   xOne(delayPoint+1:lenS) = x1(1:lenS-delayPoint);
else          % delay x2
   delayPoint = round(taulen);
   xOne = x1(1:lenS);   
   xTwo =  zeros(lenS,1);
   xTwo(delayPoint+1:lenS) = x2(1:lenS-delayPoint);   
end
 
 
 
%prepare the init parameter

    f_st = zeros(Nflt,1);
    
    f_st(1:4) =  [1,0,0,0];
    
   % f_st = [1,-2,1.5,2, 0,0,0, 0,   0,0,0, 0,   0,0,0, 0  ]';   
   % f_st = [1,-2,1.5,2, 0,0,0, 0,   0,0,0, 0,   0,0,0, 0  ]';   
   % f_st = zeros(Nflt,1); f_st(1) = 1;
  
f_st = f_st/norm(f_st);
 
C = zeros(Num_c,Nflt);
c = ones(2,1);
 
i=1;j=1;
while i<=Nflt   
    C(j:j+1,i) = c;
    i=i+1;
    j= j+2;
end
 Y_Out = zeros(lenS,1);
 Y_U = zeros(lenS,1);
 Y_D = zeros(lenS,1); 
 A_st = zeros(Nflt, 1); 
 %define the blocking matrix
 Ws = zeros(Nmic-1, Nmic);
 
 c(2) = -1;
 c= c';
 i=1;j=1;
 while i<=Nmic-1   
    Ws(i,j:j+1) = c;
    i=i+1;
    j= j+1;
 end
  
 Pest = 0;
 % down line of GSC beamformer
 Y_Block =    xOne - xTwo;

 % up line of GSC beamformer 
 Y_Up  =   (xOne+xTwo)  /2; % Y_Up 
 
 %Y_Up =  Y_Up - Y_Block;
 
 k=1;
 Yc_k =  Y_Up(k:k+Nflt-1); 
  
 Y_C = f_st' * Yc_k;
 
  
 B = zeros(Nmic-1, Nmic);
for i = 1 : Nmic - 1
    B(i, i) = 1;
    B(i, i+1) = -1;
end 
 %mu =0.001;
 
K = 256;
alpha = 0.0072;
  
 Y_Out(k:k+Nflt-1) =   Y_Up(k:k+Nflt-1) ;
  
 for k=Nflt:lenS - Nflt+1    
    % down line     
     Y_Frame_Block = Y_Block(k-Nflt+1 :k);      
     yB =  Y_Frame_Block' * Y_Frame_Block;             
     Y_Down = A_st'*Y_Frame_Block;   
     
     Y_Frame_Up = Y_Up(k-Nflt+1 :k);
     
   %  yB = Y_Frame_Up'*Y_Frame_Up;
     
     % en = yk - y'k
     yout = Y_Up(k) - Y_Down ;      
 
     %% improved nlms
     omega = 0.898; 
     Pest = omega * Pest + (1-omega) * yB;      
        
     mu = alpha ./ max(Pest,1e-10);        
      
    if(SSP<0.3 ) 
        A_st = A_st + mu*yout*Y_Frame_Block;  
    else        
        A_st = zeros(length(A_st),1);
    end
     
     Y_Out(k) = yout;    
     Y_U(k) = Y_Up(k) ;
     Y_D(k) = Y_Down;
     
 end
  
   
 Y_cmp = [Y_Up , Y_Block];
 Y_cmp2 = [Y_Out, Y_Block ];
audiowrite([fout '_OUT.wav'],Y_Out,fs);
%audiowrite([fout '_UP.wav'],Y_Up,fs);
%audiowrite([fout '_DOWN.wav'],Y_D,fs);
 audiowrite([fout '_B.wav'],Y_Block,fs);
 audiowrite([fout '_cmp.wav'],Y_cmp,fs);
 audiowrite([fout '_cmp2.wav'],Y_cmp2,fs); 
fprintf('2nd do posterfilter \n');
 
  gsc_dual_postfilterV2( Y_Up  , Y_Block  ,fs , [fout '_Y_Up']);
  gsc_dual_postfilterV2( Y_Out  ,  Y_Block ,fs, [fout '_Y_Out']);
 %gsc_dual_postfilterInv(   [fout '_B'], [fout '_OUT']);
 % omlsa(fOut);
fprintf('End of DUAL GSC\n');
