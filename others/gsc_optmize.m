
%   "A Rboust Adaptive Beamformer with a Blocking Martix Using  Coeffient-Constrained Adaptive Filters", 1999, Osamu JOUSHUYMA,...
%   Optmized GSC Beamformer. Implement in 2018/10/12  yang.yang@nanosic.com

close all; clc; clear all;
  
Nmic = 2; % mic number 
Nflt = 512; % fir length
Num_c = Nmic*Nflt;
  
LOCAL = 1;
if(LOCAL==1)
    
finl = '../voice/T7L.wav';
finr = '../voice/T7L.wav';


[x1,fs1]= audioread();  % main mic
[x2,fs2]= audioread('F:/Work/2018/Beamforming/matlab/voice/T6R.wav'); % ref mic
fs = fs1;

%define the parameter

Dmic = 0.08;
Angle =  90; % angle of signal
Speed = 340;
tau_dis =  cos(Angle*pi/180.0)*Dmic;
tau = tau_dis /Speed;  %delay 

 

taulen = fs* tau;


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
 X_martix = [xOne,xTwo] ;
 
else

% 获取输入信号 做这步纯属是为了方便理解
 xOne=sigArray(:,1);
 xTwo=sigArray(:,2);
 
 lenS = length(xOne);
 
 audiowrite('f:/Work/2018/Beamforming/matlab/GSCLMS/voice/xOne.wav',xOne,fs);
 audiowrite('f:/Work/2018/Beamforming/matlab/GSCLMS/voice/xTwo.wav',xTwo,fs);% ref mic


%X=[x1 x2 x3 x4 x5 x6].'; 
X_martix=[xOne xTwo]; 

end
 
%prepare the init parameter
f_st = zeros(Nflt,1);
f_st(1:4) =  [1,-2,1.5,2];  
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
 
 
 % down line of GSC beamformer
 Y_Block = xOne - xTwo;

 % up line of GSC beamformer 
 Y_FBF  =  (xOne+xTwo) /2; % Y_FBF 
 k=1;
 Yc_k =  Y_FBF(k:k+Nflt-1); 
  
 Y_C = f_st' * Yc_k;
  
 B = zeros(Nmic-1, Nmic);
for i = 1 : Nmic - 1
    B(i, i) = 1;
    B(i, i+1) = -1;
end 
 %mu =0.001;
 
 %% ccaf 
 Q = 5;
 wccaf0 = zeros(Nflt, 1); 
 wccaf1 = zeros(Nflt, 1); 
 wccaf_Beta = 0.1;
 ybuf0 = zeros(Nflt,1) ;
 ybuf1 = zeros(Nflt,1) ; 
 %% laf
 wlafs0= zeros(Nflt, 1); 
 wlafs1= zeros(Nflt, 1); 
 alpha_lafs = 0.2;
 gamma_lafs = 1.0e-6;
 
K = 256;
alpha = 0.1;
 
 
 Y_Out(k:k+Nflt-1) =   Y_FBF(k:k+Nflt-1) ;
  
 for k=Nflt:lenS - Nflt+1    
     
    %% ccaf based BM:  Input: 1) 2 mic data input 2) fbf output Output:     
    x0 = xOne(k-Q );
    x1 = xTwo(k-Q );  
        
    fbf = Y_FBF(k-Nflt+1:k); 
    
    O_d0 =  (sum(fbf'*fbf)); % ||D(K)||^2
    O_d1 =  (sum(fbf'*fbf));      
    
    % get output of wccaf  y(k) = x(k-P) - H*D(K)
    y0 = x0 - wccaf0' * fbf;
    y1 = x1 - wccaf1' * fbf;
    
    % update 2 wccaf coeff
    wccaf0 = wccaf0 + wccaf_Beta * y0 * fbf / O_d0; 
    wccaf1 = wccaf1 + wccaf_Beta * y1 * fbf / O_d1;
    
    ybuf0(1:Nflt-1) = ybuf0(2:Nflt);
    ybuf1(1:Nflt-1) = ybuf1(2:Nflt);    
    ybuf0(Nflt) = y0;
    ybuf1(Nflt) = y1;    
    %% LAFs
    % z(k) = d(k-Q) - sum(W(k) * Y(k))
    % W(k) = [w0, w1 ,...,wL-1]
    % Y(k) = [y(k),y(k-1), ... ,y(k-L+1)] (m = 0,1,...,M-1)
    % W(k+1) =(1 - gama) * W(k) + aplha* z(k) *Y(k) / sum(  sum(y0)^2 +  sum(y1)^2  )
      d = ( wlafs0'* ybuf0 + wlafs1'*ybuf1 );
     z = fbf(Nflt-Q) - d;  
    
     o_y0y1 =  sum(ybuf0'*ybuf0) + sum(ybuf1'*ybuf1);
     
     wlafs0 =  (1-gamma_lafs) * wlafs0 + alpha_lafs * z * ybuf0 /o_y0y1;
     wlafs1 =  (1-gamma_lafs) * wlafs1 + alpha_lafs * z * ybuf1 /o_y0y1;
    %%  output
     Y_Out(k) = z;  
     Y_U(k)   = y0;
     Y_Block(k) = y1;
     Y_D(k) = d;
 end
   
audiowrite('f:/Work/2018/Beamforming/matlab/voice/OPT_GSC_OUT.wav',Y_Out,fs);
audiowrite('f:/Work/2018/Beamforming/matlab/voice/GSC_up.wav',Y_FBF,fs);
audiowrite('f:/Work/2018/Beamforming/matlab/voice/GSC_down.wav',Y_D,fs);
audiowrite('f:/Work/2018/Beamforming/matlab/voice/GSC_B.wav',Y_Block,fs);

fOut = ['f:/Work/2018/Beamforming/matlab/voice/OPT_GSC_OUT'];
%omlsa(fOut);
fprintf('End of GSC\n');
