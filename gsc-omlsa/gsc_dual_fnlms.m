
%  目前该版本为 GSC 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 

Nmic = 2; % mic number 
Nflt = 32; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.08;
Speed = 340;
 
finl = './voice/T11L';
finr = './voice/T11R';


%  142  122  105  90,74 57 37    
Angle =  30; % angle of signal
fout = [finl '_30_fnlms'];
 
 
[x1,fs1]= audioread([finl '.wav']);  % main mic
[x2,fs2]= audioread([finr '.wav']); % ref mic
fs = fs1;
 
fprintf('1st do gsc filter \n');
 
 
tau_dis =  cos(Angle*pi/180.0)*Dmic;
tau = tau_dis /Speed;  %delay 
 
taulen = fs* tau;


tauinx =[-3,-2,-1, 0,1,2,3];
taudelay = tauinx ./ fs * Speed;

tau_angle = acos(taudelay./Dmic)  *180.0/pi;



lenS1 = length(x1);
lenS2 = length(x2);
Len = min(lenS1,lenS2);
%delay one buffer
if(taulen<0)  % delay x1  
   delayPoint = round(-taulen);     
   xTwo = x2(1:Len);  
   xOne =  zeros(Len,1);
   xOne(delayPoint+1:Len) = x1(1:Len-delayPoint);
else          % delay x2
   delayPoint = round(taulen);
   xOne = x1(1:Len);   
   xTwo =  zeros(Len,1);
   xTwo(delayPoint+1:Len) = x2(1:Len-delayPoint);   
end
  
    f_st = zeros(Nflt,1);
    
    f_st(1:4) =  [1,0,0,0];
 
  
f_st = f_st/norm(f_st);
 
C = zeros(Num_c,Nflt);
c = ones(2,1);
 
i=1;j=1;
while i<=Nflt   
    C(j:j+1,i) = c;
    i=i+1;
    j= j+2;
end
 Y_Out = zeros(Len,1);
 Y_U = zeros(Len,1);
 Y_D = zeros(Len,1); 
 A_st = zeros(Nflt, 1); 
 
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
 Y_Up  =  (xOne+xTwo)  /2; 
   
 alpha = 0.0075;
 omega = 0.36;
 NFFT = 1024;
 OVERLAP = NFFT/2;
 frame_s = 1;
 frame_e = NFFT;
 
 win = hanning(NFFT);
 G = zeros(NFFT,1);
 Pest = zeros(NFFT,1);
Nlms = zeros(NFFT,1);
A_T =   zeros(Len,1);
B_T =   zeros(Len,1);
His_A = zeros(OVERLAP,1); 
His_B = zeros(OVERLAP,1);
His_Out = zeros(OVERLAP,1);
FrameNO = 1;
 while(frame_e < Len )
    %% freq nlms
    Frame1 =   Y_Up(frame_s:frame_e); % Y_Up(InitFrame:EndFrame);
    Frame2 =   Y_Block(frame_s:frame_e); %Y_Down(InitFrame:EndFrame);
    
     X1 = fft(win .* Frame1);
     X2 = fft(win .* Frame2);
   
     X_lms =  G.* X2;          
     En = X1 - X_lms; 
     Pest = omega * Pest + (1-omega) * abs(X2).^2;      
     mu = alpha ./ Pest;    
     G = G + mu .* En .* conj( X2);
  
     A_up   = ifft(X1);
     B_down = ifft(X2);
 
     Nlms = ifft(X_lms);
     Out  = ifft(En);

     Out_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = Out(1:OVERLAP) + His_Out;
     His_Out = Out(OVERLAP+1 : NFFT);
     
     B_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = B_down(1:OVERLAP) + His_B;
     His_B = B_down(OVERLAP+1 : NFFT);     
     
    frame_s = frame_s + OVERLAP;    
    frame_e = frame_e + OVERLAP; 
    
    FrameNO = FrameNO +1;
    
 end
  
   
audiowrite([fout '_OUT.wav'],Out_T,fs);
audiowrite([fout '_Bf.wav'],B_T,fs);
audiowrite([fout '_Bt.wav'],Y_Block,fs);
%audiowrite([fout '_B.wav'],Y_Block,fs);

fprintf('2nd do posterfilter \n');
 
 gsc_dual_postfilter(  [fout '_Bf'],[fout '_OUT']);
 
 % omlsa(fOut);
fprintf('End of DUAL GSC\n');
