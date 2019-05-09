
%  目前该版本为 GSC 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 

Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.02;
Speed = 340;
 
finl = './voice/t10';
  
%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_90_gsc'];
 
 
[x,fs1]= audioread([finl '.wav']);  % main mic
fs = fs1;

x1 = x( :, 1);
x2 = x( :, 2); 
 
 
fprintf('1st do gsc filter \n');
  
tau_dis =  cos(Angle*pi/180.0)*Dmic;
tau = tau_dis /Speed;  %delay 
 
taulen = fs* tau;

tauinx =[-3,-2,-1, 0,1,2,3];
taudelay = tauinx ./ fs * Speed;

tau_angle = acos(taudelay./Dmic)  *180.0/pi;

lenS1 = length(x1);
lenS2 = length(x2);
lenS = min(lenS1,lenS2);
  
 NFFT = 512;
 
 frameshift = NFFT/2;
 
 Frame_start = 1;
 Frame_end   = NFFT;
 l = 1;
 
 
 d = 0.08;  c = 340;
 angle =  0;
 
 tau = d *cos(angle)/ 340;
 tauFloat = tau * 16000;
 tauPoint =  fix( tauFloat);
 
 
 n = (1:NFFT)';  % phase delay =  e(-jwt) = e(-j n*2*pi/N  * tau) ; tau = d*cos(angle)/c
 
 pi = 3.1415926;
 
 phase_shiftRe = cos(((n-1)*2*pi *(tauPoint) /NFFT  ) );
 
 phase_shiftIm = sin(((n-1)*2*pi *(tauPoint ) /NFFT  ));
 
 phase_shift =  phase_shiftRe - phase_shiftIm * 1i;
 
 phase_shift2 = exp(-1i*((n-1)*2*pi/NFFT) );
 

 his_up = zeros(frameshift,1);
 his_down =  zeros(frameshift,1);
 his_x_1 = zeros(frameshift,1);
 his_x_2 =  zeros(frameshift,1); 
 
 win = hanning(NFFT);
 
 
 x1_gain = x1 * 1/(sqrt(4*pi) * d);
 
 x_1_gain_array  = [ zeros(fix(tauPoint),1) ; x1_gain(1: length(x1_gain ) - fix(tauPoint) ) ];
 
 x_ds_t = [ x_1_gain_array  ,  x2 ];
 
  audiowrite([fout '_delay_orig.wav'],x_ds_t,fs);
  
  
  
  
 
 while(Frame_end < lenS)

     x_1 = x1(Frame_start:Frame_end);
     x_2 = x2(Frame_start:Frame_end);
     
     X_1F = fft(win.*x_1);
     X_2F = fft(win.*x_2);
     
     
     X_1_F   = X_1F .* phase_shift;%  + X_2F;     
     X_2_F   = X_2F ;%.* phase_shift;%  - X_2F;
       
     x_1_t = real(ifft(X_1_F));
     x_2_t = real(ifft(X_2_F));
       
     x_1_((l-1)*frameshift +1 : l*frameshift) = x_1_t(1:frameshift) + his_x_1;
     x_2_((l-1)*frameshift +1 : l*frameshift) = x_2_t(1:frameshift) + his_x_2;

     his_x_1 = x_1_t(frameshift+1 : NFFT);
     his_x_2 = x_2_t(frameshift+1 : NFFT);
     
     
     X_UP_F = (X_1_F + X_2_F)*0.5;
     X_B_F  = X_1_F - X_2_F;
     
     
     %%  use 'X_B_F' to estimate 'X_UP_F' with W
     
     
     
     
     
     %%
     
     
     x_up_t = real(ifft(X_UP_F));
     x_b_t = real(ifft(X_B_F));     
     
     x_up_((l-1)*frameshift +1 : l*frameshift)  = x_up_t(1:frameshift) + his_up;
     x_b_((l-1)*frameshift +1 : l*frameshift)   = x_b_t(1:frameshift)  + his_down;

     his_up   = x_up_t(frameshift+1 : NFFT);
     his_down = x_b_t(frameshift+1 : NFFT);
    
     
     
     l = l+1;
     Frame_start = Frame_start + frameshift;
     Frame_end   = Frame_end   + frameshift;
     
 end
  
  x_delay = [x_1_;x_2_]';  
  audiowrite([fout '_delay.wav'],x_delay,fs);
  
   x_ds = [x_up_; x_b_]';
   audiowrite([fout '_ds.wav'],x_ds,fs);
   
   x_add = 0.5*(x1 + x2 );
   
   x_sub = (x1 - x2);
   
   x_ds_t = [x_add ,  x_sub] ;
    
 audiowrite([fout '_ds_time.wav'],x_ds_t,fs);
 
%fprintf('2nd do posterfilter \n');
 
 % gsc_dual_postfilter([fout '_OUT'],  [fout '_B']);
 % gsc_dual_postfilter(   [fout '_B'], [fout '_OUT']);
 %gsc_dual_postfilterInv(   [fout '_B'], [fout '_OUT']);
 
%fprintf('End of DUAL GSC\n');
fprintf('End of freq delay\n');
