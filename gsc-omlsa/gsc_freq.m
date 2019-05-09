 

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 

Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.08;
Speed = 340;
 
finl = './voice/t10';
  
%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_freq_gsc'];
  
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
 
 N21 =  NFFT/2 + 1;
 
 frameshift = NFFT/2;
 
 Frame_start = 1;
 Frame_end   = NFFT;
 l = 1;
 
 
 d = 0.08;  c = 340;
 angle =  90;
 
 tau = d *cos( pi* angle / 180)/ 340;
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
   his_e=  zeros(frameshift,1); 
 win = hanning(NFFT);
 
 
 x1_gain = x1 * 1/(sqrt(4*pi) * d);
 
 x_1_gain_array  = [ zeros(fix(tauPoint),1) ; x1_gain(1: length(x1_gain ) - fix(tauPoint) ) ];
 
 x_ds_t = [ x_1_gain_array  ,  x2 ];
 
  audiowrite([fout '_delay_orig.wav'],x_ds_t,fs);
   
  NFIR = 16;
  x_buffer_F = zeros(N21,NFIR);
  
  w = zeros(N21, NFIR );e_F = zeros(N21,1);

  pest = zeros(N21,1);
  
  alpha = 0.0075; omega = 0.988;
 
  e_F_full = zeros(NFFT,1);
  
 while(Frame_end < lenS)

     x_1 = x1(Frame_start:Frame_end);
     x_2 = x2(Frame_start:Frame_end);
     
     X_1F = fft(win.*x_1);
     X_2F = fft(win.*x_2);
      
     X_1_F   = X_1F .* phase_shift;%  + X_2F;     
     X_2_F   = X_2F ;%.* phase_shift;%  - X_2F;
      
     
     X_UP_F = (X_1_F + X_2_F)*0.5;
     X_B_F  = X_1_F - X_2_F;  
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
  
 %  x_e_d =  [x_1_;e_t_]';
 
 
 %   audiowrite([fout '_et.wav'],  x_e_d ,   fs);
 
 % x_delay = [x_1_;x_2_]';  
%  audiowrite([fout '_delay.wav'],x_delay,fs);
  
   x_ds = [x_up_; x_b_]';
   audiowrite([fout '_ds.wav'],x_ds,fs);
   
   if(tauPoint>=0)   
     x1 = [zeros(tauPoint,1); x1(tauPoint+1:length(x1))];
   else
     x2 = [zeros(tauPoint,1); x2(tauPoint+1:length(x2))];       
   end
   
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
