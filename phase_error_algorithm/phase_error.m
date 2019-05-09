function pahs_error()
 
close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');
%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal

 
finl = './voice/t10l';
finr = './voice/t10r';
fout = [finl '_phase'];
[x1,fs1]= audioread([finl '.wav']);  % main mic
[x2,fs2]= audioread([finr '.wav']); % ref mic
 fs = fs1;

 NFFT = 1024;
 OVERLAP = NFFT/2;
 
frame_start = 1;
frame_end = NFFT;
win = hanning(NFFT);

TOA_tau  = 0;

W = zeros(NFFT,1);
l = 1;
while(frame_end < len)
     
    x1_frame = x1(frame_start:frame_end);    
    x2_frame = x2(frame_start:frame_end);    
    
    X1_TF = fft(x1_frame.*win);
    X2_TF = fft(x2_frame.*win);    
    
    X1_ang = angle(X1_TF/abs(X1_TF)); 
    X2_ang = angle(X2_TF/abs(X2_TF)); 
    
    angle_X  =  X1_ang - X2_ang - TOA_tau;
    
    W = 1./(1+ gamma * angle_X);
    
    X = W .* X1_TF;
    
    x = ifft(X);
    
    x_out( (l-1)*OVERLAP+1: l*OVERLAP  ) = x(1:OVERLAP) + his_x;    
    his_x = x(OVERLAP+1:NFFT);
    
    frame_start = frame_start + OVERLAP;
    frame_end   = frame_end   + OVERLAP;    
end

   audiowrite([fout '_OUT.wav'],e,fs);

