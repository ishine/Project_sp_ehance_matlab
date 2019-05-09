function pahs_error()
 
close all; clc; clear all;

fprintf('phase error base  \n');
%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal

 
finl = './voice/t10l';
finr = './voice/t10r';
fout = [finl '_phase'];
[x1,fs1]= audioread([finl '.wav']);  % main mic
[x2,fs2]= audioread([finr '.wav']); % ref mic
 fs = fs1;
 
 len = min( length(x1), length(x2)  );

 NFFT = 1024;
 OVERLAP = NFFT/2;
 his_x1 = zeros(OVERLAP,1);  x_out1 = zeros(len,1);
 his_x2 = zeros(OVERLAP,1);  x_out2 = zeros(len,1);  
 x_out = zeros(len,1); 
frame_start = 1;
frame_end = NFFT;
win = hanning(NFFT);

TOA_tau  = 0;
gamma = 2.5;
W1 = zeros(NFFT,1); W2 = zeros(NFFT,1);  angle_X1_ =   zeros(NFFT,1);  angle_X2_ =   zeros(NFFT,1); 
l = 1;
while(frame_end < len)
     
    x1_frame = x1(frame_start:frame_end);    
    x2_frame = x2(frame_start:frame_end);    
    
    X1_TF = fft(x1_frame.*win);
    X2_TF = fft(x2_frame.*win);    
    
    X1_ang = angle((X1_TF )); 
    X2_ang = angle((X2_TF )); 
    
    angle_X1  =  X1_ang - X2_ang - TOA_tau;    
    angle_X1_2 = angle_X1.* angle_X1;        
  %  angle_X1_  = angle_X1_ * 0.38  + angle_X1_2 *  (1-0.38);    
    W1 = 1./(1+ gamma * angle_X1_2);   
    X1 = W1 .* X1_TF;
    x1_phase = ifft(X1);
    
    angle_X2  =  X2_ang - X1_ang - TOA_tau;    
    angle_X2_2 = angle_X2.* angle_X2;        
  %  angle_X2_  = angle_X2_ * 0.38  + angle_X2_2 *  (1-0.38);    
    W2 = 1./(1+ gamma * angle_X2_2);      
    X2 = W2 .* X2_TF;    
    x2_phase = ifft(X2); 
    
    x_out1((l-1)*OVERLAP+1: l*OVERLAP  ) = x1_phase(1:OVERLAP) + his_x1;    
    his_x1 = x1_phase(OVERLAP+1:NFFT);
   
    x_out2((l-1)*OVERLAP+1: l*OVERLAP  ) = x2_phase(1:OVERLAP) + his_x2;    
    his_x2 = x2_phase(OVERLAP+1:NFFT);   
    
    angle_2000(l) = abs(angle_X1(129));
    
    l = l +1;    
    frame_start = frame_start + OVERLAP;
    frame_end   = frame_end   + OVERLAP;    
end
   
       

    x_out = x_out1 + x_out2;

    audiowrite([fout '_OUT1.wav'],x_out1,fs);
    audiowrite([fout '_OUT2.wav'],x_out2,fs);
    audiowrite([fout '_OUT.wav'],x_out,fs);
    
     
    t = 1:length(angle_2000);
    
    plot(t, angle_2000);
    hold on;
    
    
    
    