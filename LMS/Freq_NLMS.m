% freq domain NLMS algorithm
% 频域NLMS, 验证能起到降噪效果, 声音正常， 而且降噪效果要好于时域NLMS  NFFT = 1024时效果较好，低于该值效果差。

omega = 0.9;
mu = 0.073;
u = 0.5;
alpha = 0.1;
NFFT = 1024;
OVERLAP = NFFT/2;

window = hanning(NFFT);

[Y_Up,fs1]= audioread('F:/Work/2018/Beamforming/matlab/LMS/GSC_up.wav');  % main mic
[Y_Down,fs2]= audioread('F:/Work/2018/Beamforming/matlab/LMS/GSC_B.wav'); % ref mic

fs = fs1;

L1 = length(Y_Up);
L2 = length(Y_Down);

Len = min(L1,L2);

InitFrame = 1;
EndFrame = InitFrame+ NFFT-1;
Pest = zeros(NFFT,1);
G = zeros(NFFT,1);
Nlms = zeros(NFFT,1);
Nlms_T =   zeros(Len,1);

His_Nlms = zeros(OVERLAP,1);
 
Out = zeros(NFFT,1);
Out_T = zeros(Len,1);
His_Out = zeros(OVERLAP,1);

FrameNO = 1;

%Out_Cut = Y_Up + Y_Down;

% audiowrite('f:/Work/2018/Beamforming/matlab/LMS/Out_Cut.wav',Out_Cut,fs);

while (EndFrame < Len)
    
    Frame1 = Y_Up(InitFrame:EndFrame);
    Frame2 = Y_Down(InitFrame:EndFrame);
    
    X1 = fft(window .* Frame1);
    X2 = fft(window .* Frame2);
        
     X_lms =  G.* X2;
          
     En = X1 - X_lms;
     
  Pest = omega * Pest + (1-omega) * abs(X2).^2;     
    G = G + mu * En .* conj( X2)  ./ Pest;
    
     % Pest =  (alpha + abs(X2).^2);    
    % G = G + mu * En .* conj( X2) ./ Pest ;
     
     
     Nlms = ifft(X_lms);
     Out  = ifft(En);
     
     Out_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = Out(1:OVERLAP) + His_Out;
     His_Out = Out(OVERLAP+1 : NFFT);
     
     Nlms_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = Nlms(1:OVERLAP) + His_Nlms;
     His_Nlms = Nlms(OVERLAP+1 : NFFT);
     
     FrameNO = FrameNO +1;
    InitFrame = InitFrame+OVERLAP;
    EndFrame = EndFrame+ OVERLAP;
    
end

 audiowrite('f:/Work/2018/Beamforming/matlab/LMS/Out_T_1024.wav',Out_T,fs);
 audiowrite('f:/Work/2018/Beamforming/matlab/LMS/Nlms_T_1024.wav',Nlms_T,fs);
 
 fprintf('End of freq nlms\n');
