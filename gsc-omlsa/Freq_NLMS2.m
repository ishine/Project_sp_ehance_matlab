% freq domain NLMS algorithm
% 频域NLMS, 验证能起到降噪效果, 声音正常， 而且降噪效果要好于时域NLMS  NFFT = 1024时效果较好，低于该值效果差。
 
%fin0 = 'D:\F_Drive\Work\2018\Beamforming\matlab\voice\T16L';
%fin1 = 'D:\F_Drive\Work\2018\Beamforming\matlab\voice\T16R';
fin0 = '.\voice\nlms_d';
fin1 = '.\voice\nlms_x';
  
fout = [fin0 '_FreqNLMSEn'];
fout1 = [fin0 '_B'];

omega = 0.31; 
alpha = 0.098;
NFFT = 256;
OVERLAP = NFFT/2;

 

window = hanning(NFFT);

[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic

fs = fs1;
 
 Y_Up   =  Y_L ;%(Y_L + Y_R)./2;
 Y_Down =  Y_R;%Y_L - Y_R;
 
d   =  Y_L ;%(Y_L + Y_R)./2;
 
 x  =  Y_R;  %Y_L - Y_R;
blms = dsp.BlockLMSFilter(128,8);
 [yout, err] =  blms(x,d); 
   audiowrite([fout '_e.wav'],err,fs);
   audiowrite([fout '_yout.wav'],yout,fs);

    gsc_dual_postfilter(    [fout '_e'], fin1);
  %   gsc_dual_postfilterInv(   [fin1],[fout '_e']);
   
 %  gsc_dual_postfilterInv(   [fin1],[fout '_e']);
 
%Y_Up = Y_L;
%Y_Down = Y_R;
% audiowrite([fout '_down.wav'],Y_Down,fs); 
%  audiowrite([fout '_up.wav'],Y_Up,fs);
L1 = length(Y_Up);
L2 = length(Y_Down);

Len = min(L1,L2);

InitFrame = 1;
EndFrame = InitFrame+ NFFT-1;
Pest = zeros(NFFT,1);
G = zeros(NFFT,1);
Nlms = zeros(NFFT,1);
B_T =   zeros(Len,1);
A_T =   zeros(Len,1);
B_T =   zeros(Len,1);
His_A = zeros(OVERLAP,1);
His_B = zeros(OVERLAP,1);
His_B = zeros(OVERLAP,1);
 
Out = zeros(NFFT,1);
Out_T = zeros(Len,1);
His_Out = zeros(OVERLAP,1);

FrameNO = 1;
W_hat = zeros(NFFT,1);

gamma = 0.958;
P =zeros(NFFT,1);
% Out_Cut = Y_Up + Y_D own;
ek_extend = zeros(NFFT,1);
% audiowrite('f:/Work/2018/Beamforming/matlab/LMS/Out_Cut.wav',Out_Cut,fs);

while (EndFrame < Len)
    
    d =   Y_Up(InitFrame:EndFrame); % Y_Up(InitFrame:EndFrame);
    u =   Y_Down(InitFrame:EndFrame); %Y_Down(InitFrame:EndFrame);
    
    % filter
    U = (fft(u));
    Uk = diag(U);
    
    y = ifft(Uk * W_hat);
    yk = y(OVERLAP+1:NFFT);
    
    % error estimate    
    dk = d(OVERLAP+1:NFFT);    
    ek = dk - yk;    
    
    ek_extend(1:OVERLAP) =  zeros(OVERLAP,1 );
    ek_extend(OVERLAP+1:NFFT) = ek;
    
    Ek = fft(ek_extend);
    
    % Power     
    P =gamma * P + (1-gamma)*abs(U).^2;   
    P_inv = 1./(P+1e-10);
    
    Dk = diag(P_inv);
    
    eta = ifft(Dk *  Uk' * Ek);
    etak = eta(1:OVERLAP);
    
    mu = zeros(NFFT,1);
    
    etak_extend = zeros(NFFT,1);
    
    etak_extend(1:OVERLAP) = etak;
    etak_extend(OVERLAP+ 1 : NFFT)=zeros(OVERLAP,1);
    
    mu =  fft(etak_extend);
    
    W_hat = W_hat + alpha *mu; 
     
     Out_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = ek(1:OVERLAP)  ;
   %  His_Out = Out(OVERLAP+1 : NFFT);
     
       
    FrameNO = FrameNO +1;
    InitFrame = InitFrame + OVERLAP;
    EndFrame  = EndFrame  + OVERLAP;
    
end




 audiowrite([fout '_Out.wav'],Out_T,fs);



  fprintf('End of freq nlms\n');
 
 %  gsc_dual_postfilter([fout '_Out'],  [fout '_down']);
 
%  fprintf('End of dual channel omlsa\n');
 

