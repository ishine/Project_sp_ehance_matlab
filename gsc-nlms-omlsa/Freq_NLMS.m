% freq domain NLMS algorithm
% 频域NLMS, 验证能起到降噪效果, 声音正常， 而且降噪效果要好于时域NLMS  NFFT = 1024时效果较好，低于该值效果差。
 
%fin0 = 'D:\F_Drive\Work\2018\Beamforming\matlab\voice\T16L';
%fin1 = 'D:\F_Drive\Work\2018\Beamforming\matlab\voice\T16R';
fin0 = '.\voice\t10';
 
  
fout = [fin0 '_FreqNLMSEn'];
fout1 = [fin0 '_B'];

omega = 0.85; 
alpha = 0.08;
NFFT = 1024;
OVERLAP = NFFT/2;

window = hanning(NFFT);

[Y_2,fs1]= audioread([fin0 '.wav']); % main mic
%[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic

fs = fs1;
 
 Y_Up   =  (Y_2(:,1) + Y_2(:,2)) /2;
 Y_Down =   Y_2(:,1) - Y_2(:,2);
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

 % omlsa parameter
 NOMLSA = 512; q_OMLSA = NOMLSA/4;
 Y_om_in = zeros(NOMLSA,1);
 U_om_in = zeros(NOMLSA,1);
 yout_omlsa =  zeros(lenS,1);
 k_omlsa = 0;


% Out_Cut = Y_Up + Y_D own;

% audiowrite('f:/Work/2018/Beamforming/matlab/LMS/Out_Cut.wav',Out_Cut,fs);

while (EndFrame < Len)
    
    Frame1 =   Y_Up(InitFrame:EndFrame); % Y_Up(InitFrame:EndFrame);
    Frame2 =   Y_Down(InitFrame:EndFrame); %Y_Down(InitFrame:EndFrame);
    
     X1 = fft(window .* Frame1);%  X1 = fft(window .* Frame1);
     X2 = fft(window .* Frame2);%  X2 = fft(window .* Frame2);
   
     X_lms =  G.* X2;          
     En = X1 - X_lms; 
     Pest = omega * Pest + (1-omega) * abs(X2).^2;      
     mu = alpha ./ Pest;    
     G = G + mu .* En .* conj( X2);
  
     A_up = ifft(X1);
     B_down = ifft(X2);
 
     Nlms = ifft(X_lms);
     Out  = ifft(En);
     
     Out_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = Out(1:OVERLAP) + His_Out;
     His_Out = Out(OVERLAP+1 : NFFT);
     
     B_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) = B_down(1:OVERLAP) + His_B;
     His_B = B_down(OVERLAP+1 : NFFT);
     
     %omlsa 
     out_temp = Out_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP) ;
     b_temp = B_T((FrameNO-1)*OVERLAP+1 : FrameNO*OVERLAP);
      % 1)
      i_p = 0;
      Y_om_in  = [ Y_om_in(q_OMLSA+1:NOMLSA) ; out_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];
      U_om_in  = [ U_om_in(q_OMLSA+1:NOMLSA);  b_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];       
      yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );   
      yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;       
      k_omlsa = k_omlsa+1;        
      % 2)
      i_p = 1;
      Y_om_in  = [ Y_om_in(q_OMLSA+1:NOMLSA) ; out_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];
      U_om_in  = [ U_om_in(q_OMLSA+1:NOMLSA);  b_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];       
      yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );   
      yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;       
      k_omlsa = k_omlsa+1;          
      % 3)
      i_p = 2;
      Y_om_in  = [ Y_om_in(q_OMLSA+1:NOMLSA) ; out_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];
      U_om_in  = [ U_om_in(q_OMLSA+1:NOMLSA);  b_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];      
      yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );   
      yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;       
      k_omlsa = k_omlsa+1;     
      % 4)
      i_p = 3;
      Y_om_in  = [ Y_om_in(q_OMLSA+1:NOMLSA) ; out_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];
      U_om_in  = [ U_om_in(q_OMLSA+1:NOMLSA);  b_temp(i_p*q_OMLSA+1: (i_p+1)*q_OMLSA )];      
      yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );   
      yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;       
      k_omlsa = k_omlsa+1;           
    
       
    FrameNO = FrameNO +1;
    InitFrame = InitFrame + OVERLAP;
    EndFrame  = EndFrame  + OVERLAP;
    
end
 
 Out_B = [Out_T,B_T];
 
 audiowrite([fout '_Out_B.wav'],Out_B,fs);
 audiowrite([fout '_omlsa.wav'],yout_omlsa,fs); 
  fprintf('End of freq nlms\n');
 
 %  gsc_dual_postfilter([fout '_Out'],  [fout '_down']);
 
%  fprintf('End of dual channel omlsa\n');
 

