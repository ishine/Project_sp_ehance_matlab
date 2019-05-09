% TF-GSC BEAMFORMING ALGORITHM


[x1,fs1]= audioread('F:/Work/2018/Beamforming/matlab/TF-GSC/S1M1.wav'); % main mic
[x2,fs2]= audioread('F:/Work/2018/Beamforming/matlab/TF-GSC/S1M2.wav'); % ref mic
fs = fs2;



xcut = x1 - x2;

%  
Len_x1 = length(x1);
Len_x2 = length(x2);
TotalLen = max(Len_x1, Len_x2);
FrameLen = 1024;
FFT_LEN=2^nextpow2(FrameLen);
MicNum = 2;
mu = 0.1;
omega = 0.9;
forgottor = 0.85;
 
frameShift=floor(FrameLen* 0.5); %  for 50% overlap

X_matrix = [x1,x2]';

N=fix((TotalLen)/frameShift)*frameShift+frameShift;

X_array = [X_matrix, zeros(2,N -TotalLen )];

Y_TFGSC = zeros(N, 1);
Y_TFGSC_B = zeros(N, 1);
Y_TFGSC_NC = zeros(N, 1);
Y_U_TC = zeros(N, 1);
 Y_TF_GSC_OUT= zeros(N, 1);
 Y_FBF_TC =  zeros(N, 1);
 
window = hanning(FrameLen); % hanning window

% for system identify  H1 H2
H_tr = zeros(MicNum, FFT_LEN);
H_tr2= zeros(MicNum, FFT_LEN);
SEG_NUM = 13;
SEG_LEN = FFT_LEN * SEG_NUM;
 
tr_buf = zeros(SEG_LEN,2);

InitFrame  = 1;
EndFrame = InitFrame + FrameLen-1;
 
FrameNO = 1;
Ex_H_Inx = 0;
H_Inx = 1;

%BM=zeros(MicNum,MicNum-1,FFT_LEN);
BM=zeros(MicNum,FFT_LEN);

U =  zeros(FFT_LEN,1); % two mic, so just one U

W0 = zeros(MicNum, FFT_LEN);

%G = zeros( MicNum-1, FFT_LEN );
%   G_t = G';

G = zeros(FFT_LEN,1);

 Z = zeros(MicNum, FFT_LEN);
 
 YFBF = zeros(FFT_LEN,1);
 Y_NC = zeros(FFT_LEN,1);
 Pest = zeros(FFT_LEN,1); 
 
 history = zeros(frameShift, 1);
 his_FBF = zeros(frameShift, 1);
 his_U = zeros(frameShift, 1);
 historyout= zeros(frameShift, 1);
 
% main loop 
 while(EndFrame<TotalLen)
     
     % Estimate the TF ratio  [ 1, A2/A1 , A3/A1, ... , Am/A1 ]
     % when data is enough, call TF_ESTIMATE
     
     temp = mod(FrameNO * frameShift, SEG_LEN);
     
    if(mod(FrameNO * frameShift, SEG_LEN)~=0)
        H_Inx=fix(FrameNO * frameShift/SEG_LEN)+1; 
    else
        H_Inx=fix(FrameNO * frameShift / SEG_LEN);
    end
     
    if(H_Inx ~= Ex_H_Inx)
     if(H_Inx * SEG_LEN <= TotalLen )         
       tr_buf =  X_array(:, SEG_LEN * (H_Inx-1) + 1: SEG_LEN*(H_Inx) );         
      
     %H_tr  =  TF_Estimation_TD( tr_buf ,FFT_LEN  ,SEG_NUM,window);  
      H_tr  =  TF_Estimation( tr_buf ,FFT_LEN  ,SEG_NUM,window);  
       
      H_tr_T = H_tr';
     
     end      
     Ex_H_Inx = H_Inx;
    end
        
    
   HStd = H_tr';
  
     % W0
   for freIdx=1:FFT_LEN
        
        TDEN = HStd(freIdx,:);
        
        TNUM = norm(HStd(freIdx,:)).^2;
     
        W0(:,freIdx)= HStd(freIdx,:) ./ norm(HStd(freIdx,:)).^2; % MicNum*NFFT   .*timeDelay(:,freIdx).' 
    end 
    
    
     %   BM, Blocking matrix,  
     %   -A2/A1  -> -H2
     %      1    ->  1
     
     BM(1,:) = -H_tr(2,:);
     BM(2,:) = 1;
     
     % do GSC beamforming
    
     Frame = X_array(:,InitFrame:EndFrame);
     
     Z(1,:) = fft(window' .* Frame(1,:));
     Z(2,:) = fft(window' .* Frame(2,:));
     
    
     W0_t = W0';
     % YFBF = W0 * Z;
      for freIdx=1:FFT_LEN
        YFBF(freIdx) = W0_t(freIdx, : )  * Z(: , freIdx);
      end    
       
     %  Um = Zm - (Am/A1) * Z1        
     BM_t = BM'; %  ->   FFT_LEN * 2    
      for freIdx=1:FFT_LEN     
          U(freIdx) =  BM_t(freIdx,:) * Z(: , freIdx); %  FFT_LEN * MIC_NUM * MICNUM* FFT_LEN
      end
 
     Y_NC = G .* U;
 
     YOUT = YFBF - Y_NC;
    
     Pest = forgottor*Pest+(1-forgottor)* abs(U).^2 ;  % norm(Z(2,freIdx)).^2;
     G = G  + mu * (YOUT .* conj(U )) ./Pest; % (M-1)*1  G在不断的增大
        
        
  % Y_NC =  G.* U;
          
   %  YOUT = YFBF - Y_NC;
     
   %  Pest = omega * Pest + (1-omega) * abs(U).^2;
     
    % G = G + mu * YOUT .* conj( U)  ./ Pest;
        
     
 %  G_time=ifft(G.')';
   %G_time=G_time(:,1:251); % 对应论文中的251个系数的FIR滤波器(不过效果并不是很大)
  % G=fft(G_time.',FFT_LEN).'; 
    
   ytout =  (ifft(YOUT)); 
   Y_TF_GSC_OUT((FrameNO-1)*frameShift+1 : FrameNO*frameShift) =  historyout + ytout(1:frameShift);       
   historyout  = ytout(frameShift+1:FFT_LEN);   
        
   ytnCout =  (ifft(Y_NC)); 
   Y_TFGSC((FrameNO-1)*frameShift+1 : FrameNO*frameShift) =  history + ytnCout(1:frameShift);       
   history  = ytnCout(frameShift+1:FFT_LEN);
 
   yfbfout = ifft(YFBF);   
   Y_FBF_TC((FrameNO-1)*frameShift+1 : FrameNO*frameShift) = yfbfout(1:frameShift) + his_FBF;   
   his_FBF  = yfbfout(frameShift+1:FFT_LEN);
       
    yuout = ifft(U);
    Y_U_TC((FrameNO-1)*frameShift+1 : FrameNO*frameShift) = yuout(1:frameShift) + his_U;   
    his_U  = yuout(frameShift+1:FFT_LEN); 
    
   FrameNO = FrameNO +1;            
    
   InitFrame = InitFrame + frameShift;
   EndFrame  = EndFrame + frameShift;
 end
 
  audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/TF_GSC_NC.wav',Y_TFGSC,fs);
  audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/TFGSC_U.wav',Y_U_TC,fs);
  audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/TFGSC_FBF.wav',Y_FBF_TC,fs);
  audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/Y_TF_GSC_OUT.wav',Y_TF_GSC_OUT,fs);
audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/xcut.wav',xcut,fs);

  
%  audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/.wav',Y,fs);
%  audiowrite('f:/Work/2018/Beamforming/matlab/TF-GSC/LCMVOUT.wav',Y,fs);
 
  
  clear all;

  fprintf('End of TF-GSC\n');
