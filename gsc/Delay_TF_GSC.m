%Delay TF-GSC BEAMFORMING ALGORITHM


fin0 = './voice/T10L';
fin1 = './voice/T10R';
fA   = [fin0 '_TFGSC_A.wav'];
fB   = [fin0 '_TFGSC_B.wav'];
fN   = [fin0 '_TFGSC_N.wav'];
fOut = [fin0 '_TFGSC_OUT'];


[x1,fs1]= audioread([fin0 '.wav']); % main mic
[x2,fs2]= audioread([fin1 '.wav']); % ref mic


fs = fs2;
 
xcut = x1 - x2;

%  
Len_x1 = length(x1);
Len_x2 = length(x2);

x3 = zeros(Len_x1,1);

x3(3:Len_x2) = x2(1:Len_x2-2);

TotalLen = max(Len_x1, Len_x2);
FrameLen = 1024;
FFT_LEN=2^nextpow2(FrameLen);
MicNum = 2;

MArray=0.05*[-2.5;-1.5;-0.5;0.5;1.5;2.5];
alpha = 0.01; 
omega = 0.01;
 
frameShift=floor(FrameLen* 0.5); %  for 50% overlap

X_matrix = [x1,x3]';

N=fix((TotalLen)/frameShift)*frameShift+frameShift;

X_array = [X_matrix, zeros(2,N -TotalLen )];

Y_TFGSC = zeros(N, 1);
Y_TFGSC_B = zeros(N, 1);
Y_TFGSC_NC = zeros(N, 1);
Y_U_TC = zeros(N, 1);
Y_TF_GSC_OUT= zeros(N, 1);
Y_FBF_TC =  zeros(N, 1);

DesAng = 150;
Freband=linspace(0,fs/2,N/2+1);
timeDelay=zeros(MicNum,N);

%% steer vector
for freIdx=1:length(Freband)
    coefficient = 2 * pi * Freband(freIdx)/340.0;

     h=cosd(DesAng); % Ä£Äâ
     h= [cosd(DesAng),sind(DesAng)].';
    tao = MArray*h;
    timeDelay(:,freIdx)=exp(j*coefficient*tao);
end


for m=1:MicNum
    tmpDelay=timeDelay(m,1:(N/2+1)).';
    tmpConjD=timeDelay(m,2:N/2).';
    timeDelay(m,:)=[tmpDelay;flipud(conj(tmpConjD))].';
end


 
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

 Z = zeros(MicNum, FFT_LEN);
 
 YFBF = zeros(FFT_LEN,1);
 Y_NC = zeros(FFT_LEN,1);
 Pest = zeros(FFT_LEN,1); 
 G = zeros(FFT_LEN,1);
 
 
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
       
      % Um = Zm - (Am/A1) * Z1        
      BM_t = BM'; %  ->   FFT_LEN * 2    
     
      for freIdx=1:FFT_LEN     
          U(freIdx) =  BM_t(freIdx,:) * Z(: , freIdx); %  FFT_LEN * MIC_NUM * MICNUM* FFT_LEN
      end
       
     % freq nlms
     Y_NC =  G.* U;     
      
     YOUT = YFBF - Y_NC;     
     Pest = omega * Pest + (1-omega) * abs(U).^2;     
     mu = alpha ./ Pest;    
     G = G + mu .* YOUT .* conj( U)  ;      
     %  end of nlms
  
    
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
 
  audiowrite(fN,Y_TFGSC,fs);
  audiowrite(fB,Y_U_TC,fs);
  audiowrite(fA,Y_FBF_TC,fs);
  audiowrite([fOut '.wav'],Y_TF_GSC_OUT,fs);
  audiowrite('f:/Work/2018/Beamforming/matlab/GSCLMS/voice/xcut.wav',xcut,fs);  

%   omlsa(fOut);
 
  clear all;

  fprintf('End of TF-GSC\n');
