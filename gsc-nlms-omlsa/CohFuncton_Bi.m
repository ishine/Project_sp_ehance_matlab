function [enhanced_ouput]=CohFuncton_Bi(x1,x2,fs,output_fn)

% x1 , x2 : Input signals at two channels
% fs : Sampling Frequency
% output_path: address to write enhanced file
% Nima Yousefian , Sep 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: N. Yousefian and P. C. Loizou, "A Dual-Microphone Algorithm 
% for Competing Talker Scenarios", IEEE Transactions on Audio, Speech, 
% and Language Processing, Vol. 21, No. 1, pp. 145-155, Jan. 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms labelpf;
 % syms labelpf2
disp('Run Speech Enhacement ReIm');
fs = 16000;

frameLength=floor(15*fs/1000);  
%frameLength=240;

if (frameLength/2 ~= floor(frameLength/2))
      frameLength=frameLength+1;
end
 

fin0 = '..\..\voice\t20';
%fin1 = '.\voice\t12r';
 
[x,fs3]= audioread([fin0 '.wav']); % main mic
%[x4,fs4]= audioread([fin1 '.wav']); % ref mic

x3 = x(:,1);
x4 = x(:,2);
frameShift=floor(frameLength* 0.5);
window=hanning(frameLength);
FFT_LEN=2^nextpow2(frameLength);
 
lambda_x=0.68;   %Forgetting factor for smoothing power spectrum 
epsilon=10^-12;
    
lenS=min(length(x3),length(x4));
nFrame=0;
iniFrameSample=1;
endFrameSample=iniFrameSample+frameLength-1;
enhanced_ouput=zeros(lenS,1);
TotFrameNum=floor(lenS/frameShift);

  
dis= 0.095 ;%0.015;          %Microphones distance (m)
c=340;              %Sound Speed (m/s)
tau=fs * (dis /c);      
Target_DOA= 90 *pi/180;    %DOA -target angle (0 -> Endfire, 90 -> Broadside) 
     
%setting freq. bins between 0 and pi/2
f_bin=0:1/FFT_LEN:0.5;     
f_bin(1)=[];
omega=2 * pi * f_bin;
omega_tau=omega'*tau.*cos(Target_DOA);  

pri_min=10^(-9);
  
while endFrameSample<lenS
        
        %A new frame to process
        nFrame=nFrame+1;
      
 
         Frame1=x3(iniFrameSample:endFrameSample);
         Frame2=x4(iniFrameSample:endFrameSample);

        wFrame1=Frame1 .* window;
        wFrame2=Frame2 .* window;
        
        XL=fft(wFrame1,FFT_LEN);
        XR=fft(wFrame2,FFT_LEN);
        
        
        if (nFrame==1)
            PXLXL=abs(XL).^2;
            PXRXR=abs(XR).^2;
            PXRXL=XL.*conj(XR); 

        else
            %auto and cross PSD estimation
            PXLXL=lambda_x.*PXLXL+(1-lambda_x).*abs(XL).^2;
            PXRXR=lambda_x.*PXRXR+(1-lambda_x).*abs(XR).^2;
            PXRXL=lambda_x.*PXRXL+(1-lambda_x).*XR.*conj(XL);            
        end
        
        cohX=PXRXL ./ sqrt(PXLXL.*PXRXR+epsilon); 
        reCOH=real(cohX(2:FFT_LEN/2+1));
        imCOH=imag(cohX(2:FFT_LEN/2+1));
        
      %%  estimate SNR and G
      reCOH_2 = reCOH.^2;
      imCOH_2 = imCOH.^2; 
           
      W_cur = sqrt((1-reCOH_2-imCOH_2) ./ (2 - 2.* reCOH));
    %   SNR_hat =  (1-reCOH_2-imCOH_2) ./ ((1 - reCOH).^2 + imCOH_2  );        
   %     W_cur =   ((SNR_hat./(1+ SNR_hat)));        
      
    %  W_cur  =   sqrt((SNR_hat./(1+ SNR_hat)));     
      G= sqrt(W_cur);
      H=abs([G ;flipud(G)]);    %Fullband final filter   
    
       %IFFT and OLA
       enhSpeech_Frame_tmp=real(ifft( H .* XL,FFT_LEN));        
       enhSpeech_Frame=enhSpeech_Frame_tmp(1:frameLength);
       enhanced_ouput(iniFrameSample:endFrameSample)=enhSpeech_Frame + enhanced_ouput(iniFrameSample:endFrameSample);      
   
        %Update frame boundaries
        iniFrameSample=iniFrameSample+frameShift;
        endFrameSample=endFrameSample+frameShift;
end


  fout  = [fin0 'coh_bi'];

 audiowrite([fout '_Out.wav'],enhanced_ouput,16000);
  disp('Processing ... 100% Done');

    