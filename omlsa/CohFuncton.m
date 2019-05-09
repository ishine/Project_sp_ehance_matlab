function [enhanced_ouput]=ReIm_Coherence(x1,x2,fs,output_fn)

% x1 , x2 : Input signals at two channels
% fs : Sampling Frequency
% output_path: address to write enhanced file
% Nima Yousefian , Mar 2015

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
%frameLength= 240;

if (frameLength/2 ~= floor(frameLength/2))
      frameLength=frameLength+1;
end
 
[x3,fs3  ]= audioread('F:/Work/2018/Beamforming/matlab/GSCLMS/voice/1.wav');
 
 
frameShift=floor(frameLength* 0.5);
window=hanning(frameLength);
FFT_LEN=2^nextpow2(frameLength);

  fid= fopen('F:/Work/2018/MIC_ARRAY_ANC/matlabfile4/Hanning.bin','a+');
  
 fprintf(fid,'%f,\n ' ,window );       
 fclose(fid);
% fprintf('window:\n ' );
% fprintf('%f,',window);     
% fprintf(':\n ');   
      
  
lambda_x=0.68;   %Forgetting factor for smoothing power spectrum 
epsilon=10^-12;
    
lenS=min(length(x3) );
nFrame=0;
iniFrameSample=1;
endFrameSample=iniFrameSample+frameLength-1;
enhanced_ouput=zeros(lenS,1);
enhanced_ouput_fix =zeros(lenS,1,'int16');

TotFrameNum=floor(lenS/frameShift);

cohX = zeros( FFT_LEN, 1); %fftlen
DRR = zeros( FFT_LEN/2, 1);  % fftlen/2

K = zeros(FFT_LEN/2,1);      %fftlen/2
%K =1;
dis= 0.05 ;%0.015;          %Microphones distance (m)
c=340;              %Sound Speed (m/s)
tau=fs * (dis /c);      
Target_DOA= 0 *pi/180;    %DOA -target angle (0 -> Endfire, 90 -> Broadside) 
Beta = 10 *pi/180; 
%setting freq. bins between 0 and pi/2
f_bin=0:1/FFT_LEN:0.5;     
f_bin(1)=[];
omega=2 * pi * f_bin;
omega_tau=omega'*tau.*cos(Target_DOA);  
Beta_tau = omega_tau.* cos(Beta);  
%pri_min=10^(-40/10);
pri_min=10^(-9);

frame_inx = 0;

history = zeros(frameShift,1);

while endFrameSample<lenS
        
        %A new frame to process
        nFrame=nFrame+1;           
        Frame1=x3(iniFrameSample:endFrameSample);      
        wFrame1=Frame1 .* window;
        X1=fft(wFrame1,FFT_LEN);
        
      %IFFT and OLA
      enhSpeech_Frame_tmp=real(ifft(   X1,FFT_LEN));        
     % enhSpeech_Frame_tmp=real(ifft(  X1,FFT_LEN));       
      enhanced_ouput(frame_inx*frameShift +1 :(frame_inx+1) *frameShift) = enhSpeech_Frame_tmp(1:frameShift) + history;      
      enhanced_ouput_fix(frame_inx*frameShift +1 :(frame_inx+1) *frameShift) =    enhanced_ouput(frame_inx*frameShift +1 :(frame_inx+1) *frameShift).*2^15 ;
      history=enhSpeech_Frame_tmp(frameShift+1:FFT_LEN);
             
      %Update frame boundaries
      iniFrameSample=iniFrameSample+frameShift;
      endFrameSample=endFrameSample+frameShift;
      
      frame_inx = frame_inx+1;
end

    disp('Processing ... 100% Done');
   % wavwrite(enhanced_ouput,fs,16,output_fn);
 
     audiowrite('F:/Work/2018/Beamforming/matlab/GSCLMS/voice/1_coh.wav',enhanced_ouput_fix,16000, 'BitsPerSample',16);