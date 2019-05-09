Fs1 = 16000; % Units = Hz
Spec = 'Filter order and transition width';
Order = 64;
TW = 4.0e2;%4.1e2; % Units = Hz
NFFT = 1024;
OVERLAP = NFFT/2;
 framelen = NFFT
 
% Construct FIR Halfband Interpolator
halfbandInterpolator = dsp.FIRHalfbandInterpolator( ...
    'Specification',Spec, ...
    'FilterOrder',Order, ...
    'TransitionWidth',TW, ...
    'SampleRate',Fs1/2, ...
    'FilterBankInputPort',true);

% Construct FIR Halfband Decimator
halfbandDecimator = dsp.FIRHalfbandDecimator( ...
    'Specification',Spec, ...
    'FilterOrder',Order, ...
    'TransitionWidth',TW, ...
    'SampleRate',Fs1);

% Input
fin0 = '..\voice\T7L1';
fin1 = '..\voice\T7R1';
  
fout = [fin0 '_subband'];
%fout1 = [fin0 'B'];

[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic


hlflen = framelen/2;

len = length(Y_L);
Num = len /framelen;
l = 1;
Output =zeros(len,1);
Output =zeros(len,1);
Outbuf = zeros(framelen,1);

SigBuf = zeros(framelen,1);
RefBuf = zeros(framelen,1);
SigBufH = zeros(framelen,1);
RefBufH = zeros(framelen,1);
 
 while(l < Num)
   Y_Buf= Y_L((l-1)*framelen+1 :(l)*framelen );     
   Y_BufR= Y_R((l-1)*framelen+1 :(l)*framelen );  
   [Lowpass,Highpass] = halfbandDecimator(Y_Buf);
   [LowpassR,HighpassR] = halfbandDecimator(Y_BufR);
   % freq nlms    
   SigBuf(1:hlflen)=SigBuf(hlflen+1:framelen);  SigBuf(hlflen+1:framelen) = Lowpass;      
   RefBuf(1:hlflen)=RefBuf(hlflen+1:framelen);  RefBuf(hlflen+1:framelen) = LowpassR;
   NlmsOut=Freq_NLMS_call(SigBuf,RefBuf,l ,1,NFFT);
   
   SigBufH(1:hlflen)=SigBufH(hlflen+1:framelen);  SigBufH(hlflen+1:framelen) = Highpass;
   RefBufH(1:hlflen)=RefBufH(hlflen+1:framelen);  RefBufH(hlflen+1:framelen) = HighpassR;          
   NlmsOutH=Freq_NLMS_call(SigBufH,RefBufH,l,2,NFFT);
   
   Lowpass  =  NlmsOut(1:hlflen);
    Highpass = NlmsOutH(1:hlflen);
   Outbuf = halfbandInterpolator(Lowpass,Highpass);         
   Output((l-1)*framelen +1 :(l)*framelen  )  = Outbuf ;
   l = l +1;   
 end 
 
 
   audiowrite([fout '.wav'],Output,fs1);
 
  