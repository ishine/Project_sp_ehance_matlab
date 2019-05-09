Fs1 = 16000; % Units = Hz
Spec = 'Filter order and transition width';
Order = 64;
TW = 4.0e2;%4.1e2; % Units = Hz
NFFT = 256;
OVERLAP = NFFT/2;
 framelen = NFFT;
 
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
fin0 = '..\voice\T164l';
fin1 = '..\voice\T164r';
  
fout = [fin0 '_subband'];
%fout1 = [fin0 'B'];

[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic
omega = 0.075;

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
   [Lowpass,Highpass] = halfbandDecimator(Y_L);
   [LowpassR,HighpassR] = halfbandDecimator(Y_R); 
   
   len = length(Lowpass);
   En = zeros(len,1);
   EnH = zeros(len,1); 
  k = 1;
  Nflt = 1024; % fir length
  alpha = 0.13;  omega = 0.36;
    A_st = zeros(Nflt, 1);  
      A_st2 = zeros(Nflt, 1);  
 for k=Nflt:len - Nflt+1       
  
 
     %% low  
      Y_Frame_Block = LowpassR(k-Nflt+1 :k);      
      yB = ( Y_Frame_Block' * Y_Frame_Block);
      %yB = abs(Y_Frame_Block)^2;
      Y_LMS = A_st'*Y_Frame_Block;   
      err = Lowpass(k) - Y_LMS ;   
      
      Pest = Pest * 0.2 + yB*0.8;
      mu = (alpha /Pest);   
      A_st = A_st + mu*err*Y_Frame_Block;   
       
      En(k) = err;  
     
     %% high
      Y_Frame_Block = HighpassR(k-Nflt+1 :k);      
      yB = ( Y_Frame_Block' * Y_Frame_Block);           
      Y_LMS = A_st2'*Y_Frame_Block;   
      err = Highpass(k) - Y_LMS ;     
      mu = (alpha /yB);   
      A_st2 = A_st2 + mu*err*Y_Frame_Block;  
      EnH(k) = err;       
       
 end 
 
  Outbuf = halfbandInterpolator(En,EnH);     
   audiowrite([fout '.wav'],Outbuf,fs1);
 
  