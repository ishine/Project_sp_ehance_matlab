function [enhanced_ouput]=COH_FunctionV1( )

% x1 , x2 : Input signals at two channels
% fs : Sampling Frequency
% output_path: path to write enhanced file

% Nima Yousefian , Sep 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference: Nima Yousefian, Philipos C. Loizou: A Dual-Microphone 
%Speech Enhancement Algorithm Based on the Coherence Function. IEEE 
%Transactions on Audio, Speech & Language Processing 20(2): 599-609 (2012)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%syms labelpf;

disp('Run Speech Enhacement COH----');
frameLength=floor(15*16000/1000);  
if (frameLength/2 ~= floor(frameLength/2))
      frameLength=frameLength+1;
end

fin0 = '..\..\voice\11';
  
[x,fs]= audioread([fin0 '.wav']); % main mic
 

x3 = x(:,1);
x4 = x(:,2);

%[x3,fs3,nbits3]= wavread('f:/Work/2018/MIC_ARRAY_ANC/matlabfile/6CM_F3.wav');
%[x4,fs4,nbits4]= wavread('f:/Work/2018/MIC_ARRAY_ANC/matlabfile/6CM_R3.wav');
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

%Defining bands  
band1=floor(1000*FFT_LEN/16000);  %0 -> 1 KHz

if exist('labelpf','var')
fprintf('frameLength is%d\n',frameLength)
fprintf('FFT_LEN is%d\n',FFT_LEN)
fprintf('lenS is%d\n',lenS)
fprintf('band is%d\n',band1)
%fid=fopen('f:/Data.txt','a+');
end

%fprintf(fid,'%f,\n ',window)
%fclose(fid);
%Defining exponents of coherence function
P=zeros(FFT_LEN/2,1);
P(1:band1)=16;
P(band1+1:FFT_LEN/2)=2;

    
%Defining threshholds for imaganary parts to consider negative (noise)
limNeg=zeros(FFT_LEN/2,1);
limNeg(1:band1)=-0.1;
limNeg(band1+1:FFT_LEN/2)=-0.3;

%Constant value of filter when imag part is negative
negImag_ConsFilter=0.05;

 
    
while endFrameSample<lenS

 if exist('labelpf','var')
    fid= fopen('F:/Work/2018/MIC_ARRAY_ANC/matlabfile/Data.txt','a+');
   end
    
    
    %A new frame to process
    nFrame=nFrame+1;

    %Get short-time magnitude and phase spectrums for each input channel
    Frame1=x3(iniFrameSample:endFrameSample);
    Frame2=x4(iniFrameSample:endFrameSample);

    wFrame1=Frame1 .* window;
    wFrame2=Frame2 .* window;
 
 
    
    X1=fft(wFrame1,FFT_LEN);
    X2=fft(wFrame2,FFT_LEN);
 
  
    
    if (nFrame==1)
        PX1X1=abs(X1).^2;
        PX2X2=abs(X2).^2;
        PX1X2=X1.*conj(X2); 
    else
        PX1X1=lambda_x.*PX1X1+(1-lambda_x).*abs(X1).^2;
        PX2X2=lambda_x.*PX2X2+(1-lambda_x).*abs(X2).^2;
        PX1X2=lambda_x.*PX1X2+(1-lambda_x).*X1.*conj(X2);
    end
 
     cohX=PX1X2 ./ sqrt(PX1X1.*PX2X2+epsilon); 
     reCOH=real(cohX(2:FFT_LEN/2+1));
     imCOH=imag(cohX(2:FFT_LEN/2+1));
     
  
      
    G1=1-abs(reCOH).^P;   %for suppressing noise coming from angles about 90
 
 
 
    G2=ones(FFT_LEN/2,1);               %for suppressing noise coming from angles greater than 90
    ind_neg= imCOH(1:FFT_LEN/2) < limNeg; 
    G2(ind_neg)=negImag_ConsFilter;   
    
  
    G=G1.*G2;                            %Halfband final filter   
    
    G = 1 - G;
    
    H=abs([G ;flipud(G)]);    %Fullband final filter 
    
    
    HX= H .* X1;
    
    
 
    %IFFT and OLA
    enhSpeech_Frame_tmp=real(ifft(HX,FFT_LEN));  
     
 
    enhSpeech_Frame=enhSpeech_Frame_tmp(1:frameLength);
 
    enhanced_ouput(iniFrameSample:endFrameSample)=enhSpeech_Frame + enhanced_ouput(iniFrameSample:endFrameSample);   
    
 
     
 
    
     
    %Update frame boundaries
    iniFrameSample=iniFrameSample+frameShift;
    endFrameSample=endFrameSample+frameShift;
end

 fout  = [fin0 'coh_v1'];
 
 
 diff_enhanced_ouput = zeros(length(enhanced_ouput),1);
 
 diff_enhanced_ouput = enhanced_ouput - x3;

 audiowrite([fout '_Out.wav'],enhanced_ouput,16000);
 audiowrite([fout '_diff.wav'],diff_enhanced_ouput,16000);
 
 
 disp('Processing ... 100% Done');
  