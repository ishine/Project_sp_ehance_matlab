function [enhanced_ouput]=COH_RealImag( )

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
 
    
if exist('labelpf','var')
    fprintf('Frame1:\n ' );
    fprintf('%f,', Frame1);     
    fprintf(':\n ')
 end
    
    X1=fft(wFrame1,FFT_LEN);
    X2=fft(wFrame2,FFT_LEN);
 
 if exist('labelpf','var')   
     fprintf('X1:\n ' );
    fprintf('%f,', real(X1));     
    fprintf(':\n ')  
    
    fprintf('%f,', imag(X1));     
    fprintf(':\n ')    
    
     fprintf('X2:\n ' );
    fprintf('%f,', real(X2));     
    fprintf(':\n ')  
    
    fprintf('%f,', imag(X2));     
    fprintf(':\n ')      
   end    
    
    if (nFrame==1)
        PX1X1=abs(X1).^2;
        PX2X2=abs(X2).^2;
        PX1X2=X1.*conj(X2); 
    else
        PX1X1=lambda_x.*PX1X1+(1-lambda_x).*abs(X1).^2;
        PX2X2=lambda_x.*PX2X2+(1-lambda_x).*abs(X2).^2;
        PX1X2=lambda_x.*PX1X2+(1-lambda_x).*X1.*conj(X2);
    end
   
 if exist('labelpf','var')  
    fprintf('PX1X1:\n ' );
    fprintf('%f,', PX1X1);     
    fprintf(':\n ') 
    
    fprintf('PX2X2:\n ' );
    fprintf('%f,', PX2X2);     
    fprintf(':\n ')      
   end
   
     cohX=PX1X2 ./ sqrt(PX1X1.*PX2X2+epsilon); 
     reCOH=real(cohX(2:FFT_LEN/2+1));
     imCOH=imag(cohX(2:FFT_LEN/2+1));
     
  
      
    G1=1-abs(reCOH).^P;   %for suppressing noise coming from angles about 90
 
 
 
    G2=ones(FFT_LEN/2,1);               %for suppressing noise coming from angles greater than 90
    ind_neg= imCOH(1:FFT_LEN/2) < limNeg; 
    G2(ind_neg)=negImag_ConsFilter;   
    
  
    G=G1.*G2;                            %Halfband final filter   
    H=abs([G ;flipud(G)]);    %Fullband final filter 
      
   if exist('labelpf','var')
     fprintf('G:\n ' );
    fprintf('%f,',G);     
    fprintf(':\n ');   
       end
    
    HX= H .* X1;
    
    
if exist('labelpf','var')  
    fprintf('HX:\n ' );
    fprintf('%f,', real(HX));     
    fprintf(':\n ');   
 
    fprintf('%fi,', imag(HX));     
    fprintf(':\n ');
   end
    %IFFT and OLA
    enhSpeech_Frame_tmp=real(ifft(HX,FFT_LEN));  
     
    
if exist('labelpf','var')
    fprintf('enhSpeech_Frame_tmp:\n ' );
    fprintf('%f,', enhSpeech_Frame_tmp);     
    fprintf(':\n ');  
      end
    enhSpeech_Frame=enhSpeech_Frame_tmp(1:frameLength);
    
if exist('labelpf','var')
    fprintf('enhSpeech_Frame:\n ' );
    fprintf('%f,', enhSpeech_Frame );     
    fprintf(':\n '); 
       end
    enhanced_ouput(iniFrameSample:endFrameSample)=enhSpeech_Frame + enhanced_ouput(iniFrameSample:endFrameSample);   
    
    
if exist('labelpf','var')  
    fprintf('enhanced_ouput:\n ' );
    fprintf('%f,', enhanced_ouput(iniFrameSample:endFrameSample));     
    
    fprintf(fid,'%f,\n ', enhanced_ouput(iniFrameSample:endFrameSample));
    fprintf(':\n '); 

    
    fclose(fid);
   end    
    
    enhanced_ouputFix =  32767.*enhanced_ouput;
    
    
if exist('labelpf','var')  
    fprintf('enhanced_ouput:\n ' );
    fprintf('%f,', enhanced_ouputFix(iniFrameSample:endFrameSample));     
    
    fprintf(':\n '); 
   end
    
     
    %Update frame boundaries
    iniFrameSample=iniFrameSample+frameShift;
    endFrameSample=endFrameSample+frameShift;
end

 disp('Processing ... 100% Done');
 wavwrite(enhanced_ouput,16000,16,'f:/Work/2018/MIC_ARRAY_ANC/matlabfile/COH_6CM_3.wav');
