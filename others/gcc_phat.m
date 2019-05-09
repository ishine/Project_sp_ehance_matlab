function [enhanced_ouput]=gcc_phat(x3,x4,fs,output_fn)
 clear all;
disp('Run Gcc-phat----');
frameLength=floor(15*16000/1000);  
if (frameLength/2 ~= floor(frameLength/2))
      frameLength=frameLength+1;
end
   
[x3,fs3]= audioread('F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM1.wav');
[x4,fs4]= audioread('F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM2.wav');
 
  FFT_LEN=8192;
 frameLength = 4096;
 frameShift = 4096;
lenS=min(length(x3),length(x4));
nFrame=0;
iniFrameSample=1;
endFrameSample=iniFrameSample+frameLength-1;


%Defining bands  
  
  fid= fopen('F:/Work/2018/Beamforming/matlab/GSCLMS/voice/gcc_phac2.bin','a+');
  
  SegLen = floor(lenS / frameShift);
  
  delay_array = zeros(1, SegLen);
  angle_array = zeros(1, SegLen);
  arrInx = 1;
  
  Frame1 = zeros(1,FFT_LEN );
   Frame2 = zeros(1,FFT_LEN ); 
while endFrameSample<lenS
   
    %A new frame to process
    nFrame=nFrame+1;

    %Get short-time magnitude and phase spectrums for each input channel
    Frame1(1:frameLength)=x3(iniFrameSample:endFrameSample);
    Frame2(1:frameLength)=x4(iniFrameSample:endFrameSample);
    
    % local gccphat
    P = (fft(Frame1,FFT_LEN).*conj(fft(Frame2,FFT_LEN)));
    A = 1./abs(P);
    B = A.*P;  
    C = ifft(B); 
    R_est1 = fftshift(C); 
    R_est1 = R_est1(2:FFT_LEN);  
    
    % internal gccphat
    Frame1_1 = Frame1';
    Frame2_1 = Frame2';
    
   [tau, R, lag] = gccphat(Frame1_1,Frame2_1,16000) ;     
   
   
   % result
   DelayDist = tau*340;   
   angle =  acos(DelayDist/0.1);
   angle= angle*180/pi;
     
   timestart = iniFrameSample/16000.0;
   timestend = endFrameSample/16000.0
   
   maxval = max(abs(R));
   [~,maxinx] =max(R);
   
   delay_array(arrInx) = maxinx;
   angle_array(arrInx) = angle;
   
   arrInx=arrInx+1;
   
       %  fprintf('\n maxval\n'); 
    %  fprintf('%f ',maxval);
   fprintf(fid,' timestart: %10f Delay is: %10f DelayDist:%10f angle:%10f  maxval: %10f\n ', timestart ,tau,DelayDist, angle, maxval);
%   printf("time %10f Delay is %10f Angle is %10f  tdoa val is %10f  delay point %10d\n", (float )ANC_cnt*FRAME_LENGTH/32000.0 , MaxDelay, MaxAngle , tdoa_val[0], tdoa_delay[0]);

   
      %  fprintf('\n maxinx\n'); 
   %  fprintf('%f ',maxinx)
     
   %     fprintf('\n tau\n'); 
   %   fprintf('%f ',tau)
     
   %  fprintf('\nR_est1\n'); 
   %   fprintf('%f ',real(R_est1));
      
   
     % fprintf('\nR\n'); 
     % fprintf('%f ',real(R));
    %  fprintf(fid,' R_est1: %-6f\n ', R_est1);
     

   
    %Update frame boundaries
    iniFrameSample=iniFrameSample+frameShift;
    endFrameSample=endFrameSample+frameShift;
end

arrInx =  (arrInx-2)*frameShift/16000;

t = 0:(frameShift/16000):arrInx;
%figure(1);
%plot(t, delay_array);hold on
figure(1);
plot(t, angle_array); hold on

figure(2);

t=  -4:4;
y= acos(t* 340 / 0.1/16000) * 180 / 3.1415926;
stem(t,y);
text(t,y,num2str([y].','(%.2f)'))


fclose(fid);
 disp('Processing ... 100% Done');
 
 