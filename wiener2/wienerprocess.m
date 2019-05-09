function [enhanced_ouput]=wienerprocess(x1,fs1,output_fn)

disp('Run Speech Enhacement Wiener----');

[x1,fs1]= audioread('F:/Work/2018/Beamforming/matlab/WienerScalart96/Y_TF_GSC_OUT.wav');

frameLength = 2048;
frameShift = 2048;
lenS=length(x1);
IS =0.25*16000;
enhanced_ouput=zeros(lenS,1);
nFrame=0;
iniFrameSample=1;
endFrameSample=iniFrameSample+frameLength-1;

%while endFrameSample<lenS
         
  %  Frame1=x1(iniFrameSample:endFrameSample);
   
   % lensin = length(Frame1);
    
 %   enhanced_ouput=WienerScalart96(x1,fs1,IS);
    [enhanced_ouput1, enhanced_ouput2] =  WienerNoiseReduction(x1,fs1, IS);
   % lengthoutput(iniFrameSample:endFrameSample) = enhanced_ouput_tmp;
        
  %  iniFrameSample=iniFrameSample+frameShift;
   % endFrameSample=endFrameSample+frameShift;    
%end
  %  lengthoutput = length(enhanced_ouput);
 audiowrite('f:/Work/2018/Beamforming/matlab/wiener2/043Wiener2.wav',enhanced_ouput1,16000);
 audiowrite('f:/Work/2018/Beamforming/matlab/wiener2/043Wiener2_2.wav',enhanced_ouput2,16000);
 
 