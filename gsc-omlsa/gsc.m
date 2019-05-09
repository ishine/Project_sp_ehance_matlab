clear all; close all;
 
fin0 = '.\voice\t19-0';
fin1 = '.\voice\t10r';
  
fout  = [fin0 '_lcmv'];
 
fout1  = [fin0 '_gsc'];
 
Lw = 256;
%[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
%[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic

 Y =  audioread([fin0 '.wav']);
 
 Y_L = Y(:,1);
 Y_R = Y(:,2);

c = 340.0;
fs = 16.0e3;
fc = fs/2;
lam = c/fc;
%lam = 0.02;
transducer = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 20000]);
array = phased.ULA('Element',transducer,'NumElements',2,'ElementSpacing',lam);


t = 0:1/fs:.5;
signal = chirp(t,0,0.5,500);


collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c, ...
    'SampleRate',fs,'ModulatedInput',false,'NumSubbands',512);
incidentAngle = [0;0];
signal = collector(signal.',incidentAngle);
noise = 0.5*randn(size(signal));
recsignal = signal + noise;

recsignal = [Y_L, Y_R];


 frostbeamformer = phased.FrostBeamformer('SensorArray',array,'PropagationSpeed', ...
     c,'SampleRate',fs,'Direction',incidentAngle,'FilterLength',128);
 yfrost = frostbeamformer(recsignal);
audiowrite([fout '._lcmv.wav'],yfrost,fs);
%mvdrbeamformer =phased.SubbandMVDRBeamformer('SensorArray',array,'PropagationSpeed', ...
 %   c,'SampleRate',fs,'Direction',incidentAngle   );
  
%ymvdr =  mvdrbeamformer(recsignal);
 
 gscbeamformer = phased.GSCBeamformer('SensorArray',array, ...
     'PropagationSpeed',c,'SampleRate',fs,'Direction',incidentAngle, ...
    'FilterLength',256);
 ygsc = gscbeamformer(recsignal);
%plot(t*1000,recsignal(:,6),t*1000,yfrost,t,ygsc)
%xlabel('Time (ms)')
%ylabel('Amplitude')
 audiowrite([fout '_gsc.wav'],ygsc,fs);
% audiowrite([fout1 '.wav'],ygsc,fs);
