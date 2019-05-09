%% Based on 'Speech enhancement for non-stationary noise enviroments' Isral Cohen. 2001 Elsevier.
%  1) 'om-lsa' to estimate speech present probability
%  2) 'mcra' to estimate noise spectral 


%fin0 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM1';
%fin1 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM2';
fin0 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/1';

fOut = [fin0 '_OMLSA_YANG_OUT'];
 
[x1,Fs]= audioread([fin0 '.wav']); % main mic
%[x2,fs2]= audioread([fin1 '.wav']); % ref mic

Lens = length(x1);
FrameLen = 512;
FrameShift = FrameLen/2;
FFT_LEN=2^nextpow2(FrameLen);
FrameNum = Lens / FrameShift; 
M = 512; 
M21=M/2+1;

ex_gmma = ones(M21,1);

%% Noise Spectram Extimation
alpha_s = 0.8;  S = zeros(M21,1);
alpha_p1 = 0.2;  P1_est = zeros(M21,1);
delta = 5;       I = zeros(M21,1);
alpha_d = 0.95;  alpha_d1 = zeros(M21,1); lamda_d = zeros(M21 , 1);
k2_local=round(500/Fs*M+1); % 500hz
k3_local=round(3500/Fs*M+1); %3500hz
tone_flag = 1; 
%% Gain Computation
alpha_eta = 0.95; eta = zeros(M21,1);
beta = 0.7; xi = zeros(M21, 1);
qmax = 0.95;
xi_min_dB = -10; xi_max_dB = -5;  
xi_p_min_dB = 0; xi_p_max_dB = 10; P_min = 0.005;
w = 1; w_xi_local = 1; w_xi_global = 15;
GH1 =  ones(M21, 1);
G = ones(M21, 1);Gmin = 10^-25;

P =zeros(M21,1);
X_G = zeros(FFT_LEN,1);


%% parameter

window = hanning(FrameLen);
 
FrameCnt = 0;
x_frame = zeros(FrameLen,1);
history = zeros(FrameShift,1);
x_out = zeros(Lens,1);

 init_frame = 1;
 end_frame =  init_frame+FrameLen-1;

while(end_frame<Lens)

    % win and fft
    x_frame = x1(init_frame:end_frame);    
    
    X_F = fft(window.* x_frame);     
     
    x_i = ifft(X_F);    
  %  x_out(FrameCnt*FrameShift+1:(FrameCnt+1)*FrameShift) = history + x_i(1:FrameShift);     
 %   history = x_i(FrameShift+1:FFT_LEN);
 
    x_frame =   x_i;
    x_out(FrameCnt*FrameShift+1:(FrameCnt+1)*FrameShift) = history +  x_frame(1:FrameShift);
    
    history = x_frame(FrameShift+1:FFT_LEN);
     
    init_frame = init_frame + FrameShift;
    end_frame  = end_frame + FrameShift;
    FrameCnt = FrameCnt+1;    
end

 audiowrite([fOut '.wav'],x_out,fs1);

