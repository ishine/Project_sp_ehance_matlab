function [y,out]=omlsa(fin,fout)
% omlsa : Single Channel OM-LSA with IMCRA noise estimator
% ***************************************************************@
% Inputs:
%    fin,  input file name (fin.wav)
%    fout, output file name (fout.wav)
% Output:
%    in,  samples of the input file
%    out, samples of the output file
% Usage:
%    [in, out]=omlsa(fin,fout);
%    omlsa(fin);
% Defaults:
%    fout= [fin,'_omlsa'];
%
% Copyright (c) 2003. Prof Israel Cohen.
% All rights reserved. 
% ***************************************************************@
fin = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/1'; %040_ANC_PC2.wav

nin=nargin;
if nin<2
    fout=[fin,'_stft'];
end

% 1) Parameters of Short Time Fourier Analysis:
Fs_ref=16e3;		% 1.1) Reference Sampling frequency
M_ref=512;		% 1.2) Size of analysis window
Mo_ref=0.5*M_ref;	% 1.3) Number of overlapping samples in consecutive frames

% 2) Parameters of Noise Spectrum Estimate
w=1;			% 2.1)  Size of frequency smoothing window function=2*w+1
alpha_s_ref=0.9;	% 2.2)  Recursive averaging parameter for the smoothing operation
Nwin=8; 	% 2.3)  Resolution of local minima search
Vwin=15;
delta_s=1.67;		% 2.4)  Local minimum factor
Bmin=1.66;
delta_y=4.6;		% 2.4)  Local minimum factor
delta_yt=3;
alpha_d_ref=0.85;	% 2.7)  Recursive averaging parameter for the noise

% 3) Parameters of a Priori Probability for Signal-Absence Estimate
alpha_xi_ref=0.7;	% 3.1) Recursive averaging parameter
w_xi_local=1; 	% 3.2) Size of frequency local smoothing window function
w_xi_global=15; 	% 3.3) Size of frequency local smoothing window function
f_u=10e3; 		% 3.4) Upper frequency threshold for global decision
f_l=50; 		% 3.5) Lower frequency threshold for global decision
P_min=0.005; 		% 3.6) Lower bound constraint
xi_lu_dB=-5; 	% 3.7) Upper threshold for local decision
xi_ll_dB=-10; 	% 3.8) Lower threshold for local decision
xi_gu_dB=-5; 	% 3.9) Upper threshold for global decision
xi_gl_dB=-10; 	% 3.10) Lower threshold for global decision
xi_fu_dB=-5; 	% 3.11) Upper threshold for local decision
xi_fl_dB=-10; 	% 3.12) Lower threshold for local decision
xi_mu_dB=10; 	% 3.13) Upper threshold for xi_m
xi_ml_dB=0; 		% 3.14) Lower threshold for xi_m
q_max=0.998; 		% 3.15) Upper limit constraint

% 4) Parameters of "Decision-Directed" a Priori SNR Estimate
alpha_eta_ref=0.95;	% 4.1) Recursive averaging parameter
eta_min_dB=-18;	% 4.2) Lower limit constraint

% 5) Flags
broad_flag=1;               % broad band flag   % new version
tone_flag=1;                % pure tone flag   % new version
nonstat='medium';                %Non stationarity  % new version

% Read input data
[y,Fs]=audioread([fin,'.wav'])  % read size of input data, Fs and NBITS
N = length(y);

% Adjust parameters according to the actual sampling frequency
if Fs~=Fs_ref
    M=2^round(log2(Fs/Fs_ref*M_ref));
    Mo=Mo_ref/M_ref*M;
   
else
    M=M_ref;
    Mo=Mo_ref;
 
end
 
% window function
win=hamming(M);
%win=hanning(M);
% find a normalization factor for the window
win2=win.^2;
Mno=M-Mo;
W0=win2(1:Mno);
for k=Mno:Mno:M-1
    swin2=lnshift(win2,k);
    W0=W0+swin2(1:Mno);
end
W0=mean(W0)^0.5;
%win=win/W0;
Cwin=sum(win.^2)^0.5;
%win=win/Cwin;

Nframes=fix((N-Mo)/Mno);   %  number of frames
out=zeros(M,1);
 
M21=M/2+1;
 
 
l_fnz=1;      % counter for the first frame which is non-zero    % new version omlsa3
fnz_flag=0;     % flag for the first frame which is non-zero    % new version omlsa3
zero_thres=1e-10;      % new version omlsa3
% zero_thres is a threshold for discriminating between zero and nonzero sample.
% You may choose zero_thres=0, but then the program  handles samples which are identically zero (and not ¡°epsilon¡± values).

finID = fopen([fin '.wav'],'r');
for l=1:Nframes
    if l==1
        %[y,Fs,~,fidx]=readwav(fin,'rf',M,0);    % open input file and read one frame of data
        initfin=fread(finID,22,'int16');
        y=fread(finID,M,'int16')/2^15;
    else
        %[y0,Fs,~,fidx]=readwav(fidx,'rf',Mno);
        y0=fread(finID,Mno,'int16')/2^15;
        y=[y(Mno+1:M,:) ; y0];      % update the frame of data
    end
    
    if (~fnz_flag && abs(y(1))>zero_thres) ||  (fnz_flag && any(abs(y)>zero_thres))       % new version omlsa3
        fnz_flag=1;     % new version omlsa3

        % 1. Short Time Fourier Analysis
        Y=fft(win.*y);     
        X = [zeros(3,1);Y(4:M21-1); 0 ];      
        X(M21+1:M)=conj(X(M21-1:-1:2)); %extend the anti-symmetric range of the spectum
%        x=Cwin^2  *win.*real(ifft(X));
        x=   real(ifft(Y));
        out=out+x;
    else        % new version omlsa3
        if ~fnz_flag        % new version omlsa3
            l_fnz=l_fnz+1;        % new version omlsa3
        end         % new version omlsa3
    end         % new version omlsa3
    if l==1
        foutID = fopen([fout '.wav'],'w');      % open output file and write first output samples
        fwrite(foutID,[initfin; out(1:Mno)*2^15],'int16');
    else
        fwrite(foutID,out(1:Mno)*2^15,'int16');   % write to output file
    end
    out=[out(Mno+1:M); zeros(Mno,1)];   % update output frame
    
end
 
fwrite(foutID,out(1:M-Mno)*2^15,'int16');  % write to output file the last output samples
fclose(finID);    % close input file
fclose(foutID);   % close output file
 

function y = lnshift(x,t)
% lnshift -- t circular left shift of 1-d signal
%  Usage
%    y = lnshift(x,t)
%  Inputs
%    x   1-d signal
%  Outputs
%    y   1-d signal
%        y(i) = x(i+t) for i+t < n
%	 		y(i) = x(i+t-n) else
%
% Copyright (c) 2000. Prof Israel Cohen.
% All rights reserved. 
% ***************************************************************@
szX=size(x);
if szX(1)>1
    n=szX(1);
    y=[x((1+t):n); x(1:t)];
else
    n=szX(2);
    y=[x((1+t):n) x(1:t)];
end