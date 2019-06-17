%function [ sigpower,noisePowerArray] = noisePowProposed(noisy,fs)
function [noise_power] =   noisePowProposedSingal( Ya2, initflag )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% propose SPP algorithm to estimate the spectral noise power
%%%% papers: "Unbiased MMSE-Based Noise Power Estimation with Low Complexity and Low Tracking Delay", IEEE TASL, 2012 
%%%% "Noise Power Estimation Based on the Probability of Speech Presence", Timo Gerkmann and Richard Hendriks, WASPAA 2011
%%%% Input :        noisy:  noisy signal
%%%%                   fs:  sampling frequency
%%%%                   
%%%% Output:  noisePowMat:  matrix with estimated noise power for each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Copyright (c) 2011, Timo Gerkmann
%%%%%%%%%%%%%%%%%%%%%% Author: Timo Gerkmann and Richard Hendriks
%%%%%%%%%%%%%%%%%%%%%% Universitaet Oldenburg
%%%%%%%%%%%%%%%%%%%%%% KTH Royal Institute of Technology
%%%%%%%%%%%%%%%%%%%%%% Delft university of Technology
%%%%%%%%%%%%%%%%%%%%%% Contact: timo.gerkmann@uni-oldenburg.de
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Universitaet Oldenburg, Delft university, KTH 
%	Royal Institute of Technology nor the names of its contributors may be 
% 	used to endorse or promote products derived from this software without 
% 	specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Version v1.0 (October 2011)
%finl = '../../voice/161';
     
%fout = [finl 'signal_est'];  
%fout1 = [finl 'noise_est'];  
%[noisy,fs]= audioread([finl '.wav']);  % main mic

% global declear
global q priorFact xiOptDb xiOpt logGLRFact GLRexp GLR PH_est  noisePow PH_estmean estimate alphaPH_estmean;

%% some constants

fs = 16000;
frLen   = 32e-3*fs;  % frame size
fShift  = frLen/2;
 
M21 = frLen/2+1;
anWin  = hanning(frLen,'periodic'); %analysis window

%% allocate some memory
%noisePowMat = zeros(frLen/2+1,nFrames);
%noisePowerArray = zeros(length(noisy),1);
%sigpower= zeros(length(noisy),1);
%% initialize
%noisePow = init_noise_tracker_ideal_vad(noisy,frLen,frLen,fShift, anWin); % This function computes the initial noise PSD estimate. It is assumed that the first 5 time-frames are noise-only.
%noisePowMat(:,1)=noisePow;
% 
 
alphaPSD = 0.8;
 
if initflag == 1
%constants for a posteriori SPP
alphaPH_estmean = 0.9;
q          = 0.5; % a priori probability of speech presence:
priorFact  = q./(1-q);
xiOptDb    = 15; % optimal fixed a priori SNR for SPP estimation
xiOpt      = 10.^(xiOptDb./10);
xi_inv     =  1./(1+xiOpt);


logGLRFact = log(1./(1+xiOpt));
GLRexp     = xiOpt./(1+xiOpt);

GLR     = zeros(M21,1);
PH_est     = zeros(M21,1);
PH_estmean = zeros(M21,1) + 0.5;
estimate =  zeros(M21,1);
 
 
end

if(initflag< 5)
     
    noisePow =  Ya2;%abs( noisyDftFrame) .^2;    
 
end
 
%for indFr = 1:nFrames
   % indices       = (indFr-1)*fShift+1:(indFr-1)*fShift+frLen;%  50% overlap frame
   % noisy_frame   = anWin.*noisy(indices);
  %  noisy_frame   = anWin.*noisey ;
  %  noisyDftFrame = fft(noisy_frame,frLen);
  %  noisyDftFrame = noisyDftFrame(1:frLen/2+1);  % get the noise signal in freq zone
	
   % noisyPer = noisyDftFrame.*conj(noisyDftFrame); % cal power like Y_2
   % noisyPer =abs( noisyDftFrame) .^2;  
   
    noisyPer = Ya2;
    snrPost1 =  noisyPer./(noisePow);% a posteriori SNR based on old noise power estimate
    %% noise power estimation
	% GLR     = priorFact .* exp(min(logGLRFact + GLRexp.*snrPost1,200));% gerneralized likehood ratio
	GLR = priorFact .* min(  exp(200),  xi_inv.* exp(GLRexp.*snrPost1)); 
    
    PH_est     = GLR./(1+GLR); % a posteriori speech presence probability

	PH_estmean  = alphaPH_estmean * PH_estmean + (1-alphaPH_estmean) * PH_est;
	stuckInd = PH_estmean > 0.99;
	PH_est(stuckInd) = min(PH_est(stuckInd),0.99);
	estimate =  PH_est .* noisePow + (1-PH_est) .* noisyPer ;% eq 22
	noisePow = alphaPSD *noisePow+(1-alphaPSD)*estimate;
    
    noise_power  = noisePow;
 
  
    end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% EOF