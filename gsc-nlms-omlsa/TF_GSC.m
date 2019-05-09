 

close all; clc; clear all;
%*************************************************************************************************************************** 
%*  tf-gsc: signal enhancement using beamforming and nonstationarity with applications to speech; sharon gannot
%*  combine gsc and omlsa frame by frame
%***************************************************************************************************************************
fprintf('tf gsc + posterfilter \n');

global SPP;

alpha_gsc = 0.1;

fin = './voice/t10';
 fout = [fin , '-tf-gsc'];
[x,fs]= audioread([fin '.wav']);  % main mic
 
x1 = x(:,1);
x2 = x(:,2);
lenS = length(x1);

pi = 3.1415926535;
angle = 90; d = 0.08; c = 343;

tau = d*cos(angle*pi/180) / c;
tauPointFloat = tau*16000;
tauPoint = fix(tauPointFloat);

% parameter 
NFFT = 512; frameshift = NFFT/2;
win = hanning(NFFT);

n = (1:NFFT)';
%phaseshift  =  exp(-1i*((n-1)*2*pi* tauPointFix/NFFT) );  
 phase_shiftRe = cos(((n-1)*2*pi *(tauPoint) /NFFT  ) ); 
 phase_shiftIm = sin(((n-1)*2*pi *(tauPoint ) /NFFT  )); 
 phase_shift =  phase_shiftRe - phase_shiftIm * 1i;


l = 1 ;  frameHead = 1; frameEnd = NFFT;
  
hisout = zeros(frameshift,1);
out_T = zeros(lenS,1);  

out_z1 = zeros(lenS,1);   out_z2 = zeros(lenS,1);
hisoutz1 = zeros(frameshift,1);  hisoutz2= zeros(frameshift,1);

Pest = zeros(NFFT, 1);  omega = 0.985; mu = 0.001;
Gm = zeros(NFFT,1);

W_fbf = [0.5,0.5];
W_bk  = [1,-1];

N_OMLSA = 512;
q_OMLSA = N_OMLSA/4;
q2_OMLSA = 2* q_OMLSA;
q3_OMLSA = 3* q_OMLSA;
Y_om_in = zeros(N_OMLSA,1);
U_om_in = zeros(N_OMLSA,1);
yout_omlsa = zeros(lenS,1);  yout_omlsa2 = zeros(lenS,1);
% do tf-gsc
his_zeta = zeros(frameshift,1);  zeta_t = zeros(lenS,1);
k = 1;
while(frameEnd < lenS)
    
    % input data read 
    x_1 = x1(frameHead:frameEnd);
    x_2 = x2(frameHead:frameEnd);
    Z1 = fft(win.*x_1);   Z1 = Z1 .* phase_shift ; % delay
    Z2 = fft(win.*x_2);
   
    % Wfbf   Wbk    
    Y_FBF   = (Z1 + Z2)*0.5;    
    U       =  Z1 - Z2;
    
    ZETA = Y_FBF - U;
    zeta  = ifft(ZETA);
    zeta_t((l-1)*frameshift+1 :(l)*frameshift)  = zeta(1:frameshift) + his_zeta ;
    his_zeta =  zeta( frameshift+ 1 : NFFT);
    % do freq nlms
    
 %   Y = Y_FBF - Gm.* U;          
 %   Pest = omega .* Pest  + (1-omega).*( abs(Z1).^2  +  abs(Z2).^2 );  
 %   gain = mu .* U .* conj(Y) ./ (Pest+1e-10);    
 %   Gm = Gm + gain;
    
  %  X_lms =  Gm.* U;          
    X_lms =  Gm.* ZETA;        
    Y = Y_FBF - X_lms; 
    Pest = omega * Pest + (1-omega) * ( abs(Z1).^2  +  abs(Z2).^2 );   
    mu = alpha_gsc ./ (Pest + 1e-10);    
    Gm_hat = Gm + mu .* Y .* conj( U);
    
    G_hat_f = fft(Gm_hat);        
    G_hat_f  = [  zeros(130,1); G_hat_f(131:382)  ;  zeros(130:1)]; % 512 ->       
    Gm = ifft(G_hat_f,512);
     
    yout = real(ifft(Y));
     
    out_T((l-1)*frameshift+1 :(l)*frameshift) = yout(1:frameshift) + hisout;    
    hisout = yout(frameshift+1: NFFT);
     
    
    %% just for test
    z1out = real(ifft(X_lms));    
    out_z1((l-1)*frameshift+1 :(l)*frameshift) = z1out(1:frameshift) + hisoutz1;    
    hisoutz1 = z1out(frameshift+1: NFFT);
    
    z2out = real(ifft(U));  
    out_z2((l-1)*frameshift+1 :(l)*frameshift) = z2out(1:frameshift) + hisoutz2;    
    hisoutz2 = z2out(frameshift+1: NFFT);
    

%%   apply omlsa 
     out_T_temp =  out_T((l-1)*frameshift+1 :(l)*frameshift);    
     out_z2_temp = out_z2((l-1)*frameshift+1 :(l)*frameshift);     
       
    Y_om_in  =  [ Y_om_in(q_OMLSA+1 : N_OMLSA) ; out_T_temp(1:q_OMLSA) ];      
    U_om_in  =  [ U_om_in(q_OMLSA+1 : N_OMLSA) ; out_z2_temp(1:q_OMLSA) ];  
    
    %  1st omlsa 128 input
     yout_omlsa_t =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,l );
     yout_omlsa((k-1)* q_OMLSA+1: k*q_OMLSA) = yout_omlsa_t;
      
     Y_om_in  =  [ Y_om_in(q_OMLSA+1 : N_OMLSA) ; out_T_temp(  q_OMLSA+1: q2_OMLSA) ];      
     U_om_in  =  [ U_om_in(q_OMLSA+1 : N_OMLSA) ; out_z2_temp( q_OMLSA+1 :q2_OMLSA) ];  
     
    % 2nd omlsa 128 sample 
     yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,l );
     yout_omlsa( k*q_OMLSA+1 :(k+1)*q_OMLSA ) = yout_omlsa_t;
     
     k = k +2;
     l = l+1;
     frameHead = frameHead + frameshift; frameEnd = frameEnd + frameshift;

end


 out2 = [ out_T , out_z2];
 
audiowrite([fout '-yout.wav'],out2,fs);

outz = [out_z1, out_z2];

audiowrite([fout '-z.wav'],outz,fs); 
 
  audiowrite([fout '-omlsa.wav'],yout_omlsa,fs);

   outzeta = [ zeta_t , out_z2];
  
audiowrite([fout '-zeta.wav'],outzeta,fs);

%audiowrite([fout '-omlsa2.wav'],yout_omlsa2,fs); 
fprintf('End of freq delay\n');



