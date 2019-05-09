 

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 
global SPP;

 
NFIR = 512;
Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.08;
Speed = 340;
 
finl = './voice/t15';
%finr = './voice/t14r';


%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_90_gsc2'];
 
 
[x,fs1]= audioread([finl '.wav']);  % main mic
%[x2,fs2]= audioread([finr '.wav']); % ref mic

x1 = x(:,1);
x2 = x(:,2);

fs = fs1;
 
 
fprintf('1st do gsc filter \n');
  
tau_dis =  cos(Angle*pi/180.0)*Dmic;
tau = tau_dis /Speed;  %delay 
 
taulen = fs* tau;

tauinx =[-3,-2,-1, 0,1,2,3];
taudelay = tauinx ./ fs * Speed;

tau_angle = acos(taudelay./Dmic)  *180.0/pi;

lenS1 = length(x1);
lenS2 = length(x2);
lenS = min(lenS1,lenS2);
%delay one buffer
if(taulen>0)  % delay x1  
   delayPoint = fix(taulen);     
   xTwo = x2(1:lenS);  
   xOne =  zeros(lenS,1);
   xOne(delayPoint+1:lenS) = x1(1:lenS-delayPoint);
else          % delay x2
   delayPoint = fix(-taulen);
   xTwo = x1(1:lenS);   xTwo = xTwo.* (1+  delayPoint * Speed/fs);
   xOne =  zeros(lenS,1);
   xOne(delayPoint+1:lenS) = x2(1:lenS-delayPoint);   
end
  
 
 Y_Out = zeros(lenS,1);
 Y_U = zeros(lenS,1);
 Y_D = zeros(lenS,1); 
 A_st = zeros(Nflt, 1); 
  
 % blocking of GSC beamformer
 Y_Block = xOne - xTwo;
 % fbf of GSC beamformer 
 Y_Up  = 0.5*(xOne+xTwo);       
 om_cnt = 1;
 
 % omlsa parameter
 l_OMLSA = 512; q_OMLSA = l_OMLSA/4;
 Y_om_in = zeros(NFIR,1);
 U_om_in = zeros(NFIR,1);
 yout_omlsa =  zeros(lenS,1);
 k_omlsa = 0;
 
 % gsc + dual omlsa  process
 % block flms
 
 fstart =  1;
 fend   =  NFIR;
 
 err = zeros(lenS,1);
 
 
 l = 1;
 
 BLK = lenS / NFIR  -1 ;
 
 alpha_flms = 0.998;
 omega_flms = 0.85;
 
 estimated_power = zeros(2*NFIR,1);
 
 %input_signal = Y_Block; 
 %desired_signal = Y_Up;
 
 %[error,~ ]= FLMS(alpha_flms,NFIR,input_signal,desired_signal,omega_flms,estimated_power ) ;
 
 %error_pend = [ zeros(length(Y_Block)- length(error),1)  ;  error];
 
% O_B1 = [error_pend, Y_Block];
% audiowrite([fout '_O_B1.wav'],O_B1,fs); 
 
  for blk = 1:BLK
          
    
     d = Y_Up(blk*NFIR+1 :(blk+1)*NFIR);       
     
     if(blk==1)     
       x_old = Y_Block((blk-1)*NFIR+1 :(blk)*NFIR);    
       x_new = Y_Block((blk)*NFIR+1   :(blk+1)*NFIR);         
       err(fstart  :fend ) = flms_frame( [x_old; x_new],d ,NFIR, alpha_flms,omega_flms, fstart);  %[err] = flms_frame( x_new, d , Nfft , initflag )     
     else
       x_new = Y_Block((blk)*NFIR+1 :(blk+1)*NFIR);    
       err(fstart  :fend ) = flms_frame( x_new,d ,NFIR, alpha_flms,omega_flms, fstart);  %[err] = flms_frame( x_new, d , Nfft , initflag )              
     end
      
     
         err_temp = err(fstart  :fend ) ;
        
        for tmp = 1:4                
               
            Y_om_in  = [ Y_om_in(q_OMLSA+1:l_OMLSA); err_temp((tmp-1)*q_OMLSA+1:tmp*q_OMLSA)];     
            U_om_in  = [ U_om_in(q_OMLSA+1:l_OMLSA); x_new((tmp-1)*q_OMLSA+1:tmp*q_OMLSA)];
 
            yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );      
            yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;   k_omlsa = k_omlsa+1;            
        end
 
      
 
     
     fstart = fstart + NFIR;
     fend   = fend   + NFIR;
     
 end
  
 Y_B = [err, Y_Block ];
 audiowrite([fout '_O_B.wav'],Y_B,fs); 
  
 U_B = [Y_Up , Y_Block];
 audiowrite([fout '_U_B.wav'],U_B,fs);
 %Y_O_B = [Y_Out, Y_Block ];
 Y_U_D = [Y_Up, Y_D ];
 
 Y_U_OM = [Y_Up , yout_omlsa, ];
 audiowrite([fout '_omlsa.wav'],Y_U_OM,fs); 
 %audiowrite([fout '_OUT.wav'],Y_Out,fs); 
 %
% audiowrite([fout '_O_B.wav'],Y_O_B,fs); 
% audiowrite([fout '_U_D.wav'],Y_U_D,fs);  
%fprintf('2nd do posterfilter \n');

%audiowrite([fout '_omlsa.wav'],yout_omlsa,fs); 
 
 
fprintf('End of DUAL GSC\n');




