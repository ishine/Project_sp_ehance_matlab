 

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 
global SPP;

   SPP = zeros(1,1);

Nmic = 2; % mic number 
Nflt = 128; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.095;
Speed = 340;
 
finl = './voice/t19';
%finr = './voice/t14r';


%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_90_gsc3'];
 
 
[x,fs1]= audioread([finl '.wav']);  % main mic
%[x2,fs2]= audioread([finr '.wav']); % ref mic

x1 = x(:,1);
x2 = x(:,2);

fs = fs1;


fprintf('1st do gsc filter \n');
%x3 = x1 + x2;

%x4 = x1-x2;

%define the parameter


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
 

x12 = [ xTwo,xOne];
 audiowrite([fout '_x12.wav'],x12,fs); 
 
%prepare the init parameter

    f_st = zeros(Nflt,1);
    
    f_st(1:4) =  [1,0,0,0];
    
   % f_st = [1,-2,1.5,2, 0,0,0, 0,   0,0,0, 0,   0,0,0, 0  ]';   
   % f_st = [1,-2,1.5,2, 0,0,0, 0,   0,0,0, 0,   0,0,0, 0  ]';   
   % f_st = zeros(Nflt,1); f_st(1) = 1;
  
f_st = f_st/norm(f_st);
 
C = zeros(Num_c,Nflt);
c = ones(2,1);
 
i=1;j=1;
while i<=Nflt   
    C(j:j+1,i) = c;
    i=i+1;
    j= j+2;
end
 Y_Out = zeros(lenS,1);
 Y_U = zeros(lenS,1);
 Y_D = zeros(lenS,1); 
 A_st = zeros(Nflt, 1); 
 %define the blocking matrix
 Ws = zeros(Nmic-1, Nmic);
 
 c(2) = -1;
 c= c';
 i=1;j=1;
 while i<=Nmic-1   
    Ws(i,j:j+1) = c;
    i=i+1;
    j= j+1;
 end
  
 Pest = 0;
 % down line of GSC beamformer
 Y_Block = xOne - xTwo;

 % up line of GSC beamformer 
 Y_Up  = (xOne+xTwo)  /2; % Y_Up 
 
 
% lms = dsp.LMSFilter(   256);
% [~,e] = lms( Y_Block,Y_Up); 
% outfx = [e,Y_Block]; 
% audiowrite([fout 'xlms.wav'],outfx,fs); 
%  [y,out]=gsc_dual_postfilter([fout 'xlms'] );
% outfx = [e, out];   
% audiowrite([fout 'xlms_om.wav'],outfx,fs); 

K = 512;
alpha = 0.073;  
om_cnt = 1;
 
 % omlsa parameter
 NFFT = 512; q_OMLSA = NFFT/4;
 Y_om_in = zeros(NFFT,1);
 U_om_in = zeros(NFFT,1);
 yout_omlsa =  zeros(lenS,1);
 k_omlsa = 0;
 
 yout_omlsa_t_vec =  zeros(Nflt,1);
 Pest_om = 0;
 Pest_E  = 0;
 
 errork = zeros(lenS,1);
 yout_omlsa_t_vec_array =  zeros(lenS,1);
 power1 =  zeros(lenS,1);
 power2 =  zeros(lenS,1);
 power_error =  zeros(lenS,1);
 power_om =  zeros(lenS,1);
 power_B=  zeros(lenS,1);
 
 P_ratio =  1e-8 ;
 
 muarray =  zeros(lenS,1);
 muarray_om =  zeros(lenS,1);
 
 P_ratio_array=  zeros(lenS,1);
 divergeState = 0;
 ekEn = 0; dkEn =0;
 % gsc + dual omlsa  process
 for k=Nflt:lenS - Nflt+1    
    % down line     
     Y_Frame_Block = Y_Block(k-Nflt+1 :k);       
     Y_Frame_up = Y_Up(k-Nflt+1 :k);         
     yB =  Y_Frame_up' * Y_Frame_up;     
      
     Y_Down = A_st'*Y_Frame_Block;        
     yout = Y_Up(k) - Y_Down ;        
     errork(k) = yout;
      
     omega = 0.898;
     ekEn =  omega* ekEn + (1-omega) * yout*yout;
     dkEn =  omega* dkEn + (1-omega) * Y_Up(k)*Y_Up(k);
 
       
     %% improved nlms
     omega = 0.898; 
     mode = 1;
     if(mode == 1 ) % orig nlms power estimate
    
     Pest = omega * Pest + (1-omega) * yB;          
     mu = alpha ./( Pest + 1e-2);
    
     else % omlsa assist power estimate
      alpha = 0.158;
      Pest = omega * Pest + (1-omega) * yB;          
      mu = alpha ./( Pest + 1e-2);  
      muarray(k) = mu;
      P_ratio_ =  Pest_E ./ (Pest_om + 1e-10);
      
     P_ratio  =   omega * P_ratio + (1- omega) * P_ratio_ ;
     P_ratio_num =    P_ratio  /( P_ratio + 6  );
     
     if(P_ratio_num < 0.5)
         P_ratio_num = 0;
     end
     
     mu = mu .* P_ratio_num;
     muarray_om(k) = mu;
    
     P_ratio_array(k) = P_ratio_num;
     
     Power_E =  ( errork(k-Nflt+1 :k)' *  errork(k-Nflt+1 :k));
     Pest_E =  omega *  Pest_E + (1-omega) * Power_E;     
      
     power_error(k) =  Pest_E;     
     power_B(k) =    (Pest );          
     end 
      
     A_st = A_st +  mu* yout.* (Y_Frame_Block) ;  
 
     Y_om_in  = [ Y_om_in(2:NFFT) ; yout];
     U_om_in  = [ U_om_in(2:NFFT);  Y_Block(k)];
  
     if(om_cnt==q_OMLSA)     
       yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );  om_cnt = 0;    
       %yout_omlsa_t  =  tf_gsc_dual_postfilter2(Y_om_in,U_om_in,k_omlsa+1 );  om_cnt = 0;    
        
         %yout_omlsa_t  =  tf_gsc_dual_postfilter(U_om_in,Y_om_in,k_omlsa+1 );  om_cnt = 0; 
      %  yout_omlsa_t_vec = [ yout_omlsa_t_vec(q_OMLSA+1:2*q_OMLSA);yout_omlsa_t ];
        
        for tmp = 1: q_OMLSA                
              yout_omlsa_t_vec = [ yout_omlsa_t_vec(2:Nflt);yout_omlsa_t(tmp) ];
              Power_om  = yout_omlsa_t_vec' * yout_omlsa_t_vec;
              Pest_om = omega * Pest_om + (1-omega)* Power_om;
              
              power_om(k_omlsa * q_OMLSA + tmp) =   (Pest_om );    
              
        end
         
        yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;       
        k_omlsa = k_omlsa+1;        
        
     end
     om_cnt = om_cnt+ 1;
     
     Y_Out(k) = yout;    
     Y_U(k) = Y_Up(k) ;
     Y_D(k) = Y_Down;
     
 end
  
 Y_U_B = [Y_Up , Y_Block];
 Y_O_B = [Y_Out, Y_Block ];
 Y_U_O = [Y_Up, Y_Out ];
 
 audiowrite([fout '_U_B.wav'],Y_U_B,fs);
 audiowrite([fout '_O_B.wav'],Y_O_B,fs); 
 audiowrite([fout '_U_O.wav'],Y_U_O,fs);  
  
 Y_Om_Block = [yout_omlsa, Y_Out];
  
 audiowrite([fout '_Om_D.wav'],Y_Om_Block,fs); 
 
 Y_Om_out = [ yout_omlsa,Y_Block   ];
  
audiowrite([fout '_Om_Block.wav'],Y_Om_out,fs); 
   
power_fig = 1;

if (power_fig==1)

P1 =  (   Y_Out(1:length(Y_Out)-130+1).^2);
P2 =   (  yout_omlsa(130:length(yout_omlsa)).^2);
 Y_Out_om_d =   P1./ (P2 + 1e-10);
 
 Y_Out_om_d = Y_Out_om_d ./ (Y_Out_om_d + 1e-2);
  
 Y_Om_error = [P2 , P1 ];
  
  audiowrite([fout '_Om_error.wav'],Y_Om_error,fs); 
  Y_Om_error = [P2 , Y_Out_om_d ];  
  audiowrite([fout '_Om_error2.wav'],Y_Om_error,fs); 
  
   figure(1);
   x = 1: length(Y_Out_om_d);
   x1 = 1: length(Y_Up);
   subplot(411  );   plot(  x, Y_Out_om_d ,'b' );   title('1) error omlsa power ratio ')
   subplot(412   );  plot(  x1, Y_Up.^2  ,'r'  );    title('2) up power')
   subplot(413   );  plot(  x1, Y_Out.^2  ,'r'  );   title('2) error power')
   subplot(414   );  plot(  x1, Y_Block.^2  ,'r'  ); title('2) block power')
   
%   subplot(513   );  plot(  x, muarray(1:length(Y_Out_om_d)) ,'g'  );  title('3) omlsa power')
%   subplot(514   );  plot(  x, muarray_om(1:length(Y_Out_om_d)) ,'g'  );  title('3) omlsa power')
%   subplot(515   );  plot(  x, P_ratio_array(1:length(Y_Out_om_d)) ,'g'  );  title('3) omlsa power')
 end
fprintf('End of DUAL GSC\n');




