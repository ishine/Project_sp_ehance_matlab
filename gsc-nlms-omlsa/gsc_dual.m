 

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 
global SPP;

   SPP = zeros(1,1);

Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.08;
Speed = 340;
 
finl = './voice/t15';
%finr = './voice/t14r';


%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_90_gsc'];
 
 
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
   
 
K = 512;
alpha = 0.073;  
 
 om_cnt = 1;
 
 % omlsa parameter
 NFFT = 512; q_OMLSA = NFFT/4;
 Y_om_in = zeros(NFFT,1);
 U_om_in = zeros(NFFT,1);
 yout_omlsa =  zeros(lenS,1);
 k_omlsa = 0;
 
 % gsc + dual omlsa  process
 for k=Nflt:lenS - Nflt+1    
    % down line     
     Y_Frame_Block = Y_Block(k-Nflt+1 :k);       
     Y_Frame_up = Y_Up(k-Nflt+1 :k);         
     yB =  Y_Frame_up' * Y_Frame_up;      
     Y_Down = A_st'*Y_Frame_Block;   
     
     yout = Y_Up(k) - Y_Down ;      
  
 
       
     %% improved nlms
     omega = 0.898; 
     Pest = omega * Pest + (1-omega) * yB;      
     mu = alpha ./  (Pest+1e-5);        
      
   
     %A_st = A_st + (1 - SPP) .* mu*yout.*Y_Frame_Block;
   %if(SPP <0.25) 
        A_st = A_st +  mu*yout.*Y_Frame_Block;  
   % end
     
  % end
     Y_om_in  = [ Y_om_in(2:NFFT) ; yout];
     U_om_in  = [ U_om_in(2:NFFT);  Y_Block(k)];
 
     om_cnt = om_cnt+ 1;
     
     if(om_cnt==128)     
        yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );  om_cnt = 1;       
        yout_omlsa( k_omlsa*q_OMLSA+1 :(k_omlsa+1)*q_OMLSA ) = yout_omlsa_t;       
        k_omlsa = k_omlsa+1;        
     end
      
     Y_Out(k) = yout;    
     Y_U(k) = Y_Up(k) ;
     Y_D(k) = Y_Down;
     
 end
 
  
 Y_U_B = [Y_Up , Y_Block];
 Y_O_B = [Y_Out, Y_Block ];
 Y_U_D = [Y_Up, Y_D ];
 
 audiowrite([fout '_U_B.wav'],Y_U_B,fs);
 audiowrite([fout '_O_B.wav'],Y_O_B,fs); 
 audiowrite([fout '_U_D.wav'],Y_U_D,fs);  
  
 audiowrite([fout '_omlsa.wav'],yout_omlsa,fs); 
 
fprintf('End of DUAL GSC\n');




