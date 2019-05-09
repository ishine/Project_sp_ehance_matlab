close all; clc; clear all;

fprintf('test dual posterfilter \n');

finl = './voice/t15_90_gsc3_O_B';
%finr = './voice/t14r';
fout = [finl '_domtest'];

%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
 
  
[x,fs]= audioread([finl '.wav']);  % main mic
 
x1 = x(:,1);
x2 = x(:,2);

lenS = length(x1);


Nomlsa = 512;
Nomlsa41 = Nomlsa/4;

frame_head = 1;
frame_end = Nomlsa;
Y_om_in = zeros(Nomlsa,1);
U_om_in = zeros(Nomlsa,1); 
 global k_omlsa ;
 
 k_omlsa = 0;

yout_omlsa = zeros(lenS,1);

while(frame_end<lenS)

    Y_om_in = x1(frame_head:frame_end);
    U_om_in = x2(frame_head:frame_end); 
     
    yout_omlsa_t  =  tf_gsc_dual_postfilter(Y_om_in,U_om_in,k_omlsa+1 );  
    yout_omlsa(k_omlsa*Nomlsa41+1 :(k_omlsa+1)* Nomlsa41) = yout_omlsa_t;
    k_omlsa = k_omlsa +1;
    
frame_head= frame_head+ Nomlsa41;
frame_end = frame_end+ Nomlsa41;
end
om_cnt = 0;   

 out  = [yout_omlsa, x1  ];

 audiowrite([fout,'.wav'] ,out,fs); 


