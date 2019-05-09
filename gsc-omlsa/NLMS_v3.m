
fin0 = '.\voice\NLMS';
  
fout  = [fin0 '_NLMS_v2'];
fout1 = [fin0 '_B'];
Lw = 256;
%[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
%[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic
[Y ,fs1] =  audioread([fin0 '.wav']); 
 
d =  Y(:,1);
x =  Y(:,2);
% [error,filter_coeff ] = FLMS(alpha,filter_length,input_signal,desired_signal,lambda,estimated_power,imp_response); 
estimated_power = zeros(1024,1);

[error,~ ]= FLMS(0.01, 512,x, d , 0.985,estimated_power   );

errorres = [zeros(length(x) -  length(error),1) ;error ];

 cmpf = [ errorres, x];
 
 audiowrite([fout '_flmscmp.wav'],cmpf,fs);
 
 audiowrite([fout 'err.wav'],errorres,fs);
 audiowrite([fout 'x.wav'],x,fs);
 
 
% gsc_dual_postfilter(  [fout 'x'], [fout 'err']);
 
 
 
 
 