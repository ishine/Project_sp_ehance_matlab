
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

[error,yt ]= FLMS(0.0988, 512,x, d , 0.95,  estimated_power   );

 errorres = [zeros(length(x) -  length(error),1) ;error ];
 yt = yt';
 cmpf = [  x, errorres ];
 
 audiowrite([fout '_flmscmp.wav'],cmpf,fs);
 
 audiowrite([fout 'err.wav'],errorres,fs);
 audiowrite([fout 'x.wav'],x,fs);
 
 
gsc_dual_postfilter(  [fout 'x'], [fout 'err']);
  %gsc_dual_postfilter(  [fout 'err'],[fout 'x']);
 


e = zeros(length(x), 1);

w=zeros(Lw,1);
  
mu=0.098;psi=0.001;alpha=0.995;
Pd=1; Pyhat=1; Pe=1;
flag=1;
y_t = zeros(length(x),1);

for n=length(w):length(x)
    xtdl=x(n:-1:n-Lw+1);
    yhat=w'*xtdl;  y_t(n) = yhat;
    e(n)=d(n)-yhat;
      
    w=w+mu*e(n)*xtdl/(xtdl'*xtdl+psi) ;
 
end

% e_t = Y_L - y_t;
 
 cmp = [e,x];

 audiowrite([fout '_Out.wav'],e,fs);
 audiowrite([fout '_y.wav'],y_t,fs);
  audiowrite([fout 'cmp.wav'],cmp,fs); 
 
 fprintf('End of nlmsv2\n');
 
 
 
 