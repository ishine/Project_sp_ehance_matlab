
fin0 = '.\voice\NLMS';
  
fout  = [fin0 '_ccflms_'];
fout1 = [fin0 '_B'];
Lw = 256;
%[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
%[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic
[Y ,fs1] =  audioread([fin0 '.wav']); 
 
d =  Y(:,1);
x =  Y(:,2);

d1= d.';
x1= x.';
[e] = cc2_flms(x1,d1,0.01,512);

 audiowrite([fout 'ek2.wav'],e,fs);


% nlms
NFFT = 1024;  fshift = NFFT/2;
fhead = 1; fend = NFFT;

lenS = length (d);
l = 1;
mu = 0.01;
P = zeros(NFFT,1) + 1e-10;
Wk = zeros(NFFT,1);
lambda = 0.8;
histykt  = zeros(fshift,1);  histekt  = zeros(fshift,1);
win = hanning(NFFT);

y_t =zeros(lenS,1);
e_t = zeros(lenS,1);


while(fend < lenS)
    
    d_k = d(fhead:fend);
    x_k = x(fhead:fend);
    
    Dk = fft(win.*d_k);
    Xk = fft(win.*x_k);
    
    Yk = Wk .* Xk;  
    
    Ek = Dk - Yk;
    
    P = lambda * P + (1-lambda) * abs(Dk).^2;
    
    Muk = mu ./(P + 0.00001);
    
    Wk = Wk + mu.*conj(Xk).*Ek;
    
    yk_t = ifft(Yk);
    
    ek_t = ifft(Ek);
        
    y_t( (l-1)*fshift+1: l*fshift  ) =   yk_t(1:fshift) + histykt;    
    histykt = yk_t(fshift+ 1 : NFFT) ;

    
    e_t( (l-1)*fshift+1: l*fshift  ) =   ek_t(1:fshift) + histekt;    
    histekt = ek_t(fshift+ 1 : NFFT) ;
    
    
    l = l +1;
    fhead = fhead + fshift;
    fend  = fend  + fshift;
 
end
 
 out = [y_t, d] ;
 audiowrite([fout '.wav'],out,fs);
 out = [e_t, d] ;
 audiowrite([fout 'ek.wav'],out,fs);
 
 fprintf('End of ccflms\n');
 
 function [e] = cc2_flms(x,d,u,N)
 
 M = length(x);
 
 h = zeros(1,N);
 e = zeros(1,M);
 y = zeros(1,M);
 for i = 1: fix(M/N)-1
     
     X = fft(x((i-1)*N+1: (i+1)*N));
     H =fft([h,zeros(1,N)]);
     O1 = real(ifft(X.*H));
     
     y(i*N+1:(i+1)*N) = O1(N+1:2*N);
     e(i*N+1:(i+1)*N) = d(i*N+1:(i+1)*N) - y(i*N+1:(i+1)*N);
     E = fft([zeros(1,N), e(i*N+1:(i+1)*N)]);
     O2 =real(ifft(E.*conj(X)));
     V = O2(1:N);
     h = h+2*u*V;
     
 end
 
 end
 
 
 
 
 
 
 
 
 
 
 