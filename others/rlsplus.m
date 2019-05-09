 

fin0 = '..\voice\T13L';
fin1 = '..\voice\T13R';
fout = [fin0 'NLSM_OUT_t'];
 
    
 
[mixed,fs1]   = audioread([fin0 '.wav']);  % main mic
[ref_noise,fs2] = audioread([fin1 '.wav']); % ref mic
                 
                       
%clean=s';
  
mu=1.05;M=64;espon=1e-4;
% [en,wn,yn]=lmsFunc(mu,M,ref_noise,mixed);
% [en,wn,yn]=nlmsFunc(mu,M,ref_noise,mixed,espon);
delta = 1e-7;
lambda = 1;
%[en,w]=rls(lambda,M,ref_noise,mixed,delta);

w=zeros(M,1);
P=eye(M)/delta;
u=ref_noise(:);
d=mixed(:);
% input signal length
N=length(u);

eout = zeros(N,1);

% error vector
e=d.';
% Step2: Loop, RLS
for n=M:N
    uvec=u(n:-1:n-M+1);
    e(n)=d(n)-w'*uvec;
    k=lambda^(-1)*P*uvec/(1+lambda^(-1)*uvec'*P*uvec);
    P=lambda^(-1)*P-lambda^(-1)*k*uvec'*P;
    w=w+k*conj(e(n));
end


 audiowrite([fout '_rls_e.wav'],e,fs);
