
%  目前该版本为 GSC 算法， 基于mtlab 自身模拟的阵列信号8khz/16khz下能实现较好的效果。 
%  采用2通道  16Tap/8Tap 的LCMV算法。

close all; clc; clear all;

fprintf('dual channel gsc + posterfilter \n');

%    -3  -2   -1  0  1  2  3 

Nmic = 2; % mic number 
Nflt = 256; % fir length
Num_c = Nmic*Nflt;
Dmic = 0.08;
Speed = 340;
 
finl = './voice/t12l';
finr = './voice/t12r';


%  142  122  105  90,74 57 37    
Angle = 90; % angle of signal
fout = [finl '_90_gscv2'];
 
 
[x1,fs1]= audioread([finl '.wav']);  % main mic
[x2,fs2]= audioread([finr '.wav']); % ref mic
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
if(taulen<0)  % delay x1  
   delayPoint = round(-taulen);     
   xTwo = x2(1:lenS);  
   xOne =  zeros(lenS,1);
   xOne(delayPoint+1:lenS) = x1(1:lenS-delayPoint);
else          % delay x2
   delayPoint = round(taulen);
   xOne = x1(1:lenS);   
   xTwo =  zeros(lenS,1);
   xTwo(delayPoint+1:lenS) = x2(1:lenS-delayPoint);   
end
  

    f_st = zeros(Nflt,1);
    
    f_st(1:4) =  [1,0,0,0];
    
 
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
 x = xOne - xTwo;
 % up line of GSC beamformer 
 d  =  (xOne+xTwo)  /2; % Y_Up 
 
  
 B = zeros(Nmic-1, Nmic);
for i = 1 : Nmic - 1
    B(i, i) = 1;
    B(i, i+1) = -1;
end 
 %mu =0.001;
 
 
 Lw = 256;
 
  e = zeros(length(x), 1);

  w=zeros(Lw,1);
  
mu=0.5;psi=0.1;alpha=0.995;
Pd=1; Pyhat=1; Pe=1;
flag=1;

for n=length(w):length(x)
    xtdl=x(n:-1:n-Lw+1);
    yhat=w'*xtdl;
    e(n)=d(n)-yhat;
 %   r(n)=y(n)-yhat;
    if flag==1
        Pd=alpha*Pd+(1-alpha)*d(n)*d(n);
        Pyhat=alpha*Pyhat+(1-alpha)*yhat*yhat;
        Pe=alpha*Pe+(1-alpha)*e(n)*e(n);
        mu=1-0.5*Pyhat/Pd;
        if mu>1
            mu=1;
        elseif mu<0
            mu=0;
        end
    end
    w=w+mu/(xtdl'*xtdl+psi)*xtdl*e(n);
 %   zeta(n)=(w-wo)'*(w-wo);
end
   
audiowrite([fout '_OUT.wav'],e,fs);
audiowrite([fout '_UP.wav'],d,fs); 
audiowrite([fout '_B.wav'],x,fs);

fprintf('2nd do posterfilter \n');
 
 gsc_dual_postfilter([fout '_UP'],  [fout '_B']);
 
 gsc_dual_postfilterInv(   [fout '_B'], [fout '_OUT']);

 % omlsa(fOut);
fprintf('End of DUAL GSC\n');
