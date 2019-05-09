function NLMS_func( fin0 ,  d, x   )
%fin0 = '.\voice\t10l';
%fin1 = '.\voice\t10r';
  
fout  = [fin0 '_NLMS_v2'];
fout1 = [fin0 '_B'];
Lw = 256;
%[Y_L,fs1]= audioread([fin0 '.wav']); % main mic
%[Y_R,fs2]= audioread([fin1 '.wav']); % ref mic

%d = (Y_L + Y_R )*0.5;
%x = Y_L - Y_R;
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
            m=1;
        elseif mu<0
            mu=0;
        end
    end
    w=w+mu/(xtdl'*xtdl+psi)*xtdl*e(n);
 %   zeta(n)=(w-wo)'*(w-wo);
end
 audiowrite([fout '_Out.wav'],e,fs);
 audiowrite([fout '_d.wav'],d,fs);
 audiowrite([fout '_x.wav'],x,fs);
 
 
 fprintf('End of nlmsv2\n');
 
gsc_dual_postfilter([fout '_x'],  [fout '_Out']);
 %  gsc_dual_postfilter('./voice/t10l_lcmv',  [fout '_x']);

  fprintf('End of dual channel omlsa\n');

 
 