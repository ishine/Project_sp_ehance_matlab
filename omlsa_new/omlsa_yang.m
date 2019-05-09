%% Based on 'Speech enhancement for non-stationary noise enviroments' Isral Cohen. 2001 Elsevier.
%  1) 'om-lsa' to estimate speech present probability
%  2) 'mcra' to estimate noise spectral 

clear all;
%fin0 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM1';
%fin1 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM2';
fin0 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/noise2';

fOut = [fin0 '_OMLSA_YANG_OUT'];
 
[x1,Fs]= audioread([fin0 '.wav']); % main mic
%[x2,fs2]= audioread([fin1 '.wav']); % ref mic

Lens = length(x1);
FrameLen = 512;
FrameShift = FrameLen/4;
FFT_LEN=2^nextpow2(FrameLen);
FrameNum = Lens / FrameShift; 
M = FrameLen; Mo = M * 0.75;
M21=M/2+1;
gmma= ones(M21,1);
ex_gmma = ones(M21,1);

% window function
win=hamming(M);
% find a normalization factor for the window
win2=win.^2;
Mno=M-Mo;

out = zeros(M,1);

W0=win2(1:Mno);
for k=Mno:Mno:M-1
    swin2=lnshiftyang(win2,k);
    W0=W0+swin2(1:Mno);
end
W0=mean(W0)^0.5;
win=win/W0;
Cwin=sum(win.^2)^0.5;

%% Noise Spectram Extimation
w=1;
alpha_s = 0.9;  S = zeros(M21,1);  Smin  = zeros(M21,1);  Stmp = zeros(M21,1);  S_his = zeros(M21,8); FrameMode = 0; 
alpha_p1 = 0.2;  P1_est = zeros(M21,1);
delta = 5;       I = zeros(M21,1);
alpha_d = 0.95;  alpha_d1 = zeros(M21,1); lambda_d = zeros(M21 , 1); lamda_d_long = zeros(M21,1);
alpha_d_long = 0.99; lambda_d_long = zeros(M21,1);
k2_local=round(500/Fs*M+1); % 500hz
k3_local=round(3500/Fs*M+1); %3500hz
tone_flag = 1; 

delta_s=1.67;		% 2.4)  Local minimum factor
Bmin=1.66;
delta_X2=4.6;		% 2.4)  Local minimum factor
delta_yt=3;
%% Gain Computation
alpha_eta = 0.95; eta = zeros(M21,1);
beta = 0.7; xi = zeros(M21,1);
qmax = 0.95;
xi_min_dB = -10; xi_max_dB = -5;  
xi_p_min_dB = 0; xi_p_max_dB = 10; P_min = 0.005;
w = 1; w_xi_local = 1; w_xi_global = 15;
GH1 =  ones(M21, 1);
G = ones(M21, 1); 
eta_min_dB=-18;	
eta_min=10^(eta_min_dB/10);
Gmin=eta_min^0.5;	   % Gain floor
P =zeros(M21,1);
X_G = zeros(FFT_LEN,1);
eta_2term = 1;
xi_frame = 0;
%% parameter

window = hanning(FrameLen);
b=hanning(2*w+1);
b=b/sum(b);     % normalize the window function
b_xi_local=hanning(2*w_xi_local+1);
b_xi_local=b_xi_local/sum(b_xi_local);  % normalize the window function
b_xi_global=hanning(2*w_xi_global+1);
b_xi_global=b_xi_global/sum(b_xi_global);   % normalize the window function
 
FrameCnt = 1;
x_frame = zeros(FrameLen,1);
history = zeros(FrameShift,1);
x_out = zeros(Lens,1,'int16');

 init_frame = 1;
 end_frame =  init_frame+FrameLen-1;
 
 
 zz = expint(0.3318);
 
while(end_frame<Lens)

    % win and fft
    x_frame = x1(init_frame:end_frame);    
 %   X_F = fft(window.* x_frame);  
    x_frame =win.* x_frame;
    X_F = fft( x_frame);  
    
    X_F_2 = abs(X_F(1:M21)).^2; 
    
    if FrameCnt==1     % new version omlsa3
         lambda_d=X_F_2;
    end
    
    gmma=X_F_2./max(lambda_d,1e-10);
    eta     = alpha_eta * (GH1.^2) .* ex_gmma + (1-alpha_eta) * max(gmma-1,0);
    ex_gmma = gmma;    
    v  = gmma .* eta ./ (1+eta);    
   %% Noise Spectram Extimation

   % Sf = X_F_2;
    Sf=conv(b,X_F_2);  % smooth over frequency
    Sf=Sf(w+1:M21+w);
    
     if FrameCnt==1     % new version omlsa3
            Sy=X_F_2;
            S=Sf;
            St=Sf;
            lambda_d=X_F_2;      
     else
            S=alpha_s*S+(1-alpha_s)*Sf;     % smooth over time      
     end    
    
    % find the min among L frame  i
        if(FrameCnt<15)                      
            Smin = S;               
        else                
            Smin = min(Smin, S);                   
        end
        
        
      Sft = St;
    
      I = double(X_F_2 < delta_X2 * Bmin.*Smin  &  S < delta_s * Bmin .* Smin);
      conv_I = conv(b, I);
      conv_I = conv_I(w+1:M21+w);
             
      idx = find(conv_I);
    
      if ~isempty(idx)        
        conv_Y=conv(b,I.*X_F_2);
        conv_Y=conv_Y(w+1:M21+w);
        Sft(idx)=conv_Y(idx)./conv_I(idx);                            
      end        
 
  
        if FrameCnt<14+1     % new version omlsa3
            St=S;
            Smint=St;
            SMactt=St;
        else
            St=alpha_s*St+(1-alpha_s)*Sft;
            Smint=min(Smint,St);
            SMactt=min(SMactt,St);
        end
                
        qhat=ones(M21,1);    
        phat=zeros(M21,1);
         
        Sr_X2 = X_F_2./Bmin./max(Smint,1e-10);   % like the Sr for this      
        Sr_S = S./Bmin./max(Smint,1e-10);      % new version
    
        idx=find(Sr_X2>1 & Sr_X2<delta_yt & Sr_S<delta_s);
        qhat(idx)=(delta_yt-Sr_X2(idx))/(delta_yt-1);
        phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));
        phat(Sr_X2>=delta_yt | Sr_S>=delta_s)=1; % both '>' , set 1
      
   % P1_est = alpha_p1 * P1_est + (1 - alpha_p1) * I;    
   % alpha_dt=alpha_d+(1-alpha_d)*phat;                 % eq 31
   % lambda_dav=alpha_dt.*lambda_dav+(1-alpha_dt).*Ya2; % eq 30
  
      alpha_dt = alpha_d + (1-alpha_d) * phat;  
      lambda_d = alpha_dt .* lambda_d + (1- alpha_dt) .* X_F_2;
 
       if FrameCnt<15     % new version omlsa3
            lambda_d_long=lambda_d;
        else
            alpha_dt_long=alpha_d_long+(1-alpha_d_long)*phat;
            lambda_d_long=alpha_dt_long.*lambda_d_long+(1-alpha_dt_long).*X_F_2;
       end    
    
      % update Smin Stmp S_his
      FrameMode = FrameMode+1;      
      if(  FrameMode == 16)              
          FrameMode = 0;          
          if(FrameCnt==15)                         
               S_his=repmat(SMactt,1,8); 
          else               
               S_his=[S_his(:,2:8) SMactt];% SWt=[SWt(:,2:Nwin) SMactt];         
               Smint = min(S_his,[],2);              
               SMactt = S;
          end               
      end
      
      lambda_d1=1.4685*lambda_d ;
       
   %% Gain Computation
   
   % (1) est speech absense probability 'q'     
    xi = beta .* xi + (1-beta) .* eta;       % eq. (23)           
    xi_local=conv(xi,b_xi_local);            % eq. (24)      
    xi_local=xi_local(w_xi_local+1:M21+w_xi_local);
    xi_global=conv(xi,b_xi_global);          % eq. (24)      
    xi_global=xi_global(w_xi_global+1:M21+w_xi_global);
    
    ex_xi_frame = xi_frame;
    xi_frame = mean(xi(3:M21));
    if xi_local>0,  xi_local_dB=10*log10(xi_local)  ; else xi_local_dB=-100;  end
    if xi_global>0, xi_global_dB=10*log10(xi_global); else xi_global_dB=-100; end
    if xi_frame>0,  xi_frame_dB=10*log10(xi_frame)  ; else xi_frame_dB=-100;  end
      
  
    %  P_local 
     P_local = ones(M21,1);
     P_local(xi_local_dB<=xi_min_dB) = P_min;% use 0.005 instead of 0, costrain the value
     idx=find(xi_local_dB>xi_min_dB & xi_local_dB<xi_max_dB);
     P_local(idx) = P_min+(xi_local_dB(idx)-xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
    
     % P_global
     P_global = ones(M21,1);
     P_global(xi_global_dB <= xi_min_dB) = P_min;
     idx = find(xi_global_dB > xi_min_dB & xi_global_dB < xi_max_dB);
     P_global(idx) = P_min + (xi_global_dB(idx) - xi_min_dB)/ (xi_max_dB - xi_min_dB) * (1-P_min);
     
        m_P_local=mean(P_local(3:(k2_local+k3_local-3)));    % average probability of speech presence 500 ~3500hz
    
        if m_P_local<0.25
            P_local(k2_local:k3_local)=P_min;    % reset P_local (frequency>500Hz) for low probability of speech presence
        end         
    
       if tone_flag               % new version
           if (m_P_local<0.5) && (FrameCnt>120)
               idx=find( alpha_dt_long(8:(M21-8)) > 2.5*(alpha_dt_long(10:(M21-6))+alpha_dt_long(6:(M21-10))) );
               P_local([idx+6;idx+7;idx+8])=P_min;   % remove interfering tonals
           end
       end           
     
     % P_frame  Fig.3  eq.26 & eq.27
     if xi_frame_dB<=xi_min_dB     
         P_frame=P_min;        
     elseif xi_frame >= ex_xi_frame             
         xi_peak_dB=min(max(xi_frame_dB,xi_p_min_dB),xi_p_max_dB);            
         P_frame=1;        
     elseif xi_frame_dB>=xi_peak_dB+xi_max_dB % eq.(27)  u = 1  xi_frame > xi_peak * xi_max
         P_frame=1;        
     elseif xi_frame_dB<=xi_peak_dB+xi_min_dB  % eq.(27) u = 0  xi_frame < xi_peak * xi_min         
         P_frame=P_min;        
     else
          %  P_frame= (xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB);             
          P_frame=P_min+(xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB)*(1-P_min);         
     end      
      
      
    q = 1 - P_frame .* P_global .* P_local;     
    q = min(qmax, q);
     
    % (2) est speech present probability 'P' and gain

        gmma=X_F_2./max(lambda_d1,1e-10);                       % postier SNR eq.10 in paper
        eta=alpha_eta* eta_2term +(1-alpha_eta)*max(gmma-1,0); % proior SNR  eq.28  
        eta=max(eta,eta_min);
        v=gmma.*eta./(1+eta);                                % eq.10       
        
    % gmma=X_F_2./max(lambda_d,1e-10);
    % eta     = alpha_eta * eta_2term + (1-alpha_eta) * max(gmma-1,0);
    % ex_gmma = gmma;    
    % v  = gmma .* eta ./ (1+eta);         
        
     
     PH1=zeros(M21,1);
     idx=find(q<0.9);
     PH1(idx)=1./(1+q(idx)./(1-q(idx)).*(1+eta(idx)).*exp(-v(idx)));
    
    % Gain 
     GH1=ones(M21,1);
     idx=find(v>5);
     GH1(idx) = (eta(idx)./(1+eta(idx)));
     idx=find(v<=5 & v>0);
     GH1(idx)=eta(idx)./(1+eta(idx)).*exp(0.5*expint(v(idx)));
  
    % G = GH1.^P .* (Gmin.^(1-P));
         
     if tone_flag   % new version     
         lambda_d_global=lambda_d1;   % new version         
         lambda_d_global(4:M21-3)=min([lambda_d_global(4:M21-3),lambda_d_global(1:M21-6),lambda_d_global(7:M21)],[],2);   % new version         
         Sy=0.8*Sy+0.2*X_F_2;    % new version
         GH0=Gmin*(lambda_d_global./(Sy+1e-10)).^0.5;   % new version omlsa3        
     else   % new version            
         GH0=Gmin;   %#ok<UNRCH> % new version   
     end   % new version
        
     G=GH1.^PH1.*GH0.^(1-PH1);     
     eta_2term=GH1.^2.*gmma;
  
   %% Output
    % ifft & overlap & add
    X_G(1:M21) = X_F(1:M21).* G;        
    X_G(M21+1:FFT_LEN) = conj(X_G(M21-1:-1:2));
     
 %   x_i = real(ifft(X_G));    
  %  x_out(FrameCnt*FrameShift+1:(FrameCnt+1)*FrameShift) = (history + x_i(1:FrameShift)).*2^15;     
  %  history = x_i(FrameShift+1:FFT_LEN);
          
  x_i= win.*real(ifft(X_G));
  %     x_i= win.*real(ifft(X_F));
    out=out+x_i;        
    x_out((FrameCnt-1) * Mno+1 : (FrameCnt) * Mno  ) = out(1:Mno) .*2^15;    
    
      out=[out(Mno+1:M); zeros(Mno,1)]; 
        
    init_frame = init_frame + FrameShift;
    end_frame  = end_frame + FrameShift;
    FrameCnt = FrameCnt+1;    
end

 audiowrite([fOut '.wav'],x_out,Fs);

