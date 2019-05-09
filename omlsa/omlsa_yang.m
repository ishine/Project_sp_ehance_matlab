%% Based on 'Speech enhancement for non-stationary noise enviroments' Isral Cohen. 2001 Elsevier.
%  1) 'om-lsa' to estimate speech present probability
%  2) 'mcra' to estimate noise spectral 


%fin0 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM1';
%fin1 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM2';
fin0 = 'F:/Work/2018/Beamforming/matlab/GSCLMS/voice/SlM1';

fOut = [fin0 '_OMLSA_YANG_OUT'];
 
[x1,fs1]= audioread([fin0 '.wav']); % main mic
%[x2,fs2]= audioread([fin1 '.wav']); % ref mic

Lens = length(x1);
FrameLen = 512;
FrameShift = FrameLen/2;
FFT_LEN=2^nextpow2(FrameLen);
FrameNum = Lens / FrameShift; 
M = 512; 
M21=M/2+1;

ex_gmma = ones(M21,1);

%% Noise Spectram Extimation
alpha_s = 0.8;  S = zeros(M21,1);
alpha_p1 = 0.2;  P1_est = zeros(M21,1);
delta = 5;       I = zeros(M21,1);
alpha_d = 0.95;  alpha_d1 = zeros(M21,1); lamda_d = zeros(M21 , 1);
 
%% Gain Computation
alpha_eta = 0.95; eta = zeros(M21,1);
beta = 0.7; xi = zeros(M21, 1);
qmax = 0.95;
xi_min_dB = -10; xi_max_dB = -5;  
xi_p_min_dB = 0; xi_p_max_dB = 10; P_min = 0.005;
w = 1; w_xi_local = 1; w_xi_global = 15;
GH1 =  ones(M21, 1);
G = ones(M21, 1);Gmin = 10^-25;

P =zeros(M21,1);
X_G = zeros(FFT_LEN,1);
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
x_out = zeros(Lens,1);

 init_frame = 1;
 end_frame =  init_frame+FrameLen-1;

while(end_frame<Lens)

    % win and fft
    x_frame = x1(init_frame:end_frame);    
    X_F = fft(window.* x_frame);     
    
    X_F_2 = abs(X_F(1:M21)).^2; 
    
        if FrameCnt==1     % new version omlsa3
            lamda_d=X_F_2;
        end    
    
    gmma = X_F_2 ./lamda_d;
    eta     = alpha_eta * (GH1.^2) .* ex_gmma + (1-alpha_eta) * max(gmma-1,0);
    ex_gmma = gmma;    
    v  = gmma .* eta ./ (1+eta);    
   %% Noise Spectram Extimation

    Sf = X_F_2;
  
    % find the min among L frame  i
    if(FrameCnt<15)
        S = Sf;
        Smin = S;
        Stmp = S;
        
    else
        S = alpha_s * S + (1-alpha_s * Sf);
        Smin = min(Smin, S);
        Stmp = min(Stmp, S);
    end
     
    Sr = S ./ Smin;
    
    for i=1:M21
       if(Sr(i) > delta) 
           I(i) = 1;
       else
           I(i) = 0;
       end
    end
    
    P1_est = alpha_p1 * P1_est + (1 - alpha_p1) * I;    
    alpha_d1 = alpha_d + (1-alpha_d) * P1_est;
  
    lamda_d = alpha_d1 .* lamda_d + (1- alpha_d1) .* X_F_2;

    
   %% Gain Computation
   
   % (1) est q     
    xi = beta .* xi + (1-beta) .* eta;

    xi_local=conv(xi,b_xi_local);
    xi_local=xi_local(w_xi_local+1:M21+w_xi_local);
    xi_global=conv(xi,b_xi_global);
    xi_global=xi_global(w_xi_global+1:M21+w_xi_global);
    xi_frame = mean(xi(1:M/2+1));
    if xi_local>0,  xi_local_dB=10*log10(xi_local)  ; else xi_local_dB=-100;  end
    if xi_global>0, xi_global_dB=10*log10(xi_global); else xi_global_dB=-100; end
    if xi_frame>0,  xi_frame_dB=10*log10(xi_frame)  ; else xi_frame_dB=-100;  end
      
    
    % matlab中对数组元素的条件赋值应该用如下方式：
     P_local = ones(M21,1);
     P_local(xi_local_dB<=xi_min_dB) = P_min;% use 0.005 instead of 0, costrain the value
     idx=find(xi_local_dB>xi_min_dB & xi_local_dB<xi_max_dB);
     P_local(idx) = P_min+(xi_local_dB(idx)-xi_min_dB)/(xi_max_dB-xi_min_dB)*(1-P_min); 
    
     P_global = ones(M21,1);
     P_global(xi_global_dB <= xi_min_dB) = P_min;
     idx = find(xi_global_dB > xi_min_dB & xi_global_dB < xi_max_dB);
     P_global(idx) = P_min + (xi_global_dB(idx) - xi_min_dB)/ (xi_max_dB - xi_min_dB) * (1-P_min);
      
         if xi_frame_dB<=xi_min_dB
            P_frame=P_min;
        elseif xi_frame_dB > ex_xi_frame_dB
            xi_peak_dB=min(max(xi_frame_dB,xi_p_min_dB),xi_p_max_dB);
            P_frame=1;
        elseif xi_frame_dB>=xi_peak_dB+xi_max_dB %  'val' -> 'log' '*' to  '+'
            P_frame=1;
        elseif xi_frame_dB<=xi_peak_dB+xi_min_dB
            P_frame=P_min;
        else
            P_frame= (xi_frame_dB-xi_min_dB-xi_peak_dB)/(xi_max_dB-xi_min_dB);   
            %   P_frame=P_min+(xi_frame_dB-xi_m_dB-xi_fl_dB)/(xi_fu_dB-xi_fl_dB)*(1-P_min);
         end
     
        ex_xi_frame_dB = xi_frame_dB;
      
    q = 1 - P_frame .* P_global .* P_local;     
    q = min(qmax, q);
    
    % (2) est p gmma eta

     
    P = 1 + (q./(1-q)) .* (1+eta) .* exp(-v);
    P = 1./P;
    
    GH1 = (eta./(1+eta));
    tmp = exp(0.5*expint(v));
    GH1 = GH1.* tmp;
    
    
    G = GH1.^P .* (Gmin.^(1-P));
   % G = min(G, Gmax);
    
   %% Output
    % ifft & overlap & add
    X_G(1:M21) = X_F(1:M21).* G;    
    
    X_G(M21+1:FFT_LEN) = conj(X_G(M21-1:-1:2));
    
    x_i = ifft(X_G);    
    x_out(FrameCnt*FrameShift+1:(FrameCnt+1)*FrameShift) = history + x_i(1:FrameShift);     
    history = x_i(FrameShift+1:FFT_LEN);
     
    init_frame = init_frame + FrameShift;
    end_frame  = end_frame + FrameShift;
    FrameCnt = FrameCnt+1;    
end

 audiowrite([fOut '.wav'],x_out,fs1);

