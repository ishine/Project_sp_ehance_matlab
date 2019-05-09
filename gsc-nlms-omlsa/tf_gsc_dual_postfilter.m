%function [y,out]=tf_gsc_dual_postfilter(fin, fin2 )
function [out_om]=tf_gsc_dual_postfilter(y_frame_in , u_frame_in , initflag)

% omlsa : Single Channel OM-LSA with IMCRA noise estimator
% ***************************************************************@
% Inputs:
%    fin,  input file name (fin.wav)
%    fout, output file name (fout.wav)
% Output:
%    in,  samples of the input file
%    out, samples of the output file
% Usage:
%    [in, out]=omlsa(fin,fout);
%    omlsa(fin);
% Defaults:
%    fout= [fin,'_omlsa'];
%
% Copyright (c) 2003. Prof Israel Cohen.
% All rights reserved. 
% ***************************************************************@
 
%fin =  '..\..\voice\T7L1142_OUT';
%fin2 = '..\..\voice\T7L1142_B';
nin=nargin;
global k_omlsa;

global SPP  q_hat;
global Fs_ref;
global M_ref ; 
global Mo_ref;
global w;
global alpha_s_ref; global Nwin; global Vwin; global delta_s; global Bmin;
global delta_y; global delta_yt; global alpha_d_ref;
global alpha_xi_ref; global w_xi_local; global w_xi_global; global f_u;
global f_l; global P_min; global xi_lu_dB; global xi_ll_dB; global xi_gu_dB;
global xi_gl_dB; global xi_fu_dB ; global xi_fl_dB ; global xi_mu_dB ; global xi_ml_dB  ; global q_max;
global alpha_eta_ref; global eta_min_dB;
global broad_flag; global tone_flag; global nonstat;


  global S;
  global Sf; global Sf_u; global Sy; global Su; global lambda_dav;
  global lambda_d_hat;  global St; global S_u; global St_u; global lambda_dav_u;
  global lambda_d;  global lambda_d_u;
  global eta_2term_y; global eta_2term_u ;global  eta_2term ;
  global eta; global eta_u;
  global fnz_flag;  
  global Smin; global SMact; global Smin_u; global SMact_u;
  global Sft; global Sft_u;  
  global Smint;   global SMactt;
  global Smint_u; global SMactt_u;    
  global win; global Cwin;  
 global l_fnz; global l_omlsa;
     global lambda_dav_long; global lambda_dav_long_u;  global l_mod_lswitch;
  global Mno; global out_omlsa;  global lframe ;
  
  global SW; global  SWt; global SW_u; global SWt_u;
  
% 1) Parameters of Short Time Fourier Analysis:
Fs_ref=16e3;		% 1.1) Reference Sampling frequency
M_ref=512;		% 1.2) Size of analysis window
Mo_ref=0.75*M_ref;	% 1.3) Number of overlapping samples in consecutive frames

% 2) Parameters of Noise Spectrum Estimate
w=1;			% 2.1)  Size of frequency smoothing window function=2*w+1
alpha_s_ref=0.9;	% 2.2)  Recursive averaging parameter for the smoothing operation
Nwin=8; 	% 2.3)  Resolution of local minima search
Vwin=15;
delta_s=1.67;		% 2.4)  Local minimum factor
Bmin=1.66;
delta_y=4.6;		% 2.4)  Local minimum factor
delta_yt=3;
alpha_d_ref=0.85;	% 2.7)  Recursive averaging parameter for the noise


% 3) Parameters of a Priori Probability for Signal-Absence Estimate
alpha_xi_ref=0.7;	% 3.1) Recursive averaging parameter
w_xi_local=1; 	% 3.2) Size of frequency local smoothing window function
w_xi_global=15; 	% 3.3) Size of frequency local smoothing window function
f_u=10e3; 		% 3.4) Upper frequency threshold for global decision
f_l=50; 		% 3.5) Lower frequency threshold for global decision
P_min=0.005; 		% 3.6) Lower bound constraint
xi_lu_dB=-5; 	% 3.7) Upper threshold for local decision
xi_ll_dB=-10; 	% 3.8) Lower threshold for local decision
xi_gu_dB=-5; 	% 3.9) Upper threshold for global decision
xi_gl_dB=-10; 	% 3.10) Lower threshold for global decision
xi_fu_dB=-5; 	% 3.11) Upper threshold for local decision
xi_fl_dB=-10; 	% 3.12) Lower threshold for local decision
xi_mu_dB=10; 	% 3.13) Upper threshold for xi_m
xi_ml_dB=0; 		% 3.14) Lower threshold for xi_m
q_max=0.998; 		% 3.15) Upper limit constraint



% 4) Parameters of "Decision-Directed" a Priori SNR Estimate
alpha_eta_ref=0.95;	% 4.1) Recursive averaging parameter
eta_min_dB=-18;	% 4.2) Lower limit constraint

% 5) Flags

broad_flag=1;               % broad band flag   % new version
tone_flag=1;                % pure tone flag   % new version
nonstat='medium';                %Non stationarity  % new version

 

% Adjust parameters according to the actual sampling frequency
Fs = Fs_ref;
if Fs~=Fs_ref
    M=2^round(log2(Fs/Fs_ref*M_ref));
    Mo=Mo_ref/M_ref*M;
    alpha_s=alpha_s_ref^(M_ref/M*Fs/Fs_ref);
    alpha_d=alpha_d_ref^(M_ref/M*Fs/Fs_ref);
    alpha_eta=alpha_eta_ref^(M_ref/M*Fs/Fs_ref);
    
else
    M=M_ref;
    Mo=Mo_ref;
    alpha_s=alpha_s_ref;
    alpha_d=alpha_d_ref;
    alpha_eta=alpha_eta_ref;
 
end
alpha_d_long=0.99;
eta_min=10^(eta_min_dB/10);
G_f=eta_min^0.5;	   % Gain floor
  

out_om =zeros(Mo,1);
 if(initflag ==1)
% window function
win=hamming(M);
% find a normalization factor for the window
win2=win.^2;
Mno=M-Mo;

W0=win2(1:Mno);
for k=Mno:Mno:M-1
    swin2=lnshift(win2,k);
    W0=W0+swin2(1:Mno);
end

W0=mean(W0)^0.5;
win=win/W0;
Cwin=sum(win.^2)^0.5;
win=win/Cwin;
end

b=hanning(2*w+1);
b=b/sum(b);     % normalize the window function
b_xi_local=hanning(2*w_xi_local+1);
b_xi_local=b_xi_local/sum(b_xi_local);  % normalize the window function
b_xi_global=hanning(2*w_xi_global+1);
b_xi_global=b_xi_global/sum(b_xi_global);   % normalize the window function

M21=M/2+1;
k_u=round(f_u/Fs*M+1);  % Upper frequency bin for global decision
k_l=round(f_l/Fs*M+1);  % Lower frequency bin for global decision
k_u=min(k_u,M21);
k2_local=round(500/Fs*M+1);
k3_local=round(3500/Fs*M+1);
 
% 1) init for the 1st frame
if(initflag ==1) 
   l_fnz=1; l_omlsa = 1;  fnz_flag=0;   lframe = 1;
   eta_2term=ones(M21,1);
   eta_2term_y = ones(M21,1); 
   eta_2term_u = ones(M21,1); 
   out_omlsa=zeros(M,1);     
   SPP = zeros(M21, 1);        
   l_mod_lswitch=0;   
end
 
     % counter for the first frame which is non-zero    % new version omlsa3
    % flag for the first frame which is non-zero    % new version omlsa3
zero_thres=1e-10;      % new version omlsa3
  
%while(EndFrame < len)
%if(EndFrame < len)
      y =  y_frame_in; % y = Y_BUF(InitFrame:EndFrame);
      u =  u_frame_in; % u = U_BUF(InitFrame:EndFrame);
       
    if (~fnz_flag && abs(y(1))>zero_thres) ||  (fnz_flag && any(abs(y)>zero_thres))       % new version omlsa3
        fnz_flag=1;     % new version omlsa3

     %% 1 prepare the freq data Y and U   
        Y=fft(win.*y);
        Ya2=abs(Y(1:M21)).^2;
        
        U = fft(win.*u);        
        Ua2 = abs(U(1:M21)).^2;
                 
     %% 2  Noise estimation using IMCRA     Ya2 Ua2: voice power    S  S_u: Smoothed voice power in freq and time  
        if l_omlsa==l_fnz     
            lambda_d=Ya2;  lambda_d_u=Ua2;             
        end               
        % (2.1) cal posteriori SNR ,prioir SNR , v 
        % gamma = input power / noise power; (noise power will be estimated)
        % eta = smoothed gamma -1. 
        gamma=Ya2./max(lambda_d,1e-10);
        eta=alpha_eta*   eta_2term_y +(1-alpha_eta)*max(gamma-1,0); % eta=alpha_eta*   eta_2term_y+(1-alpha_eta)*max(gamma-1,0);
        eta=max(eta,eta_min);
        v=gamma.*eta./(1+eta);
        % ref noise 
        gamma_u=Ua2./max(lambda_d_u,1e-10);
        eta_u=alpha_eta*eta_2term_u+(1-alpha_eta)*max(gamma_u-1,0);
        eta_u=max(eta_u,eta_min);
        v_u=gamma_u.*eta_u./(1+eta_u); % v_u=gamma.*eta_u./(1+eta_u);  20190507
        
        % (2.2) smooth Ya2 and Ua2, in freq and time: Ya2-> Sf -> S  Ua2 -> Sf_u->S_u
        %                                       init: Ya2 -> Sy lambda_dav lambda_d_dat
        Sf=conv(b,Ya2);    Sf=Sf(w+1:M21+w);      % smooth in freq
        Sf_u=conv(b,Ua2);  Sf_u=Sf_u(w+1:M21+w);
         
        if l_omlsa==l_fnz          % smooth in time
            Sy=Ya2;  S=Sf;  St=Sf; lambda_dav=Ya2;  lambda_d_hat =Ya2; 
            Su=Ua2;  S_u=Sf_u;  St_u=Sf_u; lambda_dav_u=Ua2;
        else
            S=alpha_s*S+(1-alpha_s)*Sf;      
            S_u=alpha_s*S_u+(1-alpha_s)*Sf_u;   
        end
        % (2.3) derive the minimum data of voice: S and S_u -> Smin Smin_u
        if l_omlsa<14+l_fnz      
            Smin=S;     SMact=S;
            Smin_u=S_u; SMact_u=S_u;
        else
            Smin=min(Smin,S);  SMact=min(SMact,S);
            Smin_u=min(Smin_u,S_u);  SMact_u=min(SMact_u,S_u);
        end
        % (2.4)  * find the freq bin may be noise ( both  Ya2 and S it smaller than threshold)
        %        * update the noise like freq in Sft     
        % Signal
        I_f=double(Ya2<delta_y*Bmin.*Smin & S<delta_s*Bmin.*Smin);                
        conv_I=conv(b,I_f);      conv_I=conv_I(w+1:M21+w);     Sft=St;        
        idx=find(conv_I);         
        if ~isempty(idx)
            if w
                conv_Y=conv(b,I_f.*Ya2);
                conv_Y=conv_Y(w+1:M21+w);
                Sft(idx)=conv_Y(idx)./conv_I(idx);  % only update noise possible freq bin in Sft
            else
                Sft(idx)=Ya2(idx);  
            end
        end
        
        
        % ref noise
        I_f_u=double(Ua2<delta_y*Bmin.*Smin_u & S_u<delta_s*Bmin.*Smin_u);
        conv_I_u=conv(b,I_f_u);  conv_I_u=conv_I_u(w+1:M21+w); Sft_u=St_u;        
        idx=find(conv_I_u);
        if ~isempty(idx)      
                conv_U=conv(b,I_f_u.*Ua2);  %  I_f  conv_U=conv(b,I_f.*Ua2); 
                conv_U=conv_U(w+1:M21+w);
                Sft_u(idx)=conv_U(idx)./conv_I_u(idx);         
        end     
        
        % (2.5) Smooth twice: Sft -> St -> Smint and Smactt.
        if l_omlsa<14+l_fnz      
            St=S;  Smint=St;  SMactt=St;
            St_u=S_u;  Smint_u=St_u;  SMactt_u=St_u;
        else
            St=alpha_s*St+(1-alpha_s)*Sft;  Smint=min(Smint,St); SMactt=min(SMactt,St);
            St_u=alpha_s*St_u+(1-alpha_s)*Sft_u;  Smint_u=min(Smint_u,St_u); SMactt_u=min(SMactt_u,St_u);
        end
        
        
        % (2.6) derive q and p again, like I,  lambda_dav lambda_dav_long
        qhat=ones(M21,1);  phat=zeros(M21,1);
        qhat_u=ones(M21,1);  phat_u=zeros(M21,1);
       
        switch nonstat    % new version
            case   'low'  % new version
                gamma_mint=Ya2./Bmin./max(Smin,1e-10);   % new version
                zetat=S./Bmin./max(Smin,1e-10);      % new version
                
                gamma_mint_u=Ua2./Bmin./max(Smin_u,1e-10);   % new version
                zetat_u=S_u./Bmin./max(Smin_u,1e-10);      % new version
            otherwise  % new version
                gamma_mint=Ya2./Bmin./max(Smint,1e-10);   % new version
                zetat=S./Bmin./max(Smint,1e-10);      % new version
                
                gamma_mint_u=Ua2./Bmin./max(Smint_u,1e-10);   % new version
                zetat_u=S_u./Bmin./max(Smint_u,1e-10);      % new version                
        end  % new version
        
        %  signal 
        idx=find(gamma_mint>1 & gamma_mint<delta_yt & zetat<delta_s);         % find noise like freq bin again
        qhat(idx)=(delta_yt-gamma_mint(idx))/(delta_yt-1);                    % qhat: probablity of voice absent
        phat(idx)=1./(1+qhat(idx)./(1-qhat(idx)).*(1+eta(idx)).*exp(-v(idx)));% phat: probablity of voice present
        phat(gamma_mint>=delta_yt | zetat>=delta_s)=1;
        alpha_dt=alpha_d+(1-alpha_d)*phat;                             
        lambda_dav=alpha_dt.*lambda_dav+(1-alpha_dt).*Ya2;                    % 引导判决法 derive the estimated noise:  lambda_dav
        
        GH1_y=ones(M21,1);
        idx=find(v>5);
        GH1_y(idx)=eta(idx)./(1+eta(idx));
        idx=find(v<=5 & v>0);
        GH1_y(idx)=eta(idx)./(1+eta(idx)).*exp(0.5*expint(v(idx)));
        eta_2term_y=GH1_y.^2.*gamma;
        
        % ref noise
        idx=find(gamma_mint_u>1 & gamma_mint_u<delta_yt & zetat_u<delta_s);
        qhat_u(idx)=(delta_yt-gamma_mint_u(idx))/(delta_yt-1);
        phat_u(idx)=1./(1+qhat_u(idx)./(1-qhat_u(idx)).*(1+eta_u(idx)).*exp(-v_u(idx)));
        phat_u(gamma_mint_u>=delta_yt | zetat_u>=delta_s)=1;
        alpha_dt_u=alpha_d+(1-alpha_d)*phat_u;
        lambda_dav_u=alpha_dt_u.*lambda_dav_u+(1-alpha_dt_u).*Ua2;         
    
        GH1_u=ones(M21,1);
        idx=find(v_u>5);
        GH1_u(idx)=eta_u(idx)./(1+eta_u(idx));
        idx=find(v_u<=5 & v_u>0);
        GH1_u(idx)=eta_u(idx)./(1+eta_u(idx)).*exp(0.5*expint(v_u(idx)));        
        eta_2term_u=GH1_u.^2.*gamma_u; % eta_2term_u=GH1_u.^2.*gamma;
        
           
      
        if l_omlsa<14+l_fnz     % new version omlsa3
            lambda_dav_long=lambda_dav;
            lambda_dav_long_u=lambda_dav_u;
        else
            alpha_dt_long=alpha_d_long+(1-alpha_d_long)*phat;
            lambda_dav_long=alpha_dt_long.*lambda_dav_long+(1-alpha_dt_long).*Ya2;     
           
            alpha_dt_long_u=alpha_d_long+(1-alpha_d_long)*phat_u;
            lambda_dav_long_u=alpha_dt_long_u.*lambda_dav_long_u+(1-alpha_dt_long_u).*Ua2;        
        end
         
        
        % (2.7) update minimum of noise power:  SMact  -> SW  -> Smin 
        %                                       SMactt -> SWt -> Smint 
        l_mod_lswitch=l_mod_lswitch+1;
        if l_mod_lswitch==Vwin
            l_mod_lswitch=0;
           
            if l_omlsa==Vwin-1+l_fnz    % new version omlsa3
                SW=repmat(S,1,Nwin);  SWt=repmat(St,1,Nwin);
                SW_u=repmat(S_u,1,Nwin);  SWt_u=repmat(St_u,1,Nwin);
            else
                SW=[SW(:,2:Nwin) SMact];        SW_u=[SW_u(:,2:Nwin) SMact_u]; % every 15 frame update 1 SMact, Smin is Smallest among 8 shift window
                Smin=min(SW,[],2);              Smin_u=min(SW_u,[],2);         %    SMact -> [ Sw0 Sw1 ... Sw ] -> Smin
                SMact=S;                        SMact_u=S_u;
                SWt=[SWt(:,2:Nwin) SMactt];     SWt_u=[SWt_u(:,2:Nwin) SMactt_u]; 
                Smint=min(SWt,[],2);            Smint_u=min(SWt_u,[],2);
                SMactt=St;                      SMactt_u=St_u;
            end
        end
        % (2.8) derive lambda_d:  lambda_dav->  lambda_d
        %     lambda_d=1.4685*lambda_dav;  % new version
        switch nonstat    % new version
            case   'high'  % new version
                lambda_d=2*lambda_dav;  % new version
                lambda_d_u=2*lambda_dav_u;  % new version
            otherwise  % new version
                lambda_d=1.4685*lambda_dav;  % new version
                lambda_d_u=1.4685*lambda_dav_u;  % new version
        end  % new version
        
        
        %% end of 1) IMCRA
        
        %% 3) add dual channel q_dual estimation routine
         LAMBDA_0 = 1.54;
         LAMBDA_Y =   S ./ max(lambda_d, 1e-10);
         LAMBDA_U = S_u ./ max(lambda_d_u, 1e-10);
          k1_psi = 9; k2_psi = 88;  PSI_0 = 0.20; gamma_0 = 4.6;
         PSI = zeros(M21,1);  OMEGA = ones(M21,1);  OMEGA_HI = 3; OMEGA_LO = 1;
         
%         idx =  LAMBDA_Y > LAMBDA_0 & LAMBDA_U < LAMBDA_0;
%         PSI(idx) = 1;         
%         idx_1 = LAMBDA_Y > LAMBDA_0 & LAMBDA_U > LAMBDA_0;
%         OMEGA(idx_1) = (S(idx_1) - lambda_d(idx_1)) ./ max((S_u(idx_1) - lambda_d_u(idx_1)), 1e-10) ;         
%         idx_1   = OMEGA   >  OMEGA_HI;         
%         PSI(idx_1) = 1;         
%         idx_0    = OMEGA   <=  OMEGA_LO;         
%         PSI(idx_0) = 0;                  
%         idx_1_2  =  OMEGA  <= OMEGA_HI  & OMEGA >OMEGA_LO;
%         PSI(idx_1_2) = (OMEGA(idx_1_2) - OMEGA_LO) ./ max((OMEGA_HI - OMEGA_LO),1e-10);

         % PSI set to be 0 default
         for idx = 1:M21
             if  (LAMBDA_Y(idx) > LAMBDA_0 ) && (LAMBDA_U(idx) <= LAMBDA_0 )
                 PSI(idx) = 1;               
             else
                 if (LAMBDA_Y(idx) > LAMBDA_0 ) && (LAMBDA_U(idx) > LAMBDA_0 )
                     
                     OMEGA(idx) =    (S(idx) - lambda_d(idx)) / max((S_u(idx) - lambda_d_u(idx)), 1e-10) ;                     
                     
                     if   OMEGA(idx) > OMEGA_HI                         
                          PSI(idx) = 1;
                     else
                         if (OMEGA(idx) > OMEGA_LO) && (OMEGA(idx) <= OMEGA_HI)
                             PSI(idx) = (OMEGA(idx) - OMEGA_LO ) / (OMEGA_HI- OMEGA_LO);
                         end                         
                     end                                          
                 end                 
             end            
         end

          
           
         PSI_global = mean(PSI(k1_psi:k2_psi));
         
         q_hat = ones(M21,1); gamma_s = ones(M21,1);
 
         gamma_s  = Ya2  ./  max( lambda_d ,1e-10);         
         gamma_s = min(gamma_0,gamma_s); % cannot bigger than gamma_0         
         idx = gamma_s >1 & PSI_global > PSI_0;      
      
         q_hat(idx) = max((gamma_0 - gamma_s(idx)) /(gamma_0 - 1), 1- PSI(idx)  );
          
        %%
        
        %  SPP = 0.988 * SPP + (1 - 0.988) * (PSI_global > PSI_0);
        % SPP = 0.56 * SPP + (1 - 0.56) * (PSI_global > PSI_0);
           SPP = 0.56 * SPP + (1 - 0.56) * (PSI > PSI_0);
        
        q =min(q_hat,q_max);%use qhat instead of q
         
        lframe = lframe +1;
        
      %  fprintf("PSI_global   %f  k_omlsa %d\n" ,  PSI_global, k_omlsa);
         
 
         gamma=Ya2./max(lambda_d_hat,1e-10);
       %   gamma=Ya2./max(lambda_d,1e-10);
        eta=alpha_eta*eta_2term+(1-alpha_eta)*max(gamma-1,0);
        eta=max(eta,eta_min);
        v=gamma.*eta./(1+eta);
        PH1=zeros(M21,1);
        idx=find(q<0.9);
        PH1(idx)=1./(1+q(idx)./(1-q(idx)).*(1+eta(idx)).*exp(-v(idx)));
        
        
 
        %lambda_d_hat = lambda_d_hat*1.468;
        
        % 7. Spectral Gain
        GH1=ones(M21,1);
        idx=find(v>5);
        GH1(idx)=eta(idx)./(1+eta(idx));
        idx=find(v<=5 & v>0);
        GH1(idx)=eta(idx)./(1+eta(idx)).*exp(0.5*expint(v(idx)));
       
        %     Smint_global=[Smint [Smint(2:M21);Smint(M21)] [Smint(3:M21);Smint(M21-1:M21)] [Smint(4:M21);Smint(M21-2:M21)] [Smint(1);Smint(1:M21-1)] [Smint(1:2);Smint(1:M21-2)] [Smint(1:3);Smint(1:M21-3)]];   % new version
        %     Smint_global=min(Smint_global,[],2);   % new version
        %     lambda_d_global=1.5*Bmin*Smint_global;   % new version
        %    Sy=0.8*Sy+0.2*Ya2;    % new version
        %     GH0=G_f*(lambda_d_global./Sy).^0.5;    % new version
        if tone_flag   % new version
            lambda_d_global=lambda_d_hat;   % new version
            lambda_d_global(4:M21-3)=min([lambda_d_global(4:M21-3),lambda_d_global(1:M21-6),lambda_d_global(7:M21)],[],2);   % new version
            Sy=0.8*Sy+0.2*Ya2;    % new version
            %             GH0=G_f*(lambda_d_global./Sy).^0.5;   % new version       % new version omlsa3
            GH0=G_f*(lambda_d_global./(Sy+1e-10)).^0.5;   % new version omlsa3
        else   % new version
            GH0=G_f;   %#ok<UNRCH> % new version
        end   % new version
        alpha_d_hat = alpha_d  + (1-alpha_d)*PH1;        
        lambda_d_hat = alpha_d_hat .* lambda_d_hat + (1- alpha_d_hat) .* Ya2;
 
        G=GH1.^PH1.*GH0.^(1-PH1);
        
       % G = ones(length(G),1);
        
        eta_2term=GH1.^2.*gamma;
 
        X=[zeros(3,1); G(4:M21-1).*Y(4:M21-1); 0];
  
        X(M21+1:M)=conj(X(M21-1:-1:2)); %extend the anti-symmetric range of the spectum
        x=Cwin^2*win.*real(ifft(X));
        out_omlsa=out_omlsa+x;
    else        
        if ~fnz_flag        
            l_fnz=l_fnz+1;       
        end        
    end       
     
    out_om = out_omlsa(1:Mno);
    
    out_omlsa=[out_omlsa(Mno+1:M); zeros(Mno,1)];   % update output frame
      
    l_omlsa = l_omlsa +1;
   
 

function y = lnshift(x,t)
% lnshift -- t circular left shift of 1-d signal
%  Usage
%    y = lnshift(x,t)
%  Inputs
%    x   1-d signal
%  Outputs
%    y   1-d signal
%        y(i) = x(i+t) for i+t < n
%	 		y(i) = x(i+t-n) else
%
% Copyright (c) 2000. Prof Israel Cohen.
% All rights reserved. 
% ***************************************************************@
szX=size(x);
if szX(1)>1
    n=szX(1);
    y=[x((1+t):n); x(1:t)];
else
    n=szX(2);
    y=[x((1+t):n) x(1:t)];
end

 


