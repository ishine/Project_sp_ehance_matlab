function [q] = q_estimate(eta , ch ,l, lambda_dav_long, Fs)
global   xi_glo  xi_frame_glo  xi_m_dB_glo;       
 M21 = 257; M = 512;
 
alpha_xi =0.7;	% 3.1) Recursive averaging parameter
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

k_u=round(f_u/Fs*M+1);  % Upper frequency bin for global decision
k_l=round(f_l/Fs*M+1);  % Lower frequency bin for global decision
k_u=min(k_u,M21);
k2_local=round(500/Fs*M+1);
k3_local=round(3500/Fs*M+1);
 
b_xi_local=hanning(2*w_xi_local+1);
b_xi_local=b_xi_local/sum(b_xi_local);  % normalize the window function
b_xi_global=hanning(2*w_xi_global+1);
b_xi_global=b_xi_global/sum(b_xi_global);   % normalize the window function

 
    
       if(l<=4) % init global parameter
           xi_glo       = zeros(M21,2);  
           xi_frame_glo = zeros(1,2);
           xi_m_dB_glo  = zeros(1,2);
       end
        
         xi = xi_glo(:,ch);
         xi_frame = xi_frame_glo(:,ch);
         xi_m_dB  = xi_m_dB_glo(:,ch);
       
        xi=alpha_xi*xi+(1-alpha_xi)*eta;
        xi_local=conv(xi,b_xi_local);
        xi_local=xi_local(w_xi_local+1:M21+w_xi_local);
        xi_global=conv(xi,b_xi_global);
        xi_global=xi_global(w_xi_global+1:M21+w_xi_global);
        dxi_frame=xi_frame;
        xi_frame=mean(xi(k_l:k_u));
        dxi_frame=xi_frame-dxi_frame;
        if xi_local>0, xi_local_dB=10*log10(xi_local);     else xi_local_dB=-100;   end
        if xi_global>0, xi_global_dB=10*log10(xi_global);  else  xi_global_dB=-100; end
        if xi_frame>0, xi_frame_dB=10*log10(xi_frame);     else xi_frame_dB=-100;   end

        P_local=ones(M21,1);
        P_local(xi_local_dB<=xi_ll_dB)=P_min;
        idx=find(xi_local_dB>xi_ll_dB & xi_local_dB<xi_lu_dB);
        P_local(idx)=P_min+(xi_local_dB(idx)-xi_ll_dB)/(xi_lu_dB-xi_ll_dB)*(1-P_min);

        P_global=ones(M21,1);
        P_global(xi_global_dB<=xi_gl_dB)=P_min;
        idx=find(xi_global_dB>xi_gl_dB & xi_global_dB<xi_gu_dB);
        P_global(idx)=P_min+(xi_global_dB(idx)-xi_gl_dB)/(xi_gu_dB-xi_gl_dB)*(1-P_min);

        m_P_local=mean(P_local(3:(k2_local+k3_local-3)));    % average probability of speech presence
        if m_P_local<0.25
            P_local(k2_local:k3_local)=P_min;    % reset P_local (frequency>500Hz) for low probability of speech presence
        end
    
            if (m_P_local<0.5) && (l>120)
                idx=find( lambda_dav_long(8:(M21-8)) > 2.5*(lambda_dav_long(10:(M21-6))+lambda_dav_long(6:(M21-10))) );
                P_local([idx+6;idx+7;idx+8])=P_min;   % remove interfering tonals
            end
                    % new version

        if xi_frame_dB<=xi_fl_dB
            P_frame=P_min;
        elseif dxi_frame>=0
            xi_m_dB=min(max(xi_frame_dB,xi_ml_dB),xi_mu_dB);
            P_frame=1;
        elseif xi_frame_dB>=xi_m_dB+xi_fu_dB
            P_frame=1;
        elseif xi_frame_dB<=xi_m_dB+xi_fl_dB
            P_frame=P_min;
        else
            P_frame=P_min+(xi_frame_dB-xi_m_dB-xi_fl_dB)/(xi_fu_dB-xi_fl_dB)*(1-P_min);
        end

         q=1-P_global.*P_local*P_frame;   % new version
          
         q = min(q,q_max);
  
          xi_glo(:,ch) = xi;
          xi_frame_glo(:,ch)= xi_frame;
          xi_m_dB_glo(:,ch) = xi_m_dB;
 
end

