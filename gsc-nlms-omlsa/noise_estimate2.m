function [e_noise2] = noise_estimate2(Ua2,flag)
%UNTITLED3 此处显示有关此函数的摘要
%  A noise-estimation algorithm for highly non-stationary environments 2006
 global P_n2  Pmin_n2  Sr_n2 p_hat_n2 alpha_n_s2 D_n2  P_n_last2  detla2   ;  
 
 M21 = 257;
 thta = 0.85;
 
 alpha_n_p = 0.2; alpha_n_d = 0.85;
 gamma = 0.998;
 beta = 0.8;
 if flag==1     
     P_n2 = Ua2; Pmin_n2 = zeros(M21,1);     
     Sr_n2 = zeros(M21,1);  p_hat_n2 = zeros(M21,1);
     alpha_n_s2 = zeros(M21,1);  D_n2= Ua2; 
     
     LF = fix(M21 / 8000 * 1000);
     MF = fix(M21 / 8000 * 3000);
     
     detla2(1:LF)    = 2;
     detla2(LF+1:MF) = 2;
     detla2(MF+1:M21) = 5;
 end
 
 % update Power of signal
 P_n_last2 = P_n2;
 P_n = thta * P_n2 + (1-thta) *Ua2; 

 
  % update min power of signal
  Pmin_n_bak = Pmin_n2;
  Pmin_n2 = P_n; 
  idx = Pmin_n2 < P_n;  
  Pmin_n2(idx) = gamma * Pmin_n_bak(idx) +  (1-gamma)/(1-beta) * (P_n(idx)  - beta * P_n_last2(idx) );
  
  % cal probability   
  Sr_n2 = P_n ./ max(Pmin_n2, 1e-10);

  I = zeros(M21,1);  
  idx = Sr_n2 > detla2;
  I(idx) = 1;
  p_hat_n2 = alpha_n_p * p_hat_n2 + (1-alpha_n_p) * I;
  
  alpha_n_s2 = alpha_n_d + (1-alpha_n_d)* p_hat_n2;
   
  D_n2 = alpha_n_s2 .* D_n2 + (1-alpha_n_s2) .* Ua2; 

  e_noise2 = D_n2;
end

