function [e_noise] = noise_estimate(Ya2,flag)
%UNTITLED3 此处显示有关此函数的摘要
%  A noise-estimation algorithm for highly non-stationary environments 2006
 global P_n  Pmin_n  Sr_n p_hat_n alpha_n_s D_n  P_n_last  detla   ;  
 
 M21 = 257;
 thta = 0.85;
 
 alpha_n_p = 0.2; alpha_n_d = 0.85;
 gamma = 0.998;
 beta = 0.8;
 if flag==1     
     P_n = Ya2; Pmin_n = zeros(M21,1);     
     Sr_n = zeros(M21,1);  p_hat_n = zeros(M21,1);
     alpha_n_s = zeros(M21,1);  D_n= Ya2; 
     
     LF = fix(M21 / 8000 * 1000);
     MF = fix(M21 / 8000 * 3000);
     
     detla(1:LF)    = 2;
     detla(LF+1:MF) = 2;
     detla(MF+1:M21) = 5;
 end
 
 % update Power of signal
 P_n_last = P_n;
 P_n = thta * P_n + (1-thta) *Ya2; 

 
  % update min power of signal
  Pmin_n_bak = Pmin_n;
  Pmin_n = P_n; 
  idx = Pmin_n < P_n;  
  Pmin_n(idx) = gamma * Pmin_n_bak(idx) +  (1-gamma)/(1-beta) * (P_n(idx)  - beta * P_n_last(idx) );
  
  % cal probability   
  Sr_n = P_n ./ max(Pmin_n, 1e-10);

  I = zeros(M21,1);  
  idx = Sr_n > detla;
  I(idx) = 1;
  p_hat_n = alpha_n_p * p_hat_n + (1-alpha_n_p) * I;
  
  alpha_n_s = alpha_n_d + (1-alpha_n_d)* p_hat_n;
   
  D_n = alpha_n_s .* D_n + (1-alpha_n_s) .* Ya2; 

  e_noise = D_n;
end

