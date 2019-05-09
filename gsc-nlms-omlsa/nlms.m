function [e]=   nlms(d, u, Pest,A_st, N ) 

     Nflt = 256; omega = 0.898;  alpha = 0.072;
      e = zeros(N,1);
     
     for k = 1: N

     Y_Frame_Block = u(k:Nflt+k-1);      
     Y_Frame_up = d(k);   
     yB =  Y_Frame_Block' * Y_Frame_Block;  
    %% apply the coeff             
     Y_Down = A_st'*Y_Frame_Block;        
     yout = Y_Frame_up - Y_Down ;        
     e(k) = yout;  
    %% update the coeff  

     Pest = omega * Pest + (1-omega) * yB;          
     mu = alpha ./( Pest + 1e-2);     
     A_st = A_st +  mu*yout.* (Y_Frame_Block) ;  
     
     end
end


  
 
