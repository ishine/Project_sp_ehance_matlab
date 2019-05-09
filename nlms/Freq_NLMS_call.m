% freq domain NLMS algorithm
% 频域NLMS, 验证能起到降噪效果, 声音正常， 而且降噪效果要好于时域NLMS  NFFT = 1024时效果较好，低于该值效果差。


function Out = Freq_NLMS_call(Sig,Ref, l, ch,NFFT ) 
   
      OVERLAP = NFFT/2;
global Pest ; 
global G;
global Out;
global PestH ; 
global GH;

global His_Out1; 
global His_Out2;
 Out= zeros(NFFT,1);
  omega = 0.36; 
  alpha = 0.075 ;
     if(l==1)
           Pest = zeros(NFFT,1);         
           G = zeros(NFFT,1);         
           PestH = zeros(NFFT,1);         
           GH = zeros(NFFT,1);                        
           His_Out1 = zeros(OVERLAP,1);
           His_Out2 = zeros(OVERLAP,1);        
     end
 
     window = hanning(NFFT);       
     Frame1 =  Sig;% Y_L(InitFrame:EndFrame); % Y_Up(InitFrame:EndFrame);
     Frame2 =  Ref; %Y_R(InitFrame:EndFrame); %Y_Down(InitFrame:EndFrame);    
     X1 = fft(window .* Frame1);
     X2 = fft(window .* Frame2);     
      
     if(ch==1)    
        X_lms =  G.* X2;          
        En = X1 - X_lms;    
        Pest = omega * Pest + (1-omega) * abs(X2).^2;          
        mu = alpha ./ Pest;
        G = G + mu .* En .* conj( X2)  ;
        Out  = ifft(En);         
        Out(1:OVERLAP) = Out(1:OVERLAP) + His_Out1 ;
        His_Out1       = Out(OVERLAP+1:NFFT) ;
     else
        GH = zeros(NFFT,1);   
        X_lms =  GH.* X2;          
        En = X1 - X_lms;    
        PestH = omega * PestH + (1-omega) * abs(X2).^2;          
        mu = alpha ./ PestH;
        GH = GH + mu .* En .* conj( X2)  ;
        Out  = ifft(En);         
        Out(1:OVERLAP) = Out(1:OVERLAP) + His_Out2 ;
        His_Out2   = Out(OVERLAP+1:NFFT) ;  
     end
end
       
     
    
 

 
 

