
% process one frame flms
function [errorout] = flms_frame( input_signal, desired_signal , filter_length ,alpha , lambda  ,initflag )

    global FILTER_COEFF estimated_power input_signal_buff  k error;
    global SPP;
    
    if(initflag==1)
        FILTER_COEFF = zeros(filter_length*2,1);
        estimated_power = zeros(filter_length*2 ,1);
        input_signal_buff = input_signal;
        k = 1;
        
        SPP = zeros(filter_length,1);
        
    else         
        input_signal_buff = [ input_signal_buff(filter_length+1:filter_length*2) ; input_signal ];        
        k = k +1;
    end
    
    

   % block k-1, k; transformed input signal U(k) 
    INPUT_SIGNAL = fft(input_signal_buff,2*filter_length); 
 
    % block k, output signal y(k), last M elements 
    output = ifft(INPUT_SIGNAL.*FILTER_COEFF); 
    output = output(filter_length+1:2*filter_length,1); 
    
    
    yout(k*filter_length+1:(k+1)*filter_length) = output;
 
    % block k; desired signal 
    desired_vec = desired_signal; 
     
    % block k, error signal   (time domain error)
    error(k*filter_length+1:(k+1)*filter_length,1) = desired_vec-output; 
 
    % transformation of estimation error 
    ERROR_VEC = fft([zeros(filter_length,1);error(k*filter_length+1:(k+1)*filter_length)],2*filter_length); 
     
    % estimated power 
    estimated_power=lambda*estimated_power+(1-lambda)*abs(INPUT_SIGNAL).^2; 
 
    % block k, inverse of power 
    DESIRED_VEC = 1./(1+estimated_power); 
 
    % estimated gradient 
    phivec = ifft(DESIRED_VEC.*conj(INPUT_SIGNAL).*ERROR_VEC,2*filter_length); 
    phivec = phivec(1:filter_length); 
 
    % update of weights 
    FILTER_COEFF = FILTER_COEFF+ alpha*fft([phivec;zeros(filter_length,1)],2*filter_length); 

    error = real(error(:)); 
    
    errorout = error(k*filter_length+1:(k+1)*filter_length);
    
end




