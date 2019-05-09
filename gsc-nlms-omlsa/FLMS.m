function[error,filter_coeff ]= FLMS(alpha,filter_length,input_signal,desired_signal,lambda,estimated_power,imp_response) 
% FASTLMS 
% Call: 
% [error,filter_coeff,WEVN] = FLMS(alpha,filter_length,input_signal,desired_signal,lambda,estimated_power,imp_response); 
% 
% ******* Input arguments ******* 
% alpha = step size, dim 1x1 
% filter length, dim 1x1 
% input signal, dim Nx1 
% desired signal, dim Nx1 
% lambda = forgetting factor, dim 1x1 
% estimated power, dim 2Mx1 
% 
% ******* Output arguments ******* 
% estimation error, dim Nx1 
% filter_coefficient, dim Mx1 
% 
 
% initialization 
FILTER_COEFF = zeros(2*filter_length,1); 
input_length = length(input_signal); 
 
%length of signal should be integer multiple of filter_length 
block_length = floor(input_length/filter_length)*filter_length; 
 
%truncate signal 
input_signal = input_signal(1:block_length); 
desired_signal = desired_signal(1:block_length); 
 
% make sure that input signal and desired signal are column vectors 
input_signal = input_signal(:); 
desired_signal = desired_signal(:); 
 
error = desired_signal; 
 
% no.of blocks 
Blocks = block_length/filter_length; 
 
% loop, FLMS 
for k=1:Blocks-1 
 
    % block k-1, k; transformed input signal U(k) 
    INPUT_SIGNAL = fft([input_signal((k-1)*filter_length+1:(k+1)*filter_length)],2*filter_length); 
 
    % block k, output signal y(k), last M elements 
    output = ifft(INPUT_SIGNAL.*FILTER_COEFF); 
    output = output(filter_length+1:2*filter_length,1); 
 
    % block k; desired signal 
    desired_vec = desired_signal(k*filter_length+1:(k+1)*filter_length); 
     
    % block k, error signal 
    error(k*filter_length+1:(k+1)*filter_length,1) = desired_vec-output; 
 
    % transformation of estimation error 
    ERROR_VEC = fft([zeros(filter_length,1);error(k*filter_length+1:(k+1)*filter_length)],2*filter_length); 
     
    % estimated power 
    estimated_power=lambda*estimated_power+(1-lambda)*abs(INPUT_SIGNAL).^2; 
 
    % block k, inverse of power 
    DESIRED_VEC = 1./(1e-6+estimated_power); 
 
    % estimated gradient 
    phivec = ifft(DESIRED_VEC.*conj(INPUT_SIGNAL).*ERROR_VEC,2*filter_length); 
    phivec = phivec(1:filter_length); 
 
    % update of weights 
    FILTER_COEFF = FILTER_COEFF+alpha*fft([phivec;zeros(filter_length,1)],2*filter_length); 
 
    % The error vector should have only real values. 
    error = real(error(:)); 
 
    % transform of final weights to time domain. 
    % make sure that filter coefficient is real-valued 
    filter_coeff = ifft(FILTER_COEFF); 
    filter_coeff = real(filter_coeff(1:length(FILTER_COEFF)/2)); 
     
%     if filter_length < 256 
%         WEVN(k) = 10*log10(sum((imp_response-[filter_coeff;zeros(length(imp_response)-filter_length,1)]).^2,1)/sum(imp_response.^2,1)); 
%     elseif filter_length == 256 
%         WEVN(k) = 10*log10(sum((imp_response-filter_coeff).^2,1)/sum(imp_response.^2,1)); 
%     elseif filter_length > 256  
%         WEVN(k) = 10*log10(sum(([imp_response;zeros(filter_length-length(imp_response),1)]-filter_coeff).^2,1)/sum(imp_response.^2,1)); 
%     end 
end 