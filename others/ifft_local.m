function [y]=ifft_local(X)

Y = fft(X);
 
y=real(ifft(Y));



