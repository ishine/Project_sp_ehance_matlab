function [x]=fft_ifft_fix(y)
 


 X=fft( y);
     
     
 x= real(ifft(X));