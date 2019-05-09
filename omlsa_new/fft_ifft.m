function [x]=fft_ifft(y)
 


 X=fft( y);
     
     
 x= real(ifft(X));