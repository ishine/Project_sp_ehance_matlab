function C= convtest( a,b)
%CONV Convolution and polynomial multiplication.
%   C = CONV(A, B) convolves vectors A and B.  The resulting vector is
%   length MAX([LENGTH(A)+LENGTH(B)-1,LENGTH(A),LENGTH(B)]). If A and B are
%   vectors of polynomial coefficients, convolving them is equivalent to
%   multiplying the two polynomials.
%
%   C = CONV(A, B, SHAPE) returns a subsection of the convolution with size
%   specified by SHAPE:
%     'full'  - (default) returns the full convolution,
%     'same'  - returns the central part of the convolution
%               that is the same size as A.
%     'valid' - returns only those parts of the convolution 
%               that are computed without the zero-padded edges. 
%               LENGTH(C)is MAX(LENGTH(A)-MAX(0,LENGTH(B)-1),0).
%
%   Class support for inputs A,B: 
%      float: double, single
%
%   See also DECONV, CONV2, CONVN, FILTER, XCORR, CONVMTX.
%
%   Note: XCORR and CONVMTX are in the Signal Processing Toolbox.


%   Copyright 1984-2013 The MathWorks, Inc.

 %a = [16384,8192, 4096, 2048, 1024, 0 -1024, -2048, -4096, -8192, -16384];
 %b = [0.25, 0.5,0.25];
 
 lenA = length(a);
 lenB = length(b);
 len = lenA + lenB;
 
 
 C = zeros(len,1);
 for l=1: len
     
     for k = 1 : lenA
     
         if ((  l- k)>0 && (l-k)<lenB)
             C(l) = a(k)*b(l-k);
         end
         
     end
     
     
 end 
 
 
% compute as if both inputs are column vectors
 
end

    
    
    
