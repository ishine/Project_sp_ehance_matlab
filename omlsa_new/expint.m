function y = expint(x)
%EXPINT Exponential integral function.
%   Y = EXPINT(X) is the exponential integral function for each
%   element of X.  The exponential integral is defined as:
%
%   EXPINT(x) = integral from x to Inf of (exp(-t)/t) dt, for x > 0.
%   
%   By analytic continuation, EXPINT is a scalar-valued function in
%   the complex plane cut along the negative real axis.
%
%   Another common definition of the exponential integral function is
%   the Cauchy principal value integral from -Inf to X of (exp(t)/t)
%   dt, for positive X.  This is denoted as Ei(x). The relationships
%   between EXPINT(x) and Ei(x) are as follows:
%
%       EXPINT(-x+i*0) = -Ei(x) - i*pi, for real x > 0
%       Ei(x) = REAL(-EXPINT(-x)), for real x > 0
%
%   Class support for input X:
%      float: double, single

%   Copyright 1984-2011 The MathWorks, Inc. 

%       For elements of X in [-38,2], EXPINT uses a series expansion
%       representation (equation 5.1.11 from Abramowitz & Stegun).
%       For all other elements of X, EXPINT uses a continued fraction
%       representation (equation 5.1.22 in A&S).
%       References:
%         [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%         Functions", Dover Publications, 1965, Ch. 5.

siz = size(x);
y = zeros(numel(x),1,superiorfloat(x));

% make input a column vector
x = x(:);

% figure out which algorithm to use by evaluating interpolating polynomial 
% at real(z)

p = [-3.602693626336023e-09 -4.819538452140960e-07 -2.569498322115933e-05 ...
     -6.973790859534190e-04 -1.019573529845792e-02 -7.811863559248197e-02 ...
     -3.012432892762715e-01 -7.773807325735529e-01  8.267661952366478e+00];
polyv = polyval(p,real(x));

% series expansion

k = find( abs(imag(x)) <= polyv );
if ~isempty(k)

   %initialization
   egamma=0.57721566490153286061;
   xk = x(k);
   yk = -egamma - log(xk);
   j = 1;
   pterm = xk;
   term = xk;

   while any(abs(term) > (eps(yk)))
      yk = yk + term;
      j = j + 1;
      pterm = -xk.*pterm/j;
      term = pterm/j;
   end % end of the while loop
 
   y(k) = yk;
end

% continued fraction

k = find( abs(imag(x)) > polyv );
if ~isempty(k)
   %   note: am1, bm1 corresponds to A(j-1), B(j-1) of recursion formulae
   %         am2, bm2 corresponds to A(j-2), B(j-2) of recursion formulae
   %         a,b      corresponds to A(j), B(j) of recursion formulae

   n = 1; % we're calculating E1(x)

   % initialization
   xk = x(k);
   am2 = zeros(size(xk));
   bm2 = ones(size(xk));
   am1 = ones(size(xk));
   bm1 = xk;
   f = am1 ./ bm1;
   oldf = Inf(size(xk));
   j = 2;

   while any(abs(f-oldf) > (100*eps(class(f)).*abs(f)))
       % calculate the coefficients of the recursion formulas for j even
       alpha = n-1+(j/2); % note: beta= 1
   
       %calculate A(j), B(j), and f(j)
       a = am1 + alpha * am2;
       b = bm1 + alpha * bm2;
   
       % save new normalized variables for next pass through the loop
       %  note: normalization to avoid overflow or underflow
       am2 = am1 ./ b;
       bm2 = bm1 ./ b;
       am1 = a ./ b;
       bm1 = 1;
   
       f = am1;
       j = j+1;
   
       % calculate the coefficients for j odd
       alpha = (j-1)/2;
       beta = xk;
       a = beta .* am1 + alpha * am2;
       b = beta .* bm1 + alpha * bm2;
       am2 = am1 ./ b;
       bm2 = bm1 ./ b;
       am1 = a ./ b;
       bm1 = 1;
       oldf = f;
       f = am1;
       j = j+1;
   
   end  % end of the while loop
    
   y(k)= exp(-xk) .* f - 1i*pi*((real(xk)<0)&(imag(xk)==0)); 
end

y = reshape(y,siz);
