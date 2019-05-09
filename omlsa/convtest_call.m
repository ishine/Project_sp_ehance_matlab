fin0 = '..\..\voice\test6';
 

[x1,fs1]= audioread([fin0 '.wav']); % main mic

len = length(x1);
f_len = 128;


l_s = 1;
l_e = f_len

b = [0.25, 0.5,0.25];

c = zeros(l_e+3, 1);

while(l_e<len)

x1_buf = x1(l_s :l_e);

c = convtest(x1_buf, b);
l_s = l_s+ f_len;
l_e =  l_e + f_len;

end





