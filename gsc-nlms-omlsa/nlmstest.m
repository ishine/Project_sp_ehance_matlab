

 
  d = randn(128,1) ;
  u = randn(384,1);  
  pest = zeros(256,1);
  A_st = zeros(256,1);
  d1=fix(d);
  u1=fix(u);
  pest1 = fix(pest);
  A_st1 = fix(A_st);
  
%e = nlms(d1,u1,  pest1, A_st1 ,128  );
%e = nlms_fixpt(d1,u1,  pest1, A_st1 ,128  );