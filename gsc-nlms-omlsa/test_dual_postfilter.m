fin = './voice/t194';
fout = [fin,'out'];


%[in, fs] = audioread(['./voice/t193.wav']);
%inL = in(:,1);
%inR = in(:,2);

%inL =  inL - inR;
%audiowrite(['./voice/t193_del.wav'],inL,fs);


[in, fs] = audioread([fin,'.wav']);

inL = in(:,1);
inR = in(:,2);

lenS = length(inL);
out = zeros(lenS,1);
Nlen = 512;
Nlen41 = Nlen/4;

Blk = fix(lenS / Nlen41);

frameHead = 1; frameEnd = Nlen;
l = 1;

global lambda_dav_u  lambda_d_hat St_u  q_hat  S_u  Smin_u  enoise  lambda_d   Ya2 ;


lambda_dav_t = zeros(lenS*3,2);

while(frameEnd < lenS) 
    
    fL = inL(frameHead:frameEnd);
    fR = inR(frameHead:frameEnd);
    
  %  fL1 = 0.5*( fL + fR);
  %  fR1 = fL - fR;
   
  %  [out_om] = tf_gsc_dual_postfilter_loizou(fL, fR , l);
    [out_om] = tf_gsc_dual_postfilter(fL, fR , l);

    out((l-1)*Nlen41+1 : l*Nlen41) = out_om;
    
     lenn = 257 ;
    
    cbnoise = [Ya2, lambda_d];
     
     lambda_dav_t((l-1)*lenn+1 : l*lenn ,: ) =  cbnoise;
    
    l = l +1;
    frameHead = frameHead + Nlen41;
    frameEnd  = frameEnd  + Nlen41;
end


  audiowrite([fout '_Om.wav'],out,fs);
  audiowrite([fout '_ncmp.wav'],lambda_dav_t,fs);
