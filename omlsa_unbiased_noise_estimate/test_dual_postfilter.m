fin = '../../voice/16';
fout = [fin,'out'];


%[in, fs] = audioread(['./voice/t193.wav']);
%inL = in(:,1);
%inR = in(:,2);

%inL =  inL - inR;
%audiowrite(['./voice/t193_del.wav'],inL,fs);

arg1 = 113; arg2 = 0.25;
[in, fs] = audioread([fin,'.wav']);

inL = in(:,2);
inR = in(:,1);

 
lenS = length(inL);
out = zeros(lenS,1);
Nlen = 512;
Nlen41 = Nlen/4;

Blk = fix(lenS / Nlen41);

frameHead = 1; frameEnd = Nlen;
l = 1;

global lambda_dav_u  lambda_d_hat St_u  q_hat  S_u  Smin_u  enoise  lambda_d   Ya2 ;


lambda_dav_t = zeros(lenS*3,2);
out = zeros(lenS,1);
while(frameEnd < lenS) 
    
    fL = inL(frameHead:frameEnd);
    fR = inR(frameHead:frameEnd);
    
 
    [out_om] = tf_gsc_dual_postfilter_ad( fL ,fR, l,  arg1, arg2);
%function [out_om]=tf_gsc_dual_postfilter_ad(y_frame_in , u_frame_in , framecnt, Fs, arg1, arg2)

    out((l-1)*Nlen41+1 : l*Nlen41) = out_om;
    
     lenn = 257 ;
    
     cbnoise = [Ya2, lambda_d];
     
    
    l = l +1;
    frameHead = frameHead + Nlen41;
    frameEnd  = frameEnd  + Nlen41;
end


  audiowrite([fout '_Om_ad.wav'],out,fs);
   
  
  frameHead = 1; frameEnd = Nlen; l = 1;
  
  out = zeros(lenS,1);
  
  while(frameEnd < lenS) 
    
    fL = inL(frameHead:frameEnd);
    fR = inR(frameHead:frameEnd);
    
    [out_om] = tf_gsc_dual_postfilter2(fL, fR , l,arg1, arg2 );

    out((l-1)*Nlen41+1 : l*Nlen41) = out_om;
    
     lenn = 257 ;
    
    cbnoise = [Ya2, lambda_d];
     
   
    l = l +1;
    frameHead = frameHead + Nlen41;
    frameEnd  = frameEnd  + Nlen41;

  end


  audiowrite([fout '_Om.wav'],out,fs);
   
  
  
