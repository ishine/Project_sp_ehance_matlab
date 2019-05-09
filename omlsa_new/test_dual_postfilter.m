fin = './t192';
fout = [fin,'out'];

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
while(frameEnd < lenS) 
    
    fL = inL(frameHead:frameEnd);
    fR = inR(frameHead:frameEnd);
    
  %  fL1 = 0.5*( fL + fR);
  %  fR1 = fL - fR;
    
    [out_om] = tf_gsc_dual_postfilter(fL, fR , l);
    out((l-1)*Nlen41+1 : l*Nlen41) = out_om;
    
    l = l +1;
    frameHead = frameHead + Nlen41;
    frameEnd  = frameEnd  + Nlen41;
end


audiowrite([fout '_Om.wav'],out,fs);
 
