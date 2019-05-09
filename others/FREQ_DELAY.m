function FREQ_DELAY()

fin =  '..\voice\test';
fin2 = '..\voice\test2';
fout = [fin '_FREQ_DELAY'];
[y,fs]=audioread([fin,'.wav']);  % read Y, beamform output
[u,fs]=audioread([fin2,'.wav']); % read U,  reference noise 

len = length(y);
lenb = length(u);
len =max(len,lenb);
Out_T = zeros(len,1);
NFFT = 1024;
OVERLAP = NFFT/2;
win=hamming(NFFT);
His_Out = zeros(OVERLAP,1);
FRAME_NUM = len / NFFT;
l=1;

s_frame = 1;
e_frame = NFFT;
G = 1;
PI = 3.14159269;

TimeDelay = -8;% number of delay

bin = (0:NFFT -1)';
fbin = bin* 2* PI*TimeDelay / NFFT; % w

c = 343;
DesAng = 180;
 
Tau_c    =exp(-i*fbin) ;

while(e_frame< len)
    
    y_frame = y(s_frame:e_frame);
    u_frame = u(s_frame:e_frame);
    Y = fft(win.* y_frame);
    U = fft(win.* u_frame);
    
    
    AY = Y./ abs(Y);
    
   %  Y = Y + U;
    
    Yout = Y.*Tau_c ;
     
    
    ang  = angle(AY) ;
    ang2 =  angle(U) ;
    
    delang = cos(ang - ang2);
    
    Out = real(ifft(Yout));
    if(l==1)
       Out_T((l-1)*OVERLAP+1 : l*OVERLAP) = 0;  
    else
     Out_T((l-1)*OVERLAP+1 : l*OVERLAP) = Out(1:OVERLAP) + His_Out;
    end
     His_Out = Out(OVERLAP+1 : NFFT);
    l = l+1;
   
    s_frame = s_frame+ OVERLAP;
    e_frame = e_frame+ OVERLAP;
 
end

Out_T_fix = [Out_T,y];

 audiowrite([fout '.wav'],Out_T,fs);

end