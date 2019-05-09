% Input
clear all;

fin0 = '..\voice\test4';
 
  
fout = [fin0 'bin'];

delete('../voice/psd.bin');
 fid= fopen('../voice/psd.bin','a+');
  
%fout1 = [fin0 'B'];

[Y_L,fs1]= audioread([fin0 '.wav']); % main mic

len = length(Y_L);

NFFT = 512;
window =  hamming(NFFT);


tarr = [1,2,3,4,5,6,7];


ptar = prod(tarr(1:4));


OVERLAP = NFFT/2;
l=1;

s_frame = 1;
e_frame = NFFT;

M = round(len/OVERLAP);

PSD2500 = zeros(M, 1);
PSD2000 = zeros(M, 1);

SFM = zeros(16,1);

SFM3 = zeros(M, 1);
SFM4 = zeros(M, 1);

W_a1 = zeros(M, 1);
W_a2 = zeros(M, 1);
value = zeros(NFFT, 1);
while(e_frame<len)

  y_frame =   Y_L(s_frame:e_frame);
    
  Y_F =  fft(window.*y_frame);
  
  Y_F_2 = abs(Y_F).^2;
   
   PSD = 20*log10(Y_F_2);
   
   fprintf( fid ,'PSD2500= %f   PSD2000= %f\n', PSD(76 ),PSD(64 )  );
   
   PSD2500(l) = PSD(15 );
   PSD2000(l) = PSD(64 );
   
   
   s_inx = 1;
   e_inx = 16;
   
   for k= 1:16
        SFM(k) =(( prod(Y_F_2(s_inx:e_inx))^(1/16))) / (sum(Y_F_2(s_inx:e_inx))/16  );
        
        s_inx = s_inx +16;
        e_inx = e_inx +16;
   end
   
   
   
   % w cal 
   %value = new if new < old
   %value = (new-old )*0.01 + old  if new >=old    
   % old = new 
   newV = Y_F_2;
   
   if(l>=2)
      inx = newV < oldV;
      value(inx) = newV(inx);      
      inx = newV > oldV;      
      value(inx) =  (newV(inx)- oldV(inx))*0.01 + oldV(inx);    
      
      W_a1(l) = value(6);
      W_a2(l) = value(15);
   end
   
   
   oldV = newV;
    
     SFM3(l) =  SFM(1) ;
     SFM4(l) =  SFM(3) ;  
     
     
     
  
   l = l+1;
  s_frame = s_frame+OVERLAP;
  e_frame = e_frame+OVERLAP;
end
figure
n = (1:M);
%plot(n,PSD2500,'R',n, PSD2000, 'G',n, SFM3 ,n,  SFM4);
plot(n, SFM3 ,'G',n,  SFM4,'R');

hold on

figure
n = (1:M);
 plot(n,PSD2500,'R',n, PSD2000, 'G' );

hold on

figure
n = (1:M);
 plot(n,W_a1,'R',n, W_a2, 'G' );
 % plot( n, W_a2, 'G' );
fclose(fid);
